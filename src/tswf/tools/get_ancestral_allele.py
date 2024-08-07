"""Infer ancestral alleles in vcf based on majority vote in OUTGROUPs.

Infer ancestral alleles in VCF file based on majority vote in
OUTGROUPs. The OUTGROUPs are sample ids in the VCF.

"""

import logging
import pathlib
import re
import sys
from typing import TextIO

import click
import cyvcf2
from tqdm import tqdm


@click.command(help=__doc__)
@click.argument("vcf", type=click.Path(exists=True))
@click.argument("outgroup", nargs=-1, required=True)
@click.option("--outfile", type=str)
@click.option("--ploidy", type=int, default=2, help="Set ploidy")
@click.option(
    "--min-alleles",
    "-n",
    type=click.IntRange(
        1,
    ),
    default=1,
    help="number of alleles that need to agree on ancestral state",
)
@click.option("--logfile", type=str, help="logfile name")
def main(  # noqa: C901
    vcf: pathlib.Path,
    outgroup: tuple[str],
    outfile: None | str,
    ploidy: int,
    min_alleles: int,
    logfile: None | str,
) -> None:
    """Infer ancestral alleles in VCF based on majority vote in OUTGROUPs."""
    if outfile is None:
        outfile = re.sub(r"(.vcf.gz)$", "_AA_vote\\1", str(vcf))
    n_vote = min_alleles
    if n_vote > len(outgroup) * ploidy:
        n_vote = len(outgroup) * ploidy
        logging.warning(
            (
                "min_alleles cannot be larger than total number of alleles in outgroups:"
                " setting to %i"
            ),
            n_vote,
        )

    cyvcf = cyvcf2.VCF(vcf)
    # Make sure outgroups in samples
    if not set(outgroup) <= set(cyvcf.samples):
        logging.error("not all outgroup names found in vcf sample header")
        sys.exit(1)

    n_pass = 0
    n_skip = 0
    n_change = 0
    outgroup_indices = [
        i for i in range(len(cyvcf.samples)) if cyvcf.samples[i] in outgroup
    ]
    outext = re.compile(r"(.vcf|.vcf.gz|.bcf)$").search(outfile)
    if outext is None:
        logging.error("please use file extension .vcf, .vcf.gz, or .bcf")
        sys.exit(1)
    mode = {".vcf": "w", ".vcf.gz": "wz", ".bcf": "wb"}[outext.group(1)]

    cyvcf.add_info_to_header(
        {
            "ID": "AA",
            "Number": 1,
            "Type": "Character",
            "Description": "Ancestral Allele",
        }
    )
    vcfwriter = cyvcf2.cyvcf2.Writer(outfile, cyvcf, mode)
    for variant in tqdm(cyvcf):
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get("AA", variant.REF)
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        if len(ordered_alleles) > 2:
            continue
        outgroup_alleles = []
        anc = None
        for i in outgroup_indices:
            outgroup_alleles.extend([variant.genotypes[i][0], variant.genotypes[i][1]])
        # More derived than vote implies ancestral is actually the ALT allele
        if sum(outgroup_alleles) >= n_vote:
            anc = variant.ALT
            n_change = n_change + 1
        # More 0 alleles than vote implies REF is actually the ancestral; no change
        elif len(outgroup_alleles) - sum(outgroup_alleles) >= n_vote:
            anc = variant.REF
        if isinstance(anc, list):
            anc = anc[0]
        # We have a good allele
        if anc is not None:
            n_pass = n_pass + 1
            try:
                variant.INFO["AA"] = anc
            except AttributeError as e:
                print(e)
                print(variant, anc)
                raise
        vcfwriter.write_record(variant)
    else:
        # We only keep alleles that actually pass the selection
        # criteria (n_vote)
        n_skip = n_skip + 1
    vcfwriter.close()

    fh: TextIO = sys.stdout
    if isinstance(logfile, str):
        fh = open(logfile, "w")

    fh.write("get_ancestral_allel.py summary\n")
    fh.write("------------------------------\n")
    fh.write(f"outgroups: {outgroup}\n")
    fh.write(f"number of concordant alleles required in outgroups (n): {n_vote}\n")
    fh.write(
        f"number of sites passing vote (either at least "
        f"{n_vote} '0' or {n_vote} '1' calls) (n): {n_pass}\n"
    )
    fh.write(f"number of sites changed ref (>={n_vote} 1 calls) (n): {n_change}\n")
    fh.write(f"number of skipped sites (n): {n_skip}\n")

    fh.close()
