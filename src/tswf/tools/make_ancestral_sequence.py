"""Recreate ancestral sequence from VCF that contains an AA INFO tag.

Make ancestral reference sequence OUTFILE given REFERENCE sequence and
a VCF file with the AA INFO tag denoting the ancestral allele.

"""

import logging
import pathlib
import sys

import click
import cyvcf2
from pyfaidx import Fasta
from tqdm import tqdm


@click.command(help=__doc__)
@click.argument("reference", type=click.Path(exists=True))
@click.argument("vcf", type=click.Path(exists=True))
@click.argument("outfile", type=click.Path())
@click.option("--chromosome", type=str, help="chromosome")
def main(
    reference: pathlib.Path,
    vcf: pathlib.Path,
    outfile: pathlib.Path,
    chromosome: None | str = None,
) -> None:
    """Generate ancestral sequence from VCF."""
    cyvcf = cyvcf2.VCF(vcf)

    if chromosome is None:
        chromosome = cyvcf.seqnames[0]

    ref = Fasta(reference)
    try:
        seq = ref[chromosome]
    except KeyError as e:
        print(e)
        raise

    n_sites = 0
    n_change = 0
    # Loop vcf and update reference at positions
    for _, variant in tqdm(enumerate(cyvcf)):
        if "AA" not in dict(variant.INFO).keys():
            logging.error(
                "No AA tag in INFO field; please add AA "
                "tag that indicates ancestral state"
            )
            sys.exit(1)
        n_sites = n_sites + 1
        print(variant.INFO["AA"], variant.REF, variant.POS)
        if variant.INFO["AA"] != variant.REF:
            n_change = n_change + 1
            print(
                f"Changing REF {variant.REF} ({variant.REF[variant.POS - 1]})"
                "to {variant.INFO['AA']}"
            )
            seq[variant.POS - 1] = variant.INFO["AA"]
    seq.write(outfile)
