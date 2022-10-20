#!/usr/bin/env python3
#
# Simple script to infer ancestral alleles from outgroups based on
# majority vote
#
import argparse
import logging
import re
import sys

import cyvcf2
from tqdm import tqdm

try:
    inputfile = snakemake.input.vcf
    outputfile = snakemake.output.vcf
    outgroup = snakemake.params.outgroup
    options = snakemake.params.options
    logfile = str(snakemake.log)
except NameError:
    inputfile = None
    outputfile = None
    outgroup = None
    options = {"n": 1, "ploidy": 2}
    logfile = sys.stdout

parser = argparse.ArgumentParser(
    description="""
Infer ancestral alleles in vcf based on majority vote in outgroups.
"""
)
parser.add_argument(
    "--outgroup",
    metavar="outgroup",
    type=str,
    help="outgroup name",
    action="append",
    default=outgroup,
)
parser.add_argument(
    "--ploidy",
    metavar="ploidy",
    type=int,
    default=options.get("ploidy", 2),
    help="set the ploidy",
)
parser.add_argument(
    "-n",
    metavar="n",
    type=int,
    default=None,
    help="number of haplotypes that need to agree on ancestral state",
)
parser.add_argument(
    "vcf", type=str, nargs="?", help="Input vcf file", default=inputfile
)
parser.add_argument(
    "outfile", type=str, nargs="?", help="Output vcf file", default=outputfile
)

args = parser.parse_args()
outgroup = args.outgroup

if outgroup is None:
    logging.error("Need at least one outgroup")
    sys.exit(1)

if args.n is None or args.n < 1 or args.n > len(outgroup) * args.ploidy:
    n_vote = len(outgroup) * args.ploidy
else:
    n_vote = args.n

vcf = cyvcf2.VCF(args.vcf)
# Make sure outgroups in samples
if not set(outgroup) <= set(vcf.samples):
    logging.error("not all outgroup names found in vcf sample header")
    sys.exit(1)


n_pass = 0
n_skip = 0
n_change = 0
outgroup_indices = [i for i in range(len(vcf.samples)) if vcf.samples[i] in outgroup]
outext = re.compile(r"(.vcf|.vcf.gz|.bcf)$").search(args.outfile)
if outext is None:
    logging.error("please use file extension .vcf, .vcf.gz, or .bcf")
    sys.exit(1)
outext = outext.group(1)
mode = {".vcf": "w", ".vcf.gz": "wz", ".bcf": "wb"}[outext]


vcf.add_info_to_header(
    {"ID": "AA", "Number": 1, "Type": "Character", "Description": "Ancestral Allele"}
)
vcfwriter = cyvcf2.cyvcf2.Writer(args.outfile, vcf, mode)
for variant in tqdm(vcf):
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

if isinstance(logfile, str):
    fh = open(logfile, "w")
else:
    fh = logfile

fh.write("get_ancestral_allel.py summary\n")
fh.write("------------------------------\n")
fh.write(f"outgroups: {outgroup}\n")
fh.write(f"number of concordant alleles required in outgroups (n): {n_vote}\n")
fh.write(
    f"number of sites passing vote (either at least"
    f"{n_vote} 0 or {n_vote} 1 calls) (n): {n_pass}\n"
)
fh.write(f"number of sites changed ref (>={n_vote} 1 calls) (n): {n_change}\n")
fh.write(f"number of skipped sites (n): {n_skip}\n")

fh.close()
