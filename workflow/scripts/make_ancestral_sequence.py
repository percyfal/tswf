#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Recreate ancestral sequence from a vcf with AA INFO tag
#
import os
import sys
import re
import argparse
import logging
from pyfaidx import Fasta
import cyvcf2
from tqdm import tqdm


try:
    inref = snakemake.input.fa
    outref = snakemake.output.fa
    vcf = snakemake.input.vcf
    logfile = str(snakemake.log)
    chromosome = snakemake.params.chromosome
except NameError as e:
    inref = None
    outref = None
    vcf = None
    logfile = sys.stdout
    chromosome = None

parser = argparse.ArgumentParser(
    description=""" Make ancestral reference sequence given reference sequence and a
vcf file with the AA INFO tag denoting the ancestral allele. """
)
parser.add_argument("reference", type=str, help="Input reference file", default=inref)
parser.add_argument("vcf", type=str, help="Input vcf file", default=vcf)
parser.add_argument("outfile", type=str, help="Output vcf file", default=outref)

args = parser.parse_args()


vcf = cyvcf2.VCF(args.vcf)

if chromosome is None:
    chromosome = vcf.seqnames[0]


ref = Fasta(args.reference)
try:
    seq = ref[chromosome]
except KeyError as e:
    print(e)
    raise

n_sites = 0
n_change = 0
# Loop vcf and update reference at positions
for variant in tqdm(vcf):
    if not "AA" in dict(variant.INFO).keys():
        logging.error(
            "No AA tag in INFO field; please add AA tag that indicates ancestral state"
        )
        sys.exit(1)
    n_sites = n_sites + 1
    print(variant.INFO["AA"], variant.REF)
    if variant.INFO["AA"] != variant.REF:
        n_change = n_change + 1
        print(
            f"Changing REF {variant.REF} ({reference[variant.POS-1]}) to {variant.INFO['AA']}"
        )
        seq[variant.POS - 1] = variant.INFO["AA"]
