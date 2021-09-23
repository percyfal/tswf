#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import os
import sys
import math
try:
    import tsinfer
except ModuleNotFoundError as e:
    raise(ModuleNotFoundError("No module named 'tsinfer'"))
import cyvcf2
import json as json
import pandas as pd
from tqdm import tqdm
import logging
FORMAT = '%(levelname)s:tsinfer:%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('tsinfer')


inputvcf = snakemake.input.vcf
inputcsi = snakemake.input.csi
df_meta = snakemake.params.samples
df_pop = snakemake.params.populations
chrom = snakemake.wildcards.chrom
chromlength = snakemake.params.length
tsout = snakemake.output.trees

print(f"...read {df_meta.shape[0]} entries\n")


def add_diploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    for variant in tqdm(vcf):  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        ancestral = variant.INFO.get("AA", variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        # Both alleles may be derived and currently seems to cause the
        # add_site function to fail; skip for now
        if len(ordered_alleles) > 2:
            continue
        allele_index = {
            old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)
        }
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:2]
        ]
        # Subsetting the vcf may mean we are looking at a monomorphic site; skip
        if variant.nucl_diversity == 0.0:
            continue
        try:
            samples.add_site(pos, genotypes=genotypes, alleles=alleles)
        except ValueError as e:
            print(pos, genotypes, alleles)
            print(e)
            raise


def add_populations(vcf, samples):
    """Add numeric population id for each sample to samples"""
    vcf_samples = pd.DataFrame({"SM": vcf.samples})
    if not all(vcf_samples.SM.isin(df_meta.index)):
        logger.warning("Some samples in vcf file not present in sample metadata")
    try:
        df_filt = (
            vcf_samples.merge(df_meta, left_on="SM", right_on="SM")
            .set_index(["SM"], drop=True)
            .loc[:, :]
        )
    except KeyError as e:
        print(e)
        raise
    except Exception as e:
        raise
    sample_pops = list(df_filt.population)
    pop_lookup = {}
    for pop in sorted(set(sample_pops)):
        md = df_pop.loc[pop]
        md = {k: v for k, v in md.iteritems()}
        md['population'] = pop
        pop_lookup[pop] = samples.add_population(metadata=md)
    return [pop_lookup[pop] for pop in sample_pops]


def add_diploid_individuals(vcf, samples, populations):
    """Add sample name and population id to samples"""
    for name, population in zip(vcf.samples, populations):
        md = {"name": name, "SM": name}
        md.update(dict(df_meta.loc[name]))
        samples.add_individual(ploidy=2, metadata=md, population=population)

def init_vcf(fn, samplenames):
    vcf = cyvcf2.VCF(fn)
    try:
        vcf.set_samples(samplenames)
    except Exception as e:
        print(e)
        raise
    return vcf


# Init vcf
vcf = init_vcf(inputvcf, df_meta.index.to_list())

samples_path = re.sub(".trees$", ".samples", tsout)

with tsinfer.SampleData(path=samples_path, sequence_length=chromlength) as samples:
    populations = add_populations(vcf, samples)
    add_diploid_individuals(vcf, samples, populations)
    add_diploid_sites(vcf, samples)

print(
    "Sample file created for {} samples ".format(samples.num_samples)
    + "({} individuals) ".format(samples.num_individuals)
    + "with {} variable sites.".format(samples.num_sites),
    flush=True,
)


ts = tsinfer.infer(samples, progress_monitor=True, num_threads=snakemake.threads)

print(
    "Inferred tree sequence `{}`: {} trees over {} Mb".format(
        "ts", ts.num_trees, ts.sequence_length / 1e6
    )
)

ts.dump(tsout)
