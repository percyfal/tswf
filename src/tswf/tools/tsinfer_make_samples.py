"""Make tsinfer sample data.

Make tsinfer sample data output file SAMPLEDATA from VCF file using
metadata in SAMPLESFILE and POPULATIONFILE files.

"""

import logging

import click
import cyvcf2
import pandas as pd
import tsinfer
from tqdm import tqdm


FORMAT = "%(levelname)s:tsinfer:%(asctime)-15s %(message)s"
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("tsinfer")


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
        # Ancestral state must be first in the allele list. NB: can be missing.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        # Both alleles may be derived and currently seems to cause the
        # add_site function to fail; skip for now. This should work.
        if len(ordered_alleles) > 2:
            continue
        # Skip missing data for now
        if any(index == -1 for row in variant.genotypes for index in row[0:2]):
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


def add_populations(vcf, samples, df_meta, df_pop):
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
    except Exception:
        raise
    sample_pops = list(df_filt.population)
    pop_lookup = {}
    for pop in sorted(set(sample_pops)):
        md = df_pop.loc[pop].to_dict()
        md = {k: v for k, v in md.items()}
        md["population"] = pop
        pop_lookup[pop] = samples.add_population(metadata=md)
    return [pop_lookup[pop] for pop in sample_pops]


def add_diploid_individuals(vcf, samples, populations, df_meta):
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


@click.command(help=__doc__)
@click.argument("vcf", type=click.Path(exists=True))
@click.argument("samplesfile", type=click.Path(exists=True))
@click.argument("populationfile", type=click.Path(exists=True))
@click.argument("sampledata", type=click.Path(exists=False))
@click.option("--threads", type=int, default=1)
@click.option("chromlength", "--length", type=int)
def main(vcf, samplesfile, populationfile, sampledata, chromlength, threads):
    df_meta = pd.read_table(samplesfile, comment="#").set_index("SM")
    if df_meta.index.duplicated().any():
        logger.error("duplicate sample ids in samplesheet; sample names must be unique")
        exit(1)
    df_pop = pd.read_table(populationfile, comment="#").set_index("population")
    if df_pop.index.duplicated().any():
        logger.error(
            "duplicate population ids in population datasheet; "
            "population identifiers must be unique"
        )
        exit(1)
    print(f"...read {df_meta.shape[0]} entries\n")
    vcf = init_vcf(vcf, df_meta.index.to_list())

    with tsinfer.SampleData(path=sampledata, sequence_length=chromlength) as samples:
        populations = add_populations(vcf, samples, df_meta, df_pop)
        add_diploid_individuals(vcf, samples, populations, df_meta)
        add_diploid_sites(vcf, samples)

    print(
        f"Sample file created for {samples.num_samples} samples "
        + f"({samples.num_individuals} individuals) "
        + f"with {samples.num_sites} variable sites.",
        flush=True,
    )
