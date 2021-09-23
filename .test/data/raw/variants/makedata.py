#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Simulate test data
#
import os
import math
import msprime
import cyvcf2

# Times are provided in years, so we convert into generations.
generation_time = 25
T_OOA = 21.2e3 / generation_time
T_AMH = 140e3 / generation_time
T_ANC = 220e3 / generation_time
T_PAN = 5.6e6 / generation_time
T_GOR = 8.6e6 / generation_time
T_APE = 18.0e6 / generation_time
# Provide some time for ancestral state to differ from current
T_PREAPE = 25e6 / generation_time

# We need to work out the starting population sizes based on
# the growth rates provided for these two populations
r_CEU = 0.004
r_CHB = 0.0055
N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
N_CHB = 510 / math.exp(-r_CHB * T_OOA)
N_AFR = 12300
N_ANC = 7300


# Add human populations
demography = msprime.Demography()
demography.add_population(
    name="YRI",
    description="Yoruba in Ibadan, Nigeria",
    initial_size=N_AFR,
)
demography.add_population(
    name="CEU",
    description=("Utah Residents (CEPH) with Northern and Western European Ancestry"),
    initial_size=N_CEU,
    growth_rate=r_CEU,
)
demography.add_population(
    name="CHB",
    description="Han Chinese in Beijing, China",
    initial_size=N_CHB,
    growth_rate=r_CHB,
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa population",
    initial_size=2100,
)
demography.add_population(
    name="AMH", description="Anatomically modern humans", initial_size=N_AFR
)
demography.add_population(
    name="ANC",
    description="Ancestral equilibrium population",
    initial_size=N_ANC,
)

# Add ape populations
demography.add_population(
    name="chimpanzee",
    description="chimpanzee",
    initial_size=N_ANC,
)
demography.add_population(
    name="PAN",
    description="ancestral population to chimps and humans",
    initial_size=N_ANC,
)
demography.add_population(
    name="gorilla",
    description="gorilla",
    initial_size=N_ANC,
)
demography.add_population(
    name="GOR",
    description="ancestral population to PAN and gorillas",
    initial_size=N_ANC,
)
demography.add_population(
    name="orangutan",
    description="orangutan",
    initial_size=N_ANC,
)
demography.add_population(
    name="APE",
    description="ancestral population to GOR and orangutan",
    initial_size=N_ANC,
)
demography.add_population(
    name="PREAPE",
    description="ancestral population to ANC",
    initial_size=N_ANC,
)

# Add OOA migration rates
demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

# Add events
demography.add_population_split(time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA")
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["YRI", "OOA"], rate=25e-5
)
demography.add_population_split(time=T_AMH, derived=["YRI", "OOA"], ancestral="AMH")
demography.add_population_split(time=T_ANC, derived=["AMH"], ancestral="ANC")

# Now add the events for the ape populations:
demography.add_population_split(
    time=T_PAN, derived=["ANC", "chimpanzee"], ancestral="PAN"
)
demography.add_population_split(time=T_GOR, derived=["PAN", "gorilla"], ancestral="GOR")
demography.add_population_split(
    time=T_APE, derived=["GOR", "orangutan"], ancestral="APE"
)
demography.add_population_split(time=T_PREAPE, derived=["APE"], ancestral="PREAPE")

# Simulate a chromosome
ts = msprime.sim_ancestry(
    {"CHB": 7, "CEU": 6, "YRI": 6, "chimpanzee": 1, "gorilla": 1, "orangutan": 1},
    demography=demography,
    random_seed=10,
    recombination_rate=1e-8,
    sequence_length=1e6,
)

# Add mutations
mts = msprime.sim_mutations(ts, rate=1e-9, random_seed=42)

print("Added ", mts.num_sites, " sites")


# Set individual names
indv_names = (
    ["tsk_REF"]
    + [f"tsk_CHB_{str(i)}" for i in range(6)]
    + [f"tsk_CEU_{str(i)}" for i in range(6)]
    + [f"tsk_YRI_{str(i)}" for i in range(6)]
    + ["tsk_chimp", "tsk_gorilla", "tsk_orangutan"]
)

curdir = os.path.abspath(os.path.dirname(__file__))
fn = os.path.join(curdir, "ooa/ooa_1_PASS.vcf")
with open(fn, "w") as fh:
    mts.write_vcf(fh, individual_names=indv_names)

# in cyvcf2: Pick one individual *haplotype* as reference: this
# individual should have only 0's, so all calls at a site with a
# derived allele must be flipped for all individuals.
refseq = cyvcf2.VCF(fn, samples=indv_names[0])
vcf = cyvcf2.VCF(fn, samples=indv_names[1:])
vcfout = f"{fn}.gz"
vcfwriter = cyvcf2.cyvcf2.Writer(vcfout, vcf, "wz")
vcfwriter.set_samples(vcf.samples)

for gt, variant in zip(refseq, vcf):
    ref = gt.genotypes[0][0]
    if ref == 1:
        # Flip states for all other genotypes - or randomly?
        for i in range(len(variant.genotypes)):
            alleles = variant.genotypes[i]
            variant.genotypes[i] = [1 - alleles[0], 1 - alleles[1], alleles[2]]
        variant.genotypes = variant.genotypes
        variant.REF = gt.ALT[0]
        variant.ALT = [gt.REF]
    vcfwriter.write_record(variant)
vcfwriter.close()
