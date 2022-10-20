# FIXME: add relate config
# RELATE_PATH = Path(cfg.relate.path)absolute()
RELATE_PATH = Path("opt/relate/relate").absolute()
RELATE_BIN = RELATE_PATH / "bin"
RELATE_SCRIPTS = RELATE_PATH / "scripts"


# Must first reorder alleles such that ancestral is first
rule relate_convert_from_vcf:
    """Convert vcf file to shapeit format"""
    output:
        haps="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}{bcf}.haps",
        sample="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}{bcf}.sample",
    input:
        bcf=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}{bcf}",
        csi=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}{bcf}.csi",
    params:
        cmd=RELATE_BIN / "RelateFileFormats",
        vcf=lambda wildcards: __RAW__
        / "variants/{dataset}/{prefix}{chrom}{suffix}".format(**wildcards),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}{bcf}.convertfromvcf.log",
    shell:
        "{params.cmd} --mode ConvertFromVcf --haps {output.haps} --sample {output.sample} -i {params.vcf} 2&>1 > {log}"


rule relate_prepare_samples:
    """Prepare samples using a repolarized version of the reference"""
    output:
        haps="foo",


# FIXME: somehow generate map file *or* simple require as input
rule relate_run:
    output:
        anc="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.anc",
        mut="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.mut",
    input:
        haps="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.haps",
        sample="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.sample",
        map=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}.map",
    params:
        #options=lambda wildcards: cfg.get_analysis(wildcards.analysis).relate.options,  #"-m 0.1e-8 -N 180000 --seed 42",
        options=cfg.ruleconf("relate_run").params("options"),
        cmd=RELATE_BIN / "Relate",
        outdir=lambda wildcards: "{interim}/{analysis}/{dataset}".format(**wildcards),
        tmpdir=lambda wildcards: "{prefix}{chrom}{suffix}".format(**wildcards),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.log",
    shell:
        "{params.cmd} --mode All {params.options} --haps {input.haps} --sample {input.sample} --map {input.map} -o {params.tmpdir}  2&>1 > {log}; "
        "mv {params.tmpdir}/*.anc {params.outdir}/;"
        "mv {params.tmpdir}/*.mut {params.outdir}/;"
        "rmdir {params.tmpdir}"


rule relate_estimate_population_size:
    output:
        pdf="{prefix}_popsize.pdf",
        coal="{prefix}_popsize.coal",
        pairwise="{prefix}_popsize.pairwise.coal",
        rate="{prefix}_popsize_avg.rate",
        anc="{prefix}_popsize.anc.gz",
        mut="{prefix}_popsize.mut.gz",
        dist="{prefix}_popsize.dist",
    input:
        anc="{prefix}.anc",
        mut="{prefix}.mut",
        poplabels="{prefix}.poplabels",
    params:
        options="-m 0.1e-8 --seed 42 --threshold 0 --years_per_gen 2",
        cmd=RELATE_SCRIPTS / "EstimatePopulationSize/EstimatePopulationSize.sh",
    log:
        "logs/{prefix}.relate.populationsize.log",
    threads: 20
    shell:
        "{params.cmd} --threads {threads} {params.options} -i {wildcards.prefix} --poplabels {input.poplabels} -o {wildcards.prefix}_popsize 2&>1 > {log}"


rule relate_detect_selection:
    output:
        lin="{prefix}_popsize_selection.lin",
        freq="{prefix}_popsize_selection.freq",
        sele="{prefix}_popsize_selection.sele",
    input:
        anc="{prefix}_popsize.anc.gz",
        mut="{prefix}_popsize.mut.gz",
    params:
        options="-m 0.1e-8 --years_per_gen 2",
        cmd=RELATE_SCRIPTS / "DetectSelection/DetectSelection.sh",
    log:
        "logs/{prefix}.relate.populationsize.selection.log",
    shell:
        "{params.cmd} {params.options} -i {wildcards.prefix}_popsize -o {wildcards.prefix}_popsize_selection"


rule relate_treeview:
    output:
        pdf="{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.{bp}.pdf",
    input:
        anc=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.anc",
        mut=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.mut",
        sample=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.sample",
        haps=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.haps",
        poplabels=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.poplabels",
    params:
        cmd=RELATE_SCRIPTS / "TreeView/TreeView.sh",
        options="--years_per_gen 2",
        outpfx=lambda wildcards: __INTERIM__
        / "{analysis}/{dataset}/{prefix}{chrom}{suffix}".format(**wildcards),
    wildcard_constraints:
        bp="[0-9]+",
    log:
        "logs/{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}_popsize.{bp}.pdf.log",
    shell:
        "{params.cmd} {params.options} --anc {input.anc} --mut {input.mut} --sample {input.sample} --haps {input.haps} --poplabels {input.poplabels} --bp_of_interest {wildcards.bp} -o {params.outpfx} 2>&1 > {log};"
        "mv {params.outpfx}.pdf {output.pdf}"


rule relate_make_poplabels:
    """Make population labels"""
    output:
        poplabels="{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.poplabels",
    input:
        sample="{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.sample",
        metadata=lambda wildcards: config[wildcards.analysis].get(
            "metadata", config["metadata"]
        ),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.poplabels.log",
    script:
        "../scripts/relate_make_poplabels.py"
