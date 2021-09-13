rule bcftools_index:
    """Index a vcf file"""
    output:
        "{prefix}.vcf.gz.{ext}",
    input:
        "{prefix}.vcf.gz",
    wildcard_constraints:
        ext="(csi|tbi)",
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        *cfg.ruleconf("bcftools_index").params("envmodules"),
    log:
        "logs/{prefix}.vcf.gz.{ext}.log",
    shell:
        "bcftools index --{wildcards.ext} {input} -o {output} > {log} 2>&1"


rule tsinfer_infer:
    """Run tsinfer.infer on a group of samples

    Generate tskit tree file.
    """
    output:
        trees="{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees",
        samples="{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.samples",
    input:
        vcf=__RAW__ / "variants/{dataset}/{prefix}_{chrom}_{suffix}.vcf.gz",
        csi=__RAW__ / "variants/{dataset}/{prefix}_{chrom}_{suffix}.vcf.gz.csi",
    params:
        length=lambda wildcards: refdict[wildcards.chrom],
        samples=lambda wildcards: cfg.get_analysis(wildcards.analysis).samples.data,
    threads: 20
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_infer").params("envmodules"),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.log",
    script:
        "../scripts/tsinfer.py"


## FIXME: need to run as script so as to simplify tree sequence to get
## rid of unary and dangling nodes. Unary are removed with
## simplify(keep_unary=False) and dangling nodes are simply removed as
## a by-product(?)
rule tsinfer_tsdate:
    """Run tsdate on trees"""
    output:
        trees="{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.tsdate.{Ne}.trees",
    input:
        trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees",
    wildcard_constraints:
        Ne="\d+",
    params:
        options="-m 1e-8 --ignore-oldest ",
    threads: 20
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_tsdate").params("envmodules"),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.tsdate.{Ne}.log",
    shell:
        "tsdate date -v 1 -p -t {threads} {params.options} {input.trees} {output.trees} {wildcards.Ne} 2>&1 > {log}"


rule tsinfer_gnn:
    """Convert tree sequence to gnn"""
    output:
        gnn="{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.csv",
        mean="{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.mean.csv",
    input:
        trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees",
    threads: 1
    resources:
        mem_mb = xx_mem_mb
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_gnn").params("envmodules"),
    log:
        "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.log",
    script:
        "../scripts/tsinfer-gnn.py"

rule tsinfer_eda:
    """Make EDA document based on bokeh plots"""
    output:
        html="{results}/{analysis}/{dataset}/eda.html"
    input:
        csv=lambda wildcards: expand(expand("{{{{results}}}}/{{{{analysis}}}}/{{{{dataset}}}}/{fmt}.gnn.csv",
                                            fmt=cfg.get_analysis(wildcards.analysis).fmt),
                                     chrom=cfg.get_analysis(wildcards.analysis).chromosomes),
        trees=lambda wildcards: expand(expand(__INTERIM__ / "{{{{analysis}}}}/{{{{dataset}}}}/{fmt}.trees",
                                              fmt=cfg.get_analysis(wildcards.analysis).fmt),
                                       chrom=cfg.get_analysis(wildcards.analysis).chromosomes)
    params:
        fmt = lambda wildcards: cfg.get_analysis(wildcards.analysis).fmt
    conda:
        "../envs/plotting.yaml"
    log:
        "logs/{results}/{analysis}/{dataset}/eda.log",
    script:
        "../scripts/tsinfer-eda.py"


rule tsinfer_sample_gnn:
    """Calculate gnn data for a single sample"""
    output:
        csv="{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.{sample}.csv",
    input:
        trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees",
    threads: 1
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_sample_gnn").params("envmodules"),
    log:
        "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.{sample}.log",
    script:
        "../scripts/tsinfer-sample-gnn.py"


rule tsinfer_gnn_plot:
    """Make plots based on gnn output

    Generic plot function that can produce
    1. normalized gnn plot with dendrograms
    2. global gnn plot
    3. sample-based gnn plot

    The trees input file is needed to extract sample and population
    metadata and could possibly be precalculated in tsinfer gnn rules.
    """
    output:
        png=report("{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.png",
                   caption="../report/meangnn.rst", category="{dataset}",
                   subcategory="{analysis}: {datatype} plot")
    input:
        csv=__RESULTS__
        / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.csv",
        samples=lambda wildcards: "data/interim/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.samples".format(
            **dict(wildcards)
        ),
    resources:
        mem_mb = xx_mem_mb
    threads: 1
    wildcard_constraints:
        datatype="(mean|gnn)",
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_gnn_plot").params("envmodules"),
    log:
        "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.log",
    script:
        "../scripts/tsinfer-plot.py"


## FIXME: replace with python plotting function
rule tsinfer_gnn_plot_R:
    """Make gnn plots"""
    output:
        png=report("{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.R.gnn.png",
                   caption="../report/gnn.rst", category="{dataset}",
                   subcategory="{analysis}: Genealogical Nearest Neighbours plot")
    input:
        gnn=__RESULTS__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.csv",
        samples=config["samples"], #?!? should subset on analysis!
        populations=cfg.populations, # DITTO
    params:
        samples=lambda wildcards: cfg.get_analysis(
            wildcards.analysis
        ).samples.data.index,
    wildcard_constraints:
        label="(|.mean)",
    resources:
        mem_mb = x_mem_mb
    threads: 1
    conda:
        "../envs/tsinfer.yaml"
    envmodules:
        *cfg.ruleconf("tsinfer_gnn_plot_R").params("envmodules"),
    log:
        "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.R.gnn.png",
    script:
        "../scripts/plot-gnn.R"
