rule bcftools_index:
    """Index a vcf file"""
    output: "{prefix}.vcf.gz.{ext}"
    input: "{prefix}.vcf.gz"
    wildcard_constraints:
        ext = "(csi|tbi)"
    conda: "../envs/bcftools.yaml"
    envmodules: *cfg.ruleconf("bcftools_index").params("envmodules")
    log: "logs/{prefix}.vcf.gz.{ext}.log"
    shell:
        "bcftools index --{wildcards.ext} {input} -o {output} > {log} 2>&1"



rule tsinfer_infer:
    """Run tsinfer.infer on a group of samples

    Generate tskit tree file.
    """
    output: trees = "{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees",
            samples = "{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.samples"
    input: vcf = __RAW__ / "variants/{dataset}/{prefix}_{chrom}_{suffix}.vcf.gz",
           csi = __RAW__ / "variants/{dataset}/{prefix}_{chrom}_{suffix}.vcf.gz.csi",
    params:
        length = lambda wildcards: refdict[wildcards.chrom],
        samples = lambda wildcards: cfg.get_analysis(wildcards.analysis).samples.data
    threads: lambda wildcards: cfg.ruleconf("map_bwa_index", analysis=wildcards.analysis).threads
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_infer").params("envmodules")
    log: "logs/{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.log"
    script:
        "../scripts/tsinfer.py"

## FIXME: need to run as script so as to simplify tree sequence to get
## rid of unary and dangling nodes. Unary are removed with
## simplify(keep_unary=False) and dangling nodes are simply removed as
## a by-product(?)
rule tsinfer_tsdate:
    """Run tsdate on trees"""
    output: trees = "{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.tsdate.{Ne}.trees"
    input: trees = __INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees"
    wildcard_constraints: Ne = "\d+"
    params:
        options = "-m 1e-8 --ignore-oldest "
    threads: lambda wildcards: cfg.ruleconf("tsinfer_tsdate", analysis=wildcards.analysis).threads
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_tsdate").params("envmodules")
    log: "logs/{interim}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.tsdate.{Ne}.log"
    shell:
        "tsdate date -v 1 -p -t {threads} {params.options} {input.trees} {output.trees} {wildcards.Ne} 2>&1 > {log}"


rule tsinfer_gnn:
    """Convert tree sequence to gnn"""
    output:
        gnn = "{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.csv",
        mean = "{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.mean.csv"
    input:
        trees = __INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees"
    threads: 1
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_gnn").params("envmodules")
    log: "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.log"
    script:
        "../scripts/tsinfer-gnn.py"


rule tsinfer_sample_gnn:
    """Calculate gnn data for a single sample"""
    output:
        csv = "{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.{sample}.csv"
    input:
        trees = __INTERIM__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees"
    threads: 1
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_sample_gnn").params("envmodules")
    log: "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.{sample}.log"
    script:
        "../scripts/tsinfer-sample-gnn.py"


rule tsinfer_gnn_plot:
    """Make plots based on gnn output

    Generic plot function that can produce
    1. normalized gnn plot with dendrograms
    2. global gnn plot
    3. sample-based gnn plot
    """
    output: png = "{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.png"
    input: csv = __RESULTS__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.csv",
           trees = lambda wildcards: "data/interim/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.trees".format(**dict(wildcards))
    threads: 1
    wildcard_constraints:
        datatype = "(mean|gnn)"
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_gnn_plot").params("envmodules")
    log: "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.{datatype}{dot}{sample}.log"
    script:
        "../scripts/tsinfer-plot.py"


## FIXME: replace with python plotting function
rule tsinfer_gnn_plot_R:
    """Make gnn plots"""
    output: png = "{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.R.gnn.png",
    input: gnn = __RESULTS__ / "{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.gnn.csv",
           samples = config["samples"],
           populations = cfg.populations
    params:
        samples = lambda wildcards: cfg.get_analysis(wildcards.analysis).samples.data.index
    wildcard_constraints:
        label = "(|.mean)"
    threads: 1
    conda: "../envs/tsinfer.yaml"
    envmodules: *cfg.ruleconf("tsinfer_gnn_plot_R").params("envmodules")
    log: "logs/{results}/{analysis}/{dataset}/{prefix}_{chrom}_{suffix}.R.gnn.png"
    script:
        "../scripts/plot-gnn.R"
