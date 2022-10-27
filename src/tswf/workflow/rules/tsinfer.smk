rule bcftools_index:
    """Index a vcf file"""
    output:
        "{prefix}{bcf}.{ext}",
    input:
        "{prefix}{bcf}",
    wildcard_constraints:
        ext="(csi|tbi)",
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        *envmodules.get("bcftools_index", []),
    log:
        "logs/{prefix}{bcf}.{ext}.log",
    shell:
        "bcftools index --{wildcards.ext} {input} -o {output} > {log} 2>&1"


rule tsinfer_infer:
    """Run tsinfer.infer on a group of samples

    Generate tskit tree file.
    """
    output:
        trees="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.trees",
        samples="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.samples",
    input:
        vcf=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}.vcf.gz",
        csi=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}.vcf.gz.csi",
        samples=str(config["samples"]),
        populations=str(config["populations"]),
    params:
        length=lambda wildcards: refdict[wildcards.chrom],
    threads: 20
    envmodules:
        *envmodules.get("tsinfer_infer", []),
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.log",
    shell:
        "tswf-tsinfer {input.vcf} {input.samples} {input.populations} {output.trees} --length {params.length}"


## FIXME: need to run as script so as to simplify tree sequence to get
## rid of unary and dangling nodes. Unary are removed with
## simplify(keep_unary=False) and dangling nodes are simply removed as
## a by-product(?)
# rule tsinfer_tsdate:
#     """Run tsdate on trees"""
#     output:
#         trees="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.tsdate.{Ne}.trees",
#     input:
#         trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.trees",
#     wildcard_constraints:
#         Ne="[0-9]+",
#     params:
#         options="-m 1e-8 --ignore-oldest ",
#     threads: 20
#     conda:
#         "../envs/tsinfer.yaml"
#     envmodules:
#         *cfg.ruleconf("tsinfer_tsdate").params("envmodules"),
#     log:
#         "logs/{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.tsdate.{Ne}.log",
#     shell:
#         "tsdate date -v 1 -p -t {threads} {params.options} {input.trees} {output.trees} {wildcards.Ne} 2>&1 > {log}"
# rule tsinfer_gnn:
#     """Convert tree sequence to gnn"""
#     output:
#         gnn="{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.gnn.{mode}.csv",
#     input:
#         trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.trees",
#     wildcard_constraints:
#         mode="(individual|population)",
#     threads: 1
#     resources:
#         mem_mb=xx_mem_mb,
#     conda:
#         "../envs/tsinfer.yaml"
#     envmodules:
#         *cfg.ruleconf("tsinfer_gnn").params("envmodules"),
#     log:
#         "logs/{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.gnn.{mode}.log",
#     script:
#         "../scripts/tsinfer-gnn.py"
# rule tsinfer_gnn_archive:
#     """Archive csv output files"""
#     output:
#         targz="{results}/{analysis}/{dataset}/gnn.{mode}.tar.gz",
#     input:
#         unpack(tsinfer_eda_input),
#     log:
#         "logs/{results}/{analysis}/{dataset}/gnn.{mode}.tar.gz.log",
#     threads: 1
#     shell:
#         "tar -zcvf {output.targz} {input.csv} > {log}"
# rule tsinfer_eda:
#     """Make EDA document based on bokeh plots"""
#     output:
#         html="{results}/{analysis}/{dataset}/{mode}.eda.html",
#     input:
#         unpack(tsinfer_eda_input),
#     params:
#         fmt=fmt,
#     conda:
#         "../envs/plotting.yaml"
#     log:
#         "logs/{results}/{analysis}/{dataset}/{mode}.eda.log",
#     threads: 1
#     script:
#         "../scripts/tsinfer-eda.py"
# rule tsinfer_sample_gnn:
#     """Calculate gnn data for a single sample"""
#     output:
#         csv="{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.gnn.{sample}.csv",
#     input:
#         trees=__INTERIM__ / "{analysis}/{dataset}/{prefix}{chrom}{suffix}.trees",
#     threads: 1
#     conda:
#         "../envs/tsinfer.yaml"
#     envmodules:
#         *cfg.ruleconf("tsinfer_sample_gnn").params("envmodules"),
#     log:
#         "logs/{results}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.gnn.{sample}.log",
#     script:
#         "../scripts/tsinfer-sample-gnn.py"
