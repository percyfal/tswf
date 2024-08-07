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
        *get_envmodules(config, "bcftools_index"),
    log:
        "logs/{prefix}{bcf}.{ext}.log",
    shell:
        "bcftools index --{wildcards.ext} {input} -o {output} > {log} 2>&1"


# FIXME: input vcf should depend on whether derive_aa has been set or
# not
rule tsinfer_make_samples:
    """Make tsinfer sample data"""
    output:
        samples="{interim}/tsinfer/{analysis}/{dataset}/samples/{chrom}/{prefix}.samples",
    input:
        vcf="{interim}/variants/{dataset}/ancestral/{chrom}/{prefix}.vcf.gz",
        csi="{interim}/variants/{dataset}/ancestral/{chrom}/{prefix}.vcf.gz.csi",
        samples=str(config["samples"]),
        populations=str(config["populations"]),
    log:
        "logs/{interim}/{analysis}/{dataset}/samples/{chrom}/{prefix}.samples.log",
    params:
        length=get_chrom_length,
    threads: 1
    shell:
        "tswf-tsinfer-make-samples {input.vcf} {input.samples} {input.populations} {output.samples} --length {params.length} > {log} 2>&1"


rule tsinfer_infer:
    """Run tsinfer.infer on a group of samples

    Generate tskit tree file.
    """
    output:
        trees="{interim}/tsinfer/{analysis}/{dataset}/infer/{chrom}/{prefix}.trees",
    input:
        samples="{interim}/tsinfer/{analysis}/{dataset}/samples/{chrom}/{prefix}.samples",
    params:
        length=get_chrom_length,
    threads: 20
    log:
        "logs/{interim}/tsinfer/{analysis}/{dataset}/infer/{chrom}/{prefix}.log",
    shell:
        "tsinfer infer {input.samples} -O {output.trees} --progress --num-threads {threads} > {log} 2>&1"


rule tsinfer_gnn:
    """Convert tree sequence to gnn"""
    output:
        gnn="{results}/tsinfer/{analysis}/{dataset}/{chrom}/{mode}/{prefix}.gnn.csv",
    input:
        trees=__INTERIM__ / "tsinfer/{analysis}/{dataset}/infer/{chrom}/{prefix}.trees",
    wildcard_constraints:
        mode="(individual|population)",
    threads: 1
    log:
        "logs/{results}/tsinfer/{analysis}/{dataset}/{chrom}/{mode}/{prefix}.gnn.log",
    shell:
        "tswf-tsinfer-gnn {input.trees} --output-file {output.gnn} --mode {wildcards.mode} > {log} 2>&1"


rule tsinfer_eda:
    """Make EDA document based on bokeh plots"""
    output:
        html="{results}/tsinfer/{analysis}/{dataset}/{mode}.eda.html",
    input:
        unpack(tsinfer_eda_input),
    params:
        gnn=lambda wildcards, input: " ".join([f"--gnn {x}" for x in input.gnn]),
        trees=lambda wildcards, input: " ".join([f"--ts {x}" for x in input.trees]),
    log:
        "logs/{results}/tsinfer/{analysis}/{dataset}/{mode}.eda.log",
    threads: 1
    shell:
        "tswf-tsinfer-eda {params.gnn} {params.trees} --output-file {output.html} > {log} 2>&1"


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
