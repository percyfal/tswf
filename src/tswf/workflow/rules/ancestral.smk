rule derive_ancestral_by_vote:
    output:
        vcf=__INTERIM__ / "variants/{dataset}/ancestral/{chrom}/{prefix}{bcf}",
    input:
        vcf=__RAW__ / "variants/{dataset}/{chrom}/{prefix}{bcf}",
        tbi=__RAW__ / "variants/{dataset}/{chrom}/{prefix}{bcf}.tbi",
    params:
        outgroup=cfg.derive_aa.outgroups,
        options=lambda wildcards: " ".join(
            [f"--{k} {v}" for k, v in dict(cfg.derive_aa.options).items()]
        ),
    log:
        "logs/variants/{dataset}/ancestral/{chrom}/{prefix}{bcf}.log",
    threads: 1
    shell:
        "tswf-get-ancestral-allele {input.vcf} {params.outgroup} --outfile {output.vcf} {params.options}"


rule ancestralize_reference_sequence:
    """NOT_IMPLEMENTED: Given a reference sequence and vcf file with AA
    tag, convert reference to ancestral state"""
    output:
        fa="{interim}/{analysis}/{dataset}/{chrom}/{prefix}.{mode}{bcf}.fasta",
    input:
        fa=cfg.ref,
        vcf=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}_AA_{mode}{bcf}",
    log:
        "logs/{interim}/{analysis}/{dataset}/{chrom}/{prefix}.{mode}{bcf}.fasta.log",
    threads: 1
    shell:
        "tswf-make-ancestral-sequence {input.fa} {input.vcf} {output.fa} --chromosome {wildcards.chrom}"
