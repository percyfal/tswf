rule derive_ancestral_by_vote:
    output:
        vcf=__RAW__ / "variants/{dataset}/{prefix}_AA_vote{bcf}",
    input:
        vcf=__RAW__ / "variants/{dataset}/{prefix}{bcf}",
        tbi=__RAW__ / "variants/{dataset}/{prefix}{bcf}.tbi",
    params:
        outgroup=cfg.derive_aa.outgroups,
        options=dict(cfg.derive_aa.options),
    conda:
        "../envs/tsinfer.yaml"
    log:
        "logs/variants/{dataset}/{prefix}_AA_vote{bcf}.log",
    threads: 1
    script:
        "../scripts/get_ancestral_allele.py"



rule ancestralize_reference_sequence:
    """Given a reference sequence and vcf file with AA tag, convert reference to ancestral state"""
    output:
        fa="{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.{mode}{bcf}.fasta"
    input:
        fa=cfg.ref,
        vcf=__RAW__ / "variants/{dataset}/{prefix}{chrom}{suffix}_AA_{mode}{bcf}"
    conda:
        "../envs/seqio.yaml"
    log:
        "logs/{interim}/{analysis}/{dataset}/{prefix}{chrom}{suffix}.{mode}{bcf}.fasta.log"
    threads: 1
    script:
        "../scripts/make_ancestral_sequence.py"
