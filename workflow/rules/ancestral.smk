rule derive_ancestral_by_vote:
    output: vcf = __RAW__ / "variants/{dataset}/{prefix}_AA_vote.vcf.gz"
    input: vcf = __RAW__ / "variants/{dataset}/{prefix}.vcf.gz"
    params:
        outgroup = cfg.derive_aa.outgroups,
        options = dict(cfg.derive_aa.options)
    conda:
        "../envs/tsinfer.yaml"
    log:
        "logs/variants/{dataset}/{prefix}_AA_vote.log"
    threads: 1
    script:
        "../scripts/get_ancestral_allele.py"
