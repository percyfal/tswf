rule derive_ancestral_by_vote:
    output:
        vcf=__RAW__ / "variants/{dataset}/{prefix}_AA_vote{bcf}",
    input:
        vcf=__RAW__ / "variants/{dataset}/{prefix}{bcf}",
        tbi=__RAW__ / "variants/{dataset}/{prefix}{bcf}.tbi",
    params:
        outgroup=cfg.get("derive_aa", {}).get("outgroups", []),
        options=dict(cfg.get("derive_aa", {}).get("options", {})),
    conda:
        "../envs/tsinfer.yaml"
    log:
        "logs/variants/{dataset}/{prefix}_AA_vote{bcf}.log",
    threads: 1
    script:
        "../scripts/get_ancestral_allele.py"
