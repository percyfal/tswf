rule vcf2zarr_explode:
    """Convert vcf file to zarr format."""
    output:
        icf="{prefix}.icf",
    input:
        vcf="{prefix}.vcf.gz",
        tbi="{prefix}.vcf.gz.tbi",
    envmodules:
        *get_envmodules(config, 'bio2zarr'),
    conda:
        "../envs/bio2zarr.yaml"
    benchmark:
        "benchmarks/vcf2zarr_explode/{prefix}.zarr.benchmark.txt"
    log:
        "logs/vcf2zarr_explode/{prefix}.zarr.log",
    threads: 1
    shell:
        """vcf2zarr explode {input.vcf} {output.icf}"""


rule vcf2zarr_encode:
    """Encode icf file to zarr format."""
    output:
        "{prefix}.zarr",
    input:
        "{prefix}.icf",
    envmodules:
        *get_envmodules(config, 'bio2zarr'),
    conda:
        "../envs/bio2zarr.yaml"
    benchmark:
        "benchmarks/vcf2zarr_encode/{prefix}.zarr.benchmark.txt"
    log:
        "logs/vcf2zarr_encode/{prefix}.zarr.log",
    threads: 1
    shell:
        """vcf2zarr encode {input.icf} {output.zarr}"""
