
rule call_variants:
    input:
        bam= join(BAM_DIR, "{sample}_recal.bam")
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra=get_call_variants_params
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/{sample}.g.vcf.gz", sample=BAM_PREFIX)
    output:
        gvcf="called/all.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        ref=config["ref"]["genome"] + ".fai"
        vcfs=rules.genotype_variants.output.vcf
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"


rule variants2zarr:
    input:
        vcf=rules.merge_variants.output.vcf
    output:
        vcf="genotyped/all.zarr"
    log:
        "logs/zarr/merged-zarr.log"
    script:
        "rules/zarr_format.py {input} {output}"
