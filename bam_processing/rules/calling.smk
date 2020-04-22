###############################################################################
rule call_variants:
    input:
        bam= join(BAM_DIR, "{sample}_recal.bam")
        ref=config["ref"]["genome"],
        # known=config["ref"]["known-variants"],
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    shell: """
        gatk HaplotypeCaller \
            --reference {input.ref} \
            --input {input.bam} \
            --output {output} \
            -ERC GVCF
    """

###############################################################################
rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("--variant called/{sample}.g.vcf.gz", sample=BAM_PREFIX)
    output:
        gvcf="called/all.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.log"
     shell: """
        gatk CombineGVCFs \
            --reference {REF_FILE} \
            {params.gvcf_string} \
            --output {output} \
    """

###############################################################################
rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.vcf.gz")
    log:
        "logs/gatk/genotypegvcfs.log"
    shell: """
		gatk GenotypeGVCFs \
            --reference {input.ref} \
            --variant {input.gvcf} \
            --output {output} \
    """

###############################################################################
rule variants2zarr:
    input:
        vcf=rules.merge_variants.output.vcf
    output:
        vcf="genotyped/all.zarr"
    log:
        "logs/zarr/merged-zarr.log"
    script:
        "rules/zarr_format.py {input} {output}"
