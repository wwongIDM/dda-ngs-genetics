###############################################################################
rule call_variants:
    input:
        bam= join(config['bam_dir'], "{sample}_recal.bam"),
        ref=config["ref"]["genome"],
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
RUN_DIR = os.getcwd()
rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/{sample}.g.vcf.gz", sample=BAM_PREFIX)
    output:
        gvcf="called/all.g.vcf.gz"
    params:
        gvcf_string = expand("--variant " + join(RUN_DIR, "called/{sample}.g.vcf.gz"), sample=BAM_PREFIX)
    log:
        "logs/gatk/combinegvcfs.log"
    shell: """
        gatk CombineGVCFs \
            --reference {input.ref} \
            {params.gvcf_string} \
            --output {output} \
    """

###############################################################################
rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf=rules.combine_calls.output
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
        vcf=rules.genotype_variants.output.vcf
    output:
        vcf="genotyped/all_merged.zarr"
    log:
        "logs/zarr/merged-zarr.log"
    script:
        "rules/zarr_format.py {input} {output}"
