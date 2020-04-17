################################################################################
# import scripts
configfile: "config.yaml"
include: "rules/calling.smk"
include: "rules/formats.smk"

################################################################################
# specify rule lists
BAM_PREFIX =  path(config['bam_dir']).glob('**_recal.bam')

################################################################################
rule all:
    input: "genotyped/all.zarr"
