
################################################################################
# import packages
from snakemake.utils import validate
from snakemake.utils import min_version
import zarr
import allel

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
