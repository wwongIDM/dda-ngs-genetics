from os.path import join
from pathlib import Path

################################################################################
# pull sample prefixes
BAM_FILES  = Path(config['bam_dir']).glob('*_recal.bam')
BAM_PREFIX = [os.path.split(f)[-1].split("_")[0] for f in BAM_FILES]

################################################################################
rule all:
	input:
		expand("called/{sample}.g.vcf.gz", sample=BAM_PREFIX),
		"genotyped/all_merged.zarr"

include: "rules/calling.smk"
