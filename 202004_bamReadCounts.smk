from os.path import join, split

################################################################################
# pull sample prefixes
BAM_FILES  = [line.strip() for line in open(config['bam_files'], 'r')]
BAM_PREFIX = [split(f)[-1].split(".bam")[0] for f in BAM_FILES]

################################################################################
rule all:
	input:
		expand("bam_counts/{sample}_counts.txt", sample=BAM_PREFIX)

###############################################################################
rule bam_index:
    input:  join(config['bam_dir'], "{sample}.bam"),
    output: join(config['bam_dir'], "{sample}.bam.bai")
    shell: """
        samtools index {input}
    """

###############################################################################
rule bam_read_counts:
	input:
		ref = config['bam_calls']['reference_genome'],
		bam = join(config['bam_dir'], "{sample}.bam"),
		bam_index = rules.bam_index.output,
		bed = config['bam_calls']['interval_bed']
	output: "bam_counts/{sample}_counts.txt"
	params:
		min_qual = config['bam_calls']['min_read_qual']
	log:
		"logs/bam-readcounts/{sample}.log"
	shell: """
		bam-readcount -w 1 -q {params.min_qual} \
			-f {input.ref} {input.bam} \
			-l {input.bed} > {output}
	"""
