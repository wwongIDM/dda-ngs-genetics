# Next-generation sequencing processing pipelines

Pipelines are set-up to take different input files (raw fastq reads or BAM aligned files) that have been provided by collaborators. To fit within the scope of computational resources, variant calling with GATK 4.1.4 is done on user provided genomic ranges. 

# Directories:
- Fastq_processing
Takes in reads directly from a sequencer and aligns to a specified reference. Requires a known set of high confidence variants to run as is. Will merge replicates and create sample level BAM and gvcf files. Only files included in the user metadata file are considered for joint genotype calling.

- bam_processing
Takes aligned reads to call variants in user defined windows. Currently set up to run and joint genotype call on all files in a user defined directory. (Can be upgraded as above, if needed).

- zarr_parsing
Contains Juptyer notebooks for VCF files in the Zarr format to extract features on any computer regardless of resources. 

# Other pipelines:
- bamReadCounts: pulls directly from alignment files from a user defined set of samples and positions in a bed file format. Currently set up so that files of interest must be in the same directory. (Can be upgraded to pull from multiple directories, if needed.)
