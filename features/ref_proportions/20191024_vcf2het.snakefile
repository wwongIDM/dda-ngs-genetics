import os

'''
The goal of this analysis is to replicate earlier work by Josh Protoctor on the Pf3k v3 release on the heterozygosity at variant sites. A higher proportion of heterozygozity scores across the genome within a sample is indicative of polyinfections.
'''

# set directories
PROJECT_DIR="/mnt/c/Users/jribado/Desktop/pf_v6"
VCF_DIR=os.path.join(PROJECT_DIR, "ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf6/Pf_6_vcf")
# OUTPUT_DIR="/mnt/c/Users/jribado/Desktop/pf_v6/vcf2het"
OUTPUT_DIR="/mnt/c/Users/jribado/Dropbox\ \(IDM\)/Data\,\ Dynamics\,\ and\ Analytics\ Folder/Projects/Pfv6\ Malaria/vcf2het"
CHR_NUMS=["{:02d}".format(item) for item in list(range(1,15))]


rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "Pf3D7_{chr}_biPass6.vcf.gz"), chr=CHR_NUMS)


################################################################################
rule filter_vcf:
    input: os.path.join(VCF_DIR, "Pf_60_public_Pf3D7_{chr}_v3.final.vcf.gz")
    output: os.path.join(OUTPUT_DIR, "Pf3D7_{chr}_biPass6.vcf.gz")
    params:
        n_alt=config["alt_allele_max"],
        vsqlod=config["vsqlod_min"]
    shell:'''
        bcftools view --threads 32 \
            --include 'FILTER="PASS" && N_ALT={params.n_alt} && CDS==1 && TYPE="snp" && VQSLOD>{params.vsqlod}' \
            --output-type z \
            --output-file {output} {input}
    '''

# bcftools view \
# --include 'FILTER="PASS" && N_ALT=1 && CDS==1 && TYPE="snp" && VQSLOD>6.0' \
# --output-type z \
# --output-file Pf3D7_05_biPass6.vcf.gz \
# Pf_60_public_Pf3D7_05_v3.final.vcf.gz


################################################################################
rule vcf_melt:
    '''Convert vcf into long format for sample level parsing using the PyVCF vcf_melt script. Had to update vcf_melt scipt to turn key into lists at the beginning to avoid odict_keys errors in Python3 from Python2. '''
    input: rules.filter_vcf.output
    output: os.path.join(OUTPUT_DIR, "Pf3D7_{chr}_biPass6_melt.txt")
    shell:'''
        zcat {input} | vcf_melt > {output}
    '''

# zcat Pf3D7_05_biPass6.vcf.gz | vcf_melt > Pf3D7_05_biPass6_melt.txt

################################################################################
rule calculate_ref_mapping:
    '''Calculate the percentage of reads that map to the reference allele. Remove variant sites where all samples have the alternative allele. '''
    input: rules.vcf_melt.output
    output:
        raw_counts  = os.path.join(OUTPUT_DIR, "snpHetCounts_chr{chr}_raw.txt"),
        uniq_counts = os.path.join(OUTPUT_DIR, "snpHetCounts_chr{chr}_uniq.txt"),
        end_counts  = os.path.join(OUTPUT_DIR, "snpHetCounts_chr{chr}.txt")
    params:
        min_counts = config["min_site_reads"]
    shell:"""
        cat {input} | cut -f 1,2,3,4,8,9 | tr ',' '\t' | awk \'FNR != 1 && $4 > {params.min_counts} {{print $0, $2/$4}}\' > {output.raw_counts}
        cat {output.raw_counts} | awk \'0 < $8 && $8 < 1 {{print $0}}\' | cut -f 7 | uniq > {output.uniq_counts};
        grep -F -f {output.uniq_counts} {output.raw_counts} > {output.end_counts}
    """

        #| tr \",\" \"\t\"

        #awk "{{OFS='\\t'}}; FNR != 1 && $4 > 10 {{print $0, $2/$4}}" > {{output.raw_counts}};
        #cat {output.raw_counts} | awk "0 < $8 && $8 < 1 {{print $0}}" | cut -f 7 | uniq > {output.uniq_counts};
        # grep -wFf  {output.uniq_counts} {output.raw_counts} > {output.end_counts}
        # awk 'FNR==NR{{a[$1,$2]=$0;next}}{if(b=a[$6,$7]){{print $0}}' {output.uniq_counts} {output.raw_counts} > {output.end_counts}

# cat Pf3D7_05_biPass6_melt.txt | cut -f 1,2,3,4,8,9 | tr "," "\t" | awk -v OFS='\t' 'FNR != 1 && $4 > 10 {print $0, $2/$4}' > snpHetCounts_chr5.txt
# cat snpHetCounts_chr5.txt | awk '0 < $8 && $8 < 1 {print $0}' | cut -f 7 | uniq  | wc -l
# awk 'FNR==NR{{a[$1,$2]=$0;next}}{if(b=a[$6,$7]){{print $0}}' {output.uniq_counts} {output.raw_counts} > {output.end_counts}
