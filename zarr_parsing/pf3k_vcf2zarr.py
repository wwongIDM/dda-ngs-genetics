# run in parallel
################################################################################
# conda activate vcf
# VCF_DIR="/mnt/md0/malariaGen_genomic/pf3k_v5/5.1"
# OUT_DIR="/mnt/md0/malariaGen_genomic/pf3k_v5/pf3k_zarr""
# parallel "python3 pf3k_vcf2zarr.py $VCF_DIR $OUT_DIR {}":::{01..14}

################################################################################
import os, sys
import zarr
import allel

VCF_DIR, OUT_DIR, CHROM = sys.argv[:]
allel.vcf_to_zarr("/".join([VCF_DIR, "_".join(["SNP_INDEL_Pf3D7", CHROM, "v3.combined.filtered.vcf.gz"])]),
                       "/".join([OUT_DIR, "_".join(["SNP_INDEL_Pf3D7", CHROM, "v3.zarr"])]),
                       fields='*', overwrite=False)
