from os.path import join
import zarr
import allel

VCF_FILE, OUT_FILE  = sys.argv[:]
allel.vcf_to_zarr(VCF_FILE, OUT_FILE, fields='*', overwrite=False)
