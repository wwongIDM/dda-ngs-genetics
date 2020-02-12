'''
Purpose: Pulls out singleton and doubleton biallelic SNPs found in a specified subset of samples.
Author: Jessica Ribado
Date: 02/2020

Example run (parrallelized for different populations):
 PF3K_DIR='/mnt/md0/malariaGen_genomic/pf3k_v5/'
 parallel python 202001_zarr2singleton.py --zarr_path $PF3K_DIR/pf3k_zarr --metafile $PF3K_DIR/pf3k_release_5_metadata_202002_regionAdd.txt --out_path $PF3K_DIR/singleton_positions/ --chrom {1}  --group {2} --grouping_var population ::: {01..14} ::: $(sed 1d /mnt/md0/malariaGen_genomic/pf3k_v5/pf3k_release_5_metadata_202002_regionAdd.txt | cut -f 7 | uniq | grep -v "LAB\|PARENT")
'''


################################################################################
from os.path import join
import zarr
import numpy as np
import pandas as pd
import allel
import dask.array as da
import argparse
import zarr_functions as zf

################################################################################
parser = argparse.ArgumentParser(description='Store variables for parsing singletons.')
parser.add_argument('--zarr_path', action="store")
parser.add_argument('--metafile', action="store", help='Tab delimted file with sample names and grouping variable.')
parser.add_argument('--out_path', action="store")
parser.add_argument('--chrom', action="store")
parser.add_argument('--grouping_var', action="store", help='Variable in metafile that specifies the sames grouped for singleton counts.')
parser.add_argument('--group', action="store")
args = parser.parse_args()

################################################################################
# read in zarr format
zarr_chr_path = join(args.zarr_path, "SNP_INDEL_Pf3D7_" + args.chrom + "_v3.zarr")
callset = zarr.open(zarr_chr_path, mode="r")
print("Loaded zarr file for chromosome " + args.chrom)

# read in metadata
meta = pd.read_csv(args.metafile, sep='\t', header=0)
meta.rename(columns={args.grouping_var:'group'}, inplace=True)

# subset the variants and samples pertaining to specified group(s)
hq_variants = zf.variant_filter(callset, 6)
group_samps = zf.group_index_list(callset, meta, args.group)
print(' '.join(["Variant filtering retained", str(np.count_nonzero(hq_variants)), " high quality variants out of a possible", str(len(hq_variants)), "variants."]))
print(' '.join(["Sample filtering by group", args.group, "retained", str(np.count_nonzero(group_samps)), "out of a possible", str(len(callset['samples'])), "samples."]))

# subset the genotype file for the different variants
gt_zarr = callset['calldata/GT']
gt_dask = allel.GenotypeDaskArray(gt_zarr)
gt_daskSub = gt_dask.subset(hq_variants, group_samps).compute()

# count singletons and doubletons
singleton_df = zf.singleton_sample_pairing(args.group, args.chrom,
    gt_daskSub,
    zf.singleton_count(gt_daskSub.count_alleles()),
    callset['variants/POS'][:][hq_variants], callset['samples'][:][meta[meta.group == args.group].callset_index.tolist()])

# print to file
singleton_df.to_csv(join(args.out_path + '_'.join([args.group, args.chrom, "singletons.txt"])), sep="\t", index=False)
print(' '.join(["Singleton results for group", args.group, "on chromosome", args.chrom, "are saved in file"]), join(args.out_path + '_'.join([args.group, args.chrom, "singletons.txt"])))
