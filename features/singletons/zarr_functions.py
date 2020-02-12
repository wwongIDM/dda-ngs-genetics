import zarr
import numpy as np
import pandas as pd
import allel
import dask.array as da
import argparse
from collections import Counter


def group_index_list(zarr_object, metafile, group):
    ''' Check that the samples in the metafile are aligned with the samples in the Zarr file. '''
    samples_list = list(zarr_object['samples'])
    metafile['callset_index'] =  [samples_list.index(s) for s in metafile['sample']]
    return(metafile[metafile.group == group].callset_index.tolist())


def singleton_count(allele_count_dask):
    return(np.logical_or(allele_count_dask.is_singleton(allele=1)[:],
                  allele_count_dask.is_doubleton(allele=1)[:]))

def singleton_sample_pairing(group, chrom, gt_daskSub, singleton_bool, pos_filtSub, sample_popSub_index):
    ''' Creates a data frame with the position of each singleton, the number of samples that had a genotype or missing genotype for fruther filtering, and the samples that had that genotype. '''
    tmp_df = pd.DataFrame()
    for pos_index, x in enumerate(singleton_bool):
        sample_list = []
        if x  == True:
            # pull allele count from row
            test = gt_daskSub[pos_index, :, 1]
            # count the number of samples with seach genotype
            genotype_counts = dict(Counter(test))
            minor_genotype = [k for k, v in genotype_counts.items() if k > 0 & v <= 1]
            # pull the samples with the minor genotype
            for index, genotype in enumerate(test):
                if genotype == minor_genotype :
                    sample_list.append(sample_popSub_index[index])
                    tmp_df = tmp_df.append({
                            'pop': group,
                            'chr': chrom,
                            'pos': int(pos_index),
                            'n_major_genotype': genotype_counts.get(0),
                            'n_minor_genotype': genotype_counts.get(1),
                            'n_missing_genotype': genotype_counts.get(-1),
                            'samples': ';'.join(sample_list)
                        }, ignore_index=True)
    return(tmp_df)

def variant_filter(zarr_object, vsqlod):
    '''Filter variants based on MalariaGen guidelines. Currently returns biallelic SNPS. '''
    quality_set   = zarr_object['variants/FILTER_PASS'][:]
    snp_set       = zarr_object['variants/is_snp'][:]
    biallelic_set = zarr_object['variants/numalt'][:] == 1
    vsqlod_set    = zarr_object['variants/VQSLOD'][:]  > vsqlod
    return(quality_set & snp_set & vsqlod_set & biallelic_set)
