#! /usr/bin/python3

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from utils import bin_data_continuous

def load_tracts(outdir, generation):
    basename = outdir + '/generation' + str(generation)
    try:
        tracts = pd.read_csv(basename + '_local_ancestry_tracts.csv', na_values='NA')
        tracts[['start_cM', 'end_cM', 'length_cM']] = tracts[['start', 'end', 'length']] * 1e-8
        tracts = tracts[~tracts['is_whole_chrom']]
        return tracts
    except:
        return None

def bin_tract_lengths(tracts, width=0.01):
    binned_tracts = {}
    for s in range(0, 2):
        source_tracts = tracts[tracts['source'] == s]
        binned_tracts[s] = bin_data_continuous(source_tracts['length_cM'],
                                               bin_size=width, max_value=2.49)
        binned_tracts[s]['source'] = s
        binned_tracts[s]['ndensity'] = (binned_tracts[s]['count'] /
                                         max(binned_tracts[s]['count']))
    binned_tracts = pd.concat(binned_tracts, axis=0).reset_index(drop=True)
    binned_tracts = binned_tracts.drop('count', axis=1)
    binned_tracts = binned_tracts[['source', 'length_cM', 'ndensity']]
    return  binned_tracts

def estimate_timing_from_tract_length(tracts, width=0.01):
    binned_tracts = bin_tract_lengths(tracts, width)
    n = sum(binned_tracts['source']==0)
    def f(x, t1, t2, m):
        def f_source0(x, t1, t2, m):
            raw = m * t1 * np.exp(-m * t1 * x)
            return raw/max(raw)
        def f_source1(x, t1, t2, m):
            raw = (1-m) * t2 * np.exp(-(1-m) * t2 * x)
            return raw/max(raw)
        return np.append(f_source0(x[:n], t1, t2, m),
                         f_source1(x[n:], t1, t2, m))
    if tracts.empty:
        p = [np.nan, np.nan, np.nan]
    else:
        p, _ = curve_fit(f, binned_tracts['length_cM'], binned_tracts['ndensity'],
                bounds=(0,(np.inf, np.inf, 1)))
    estimated_timing = pd.DataFrame({'timing_0':p[0]-1, 'timing_1':p[1]-1,
                                     'contribution_0':p[2],
                                     'contribution_1':1-p[2]}, index=[0])
    return estimated_timing

def define_snps(chroms, n, seed):
    rng = np.random.default_rng(seed)
    snps = rng.integers(chroms[-1], size=n)
    snp_chrom = pd.cut(snps, [c for c in chroms],
            include_lowest=True, right=True,
            labels=np.arange(1,len(chroms)))
    snps = pd.DataFrame({'pos': snps * 1e-8})
    return snps, snp_chrom

def count_alleles(tracts, snp, ids):
    t = tracts[(tracts['start_cM'] <= snp) & (tracts['end_cM'] >= snp)]
    allele_count = t.groupby('id')['source'].sum()
    allele_count = allele_count.reindex(ids).fillna(0)
    allele_count = 2 - allele_count
    allele_count = allele_count.to_frame().T.reset_index(drop=True)
    return allele_count

def calculate_ld(count, snps):
    ld = pd.DataFrame(np.corrcoef(np.concatenate(count)), columns=snps.pos, index=snps.pos)
    ld = ld.where(np.tril(np.ones(ld.shape)).astype(bool))
    ld.index.name = 'snp1'
    ld.columns.name = 'snp2'
    ld = ld.stack().reset_index()
    ld['dist'] = abs(ld['snp1'] - ld['snp2'])
    ld['rho'] = ld[0]
    ld = ld[['dist', 'rho']]
    return ld

def bin_ld(ld, bin_size, max_value):
    ld = ld[ld['dist'] > bin_size].copy()
    bins = np.arange(0, max_value+bin_size, bin_size)
    ld['bin'] = pd.cut(ld['dist'], bins, include_lowest=True, right=True)
    ld = ld.groupby('bin').mean('rho').sort_index().reset_index()
    ld['dist'] = pd.IntervalIndex(ld['bin']).mid.values
    ld = ld.drop('bin', axis=1)
    ld = ld[~np.isnan(ld['rho'])]
    return ld

def calculate_admixture_LD(tracts, chroms, n, seed, bin_size, max_value):
    if tracts.empty:
        return None
    snps, snp_chrom = define_snps(chroms, n, seed)
    ids = np.unique(tracts['id'])
    ld = []
    for chrom in range(1, 23):
        t = tracts[tracts['chrom'] == chrom]
        s = snps[snp_chrom == chrom]
        count = s['pos'].apply(lambda x: count_alleles(t, x, ids)).reset_index(drop=True)
        ld.append(calculate_ld(count, s))
    ld = pd.concat(ld)
    ld_decay = bin_ld(ld, bin_size, max_value)
    return ld_decay

def estimate_timing_from_admixture_LD(ld_decay):
    if ld_decay is None:
        return np.nan
    def f(x, d):
        return np.exp(-d * x)
    p,_ = curve_fit(f, ld_decay['dist'], ld_decay['rho'])
    return p[0]

def estimate_timing(tracts, chroms, seed, generation, n=10000):
    estimated_timing = estimate_timing_from_tract_length(tracts)
    ld_decay = calculate_admixture_LD(tracts, chroms, n=n,
            seed=seed, bin_size=0.005, max_value=28.8)
    estimated_timing['g_LD'] = estimate_timing_from_admixture_LD(ld_decay)
    estimated_timing['timing_discrepancy'] = generation - estimated_timing['timing_0']
    return estimated_timing, ld_decay

def calculate_expected_quantile(estimated_timing, quantile):
    l = ((estimated_timing.loc[0,'timing_0']+1) *
            estimated_timing.loc[0,'contribution_0'])
    if not isinstance(quantile, list):
        quantile = [quantile]
    q = [-np.log(1-x)/l if x > 0 else 0 for x in quantile]
    return q

def calculate_expected_length_dist(tracts, estimated_timing):
    tracts = tracts[tracts['source'] == 0].reset_index(drop=True)
    l_rate = ((estimated_timing.loc[0,'timing_0']+1) *
            estimated_timing.loc[0,'contribution_0'])
    q = calculate_expected_quantile(estimated_timing, [0, 0.25, 0.5, 0.75, 0.95])
    iqr = np.log(3)/l_rate
    q.append(q[3] + (1.5 * iqr))
    q.append(1/l_rate)
    outliers =  tracts[tracts['length_cM'] > q[5]]['length_cM']
    outliers = outliers.reset_index(drop=True)
    try:
        p_outliers = len(outliers.index) / len(tracts.index)
    except:
        p_outliers = np.nan
    q.append(p_outliers)
    q = pd.DataFrame([q], columns = ['tracts_exp_min', 'tracts_exp_q25',
        'tracts_exp_q50', 'tracts_exp_q75', 'tracts_exp_q95', 'tracts_exp_max',
        'tracts_exp_mean', 'p_outliers'])
    return q*100, outliers*100
