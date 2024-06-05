#! /usr/bin/python3

import numpy as np
import pandas as pd
from scipy import stats
from utils import append_labels, bin_data_continuous, calculate_summary_stats

def load_pedigree(outdir, generation):
    basename = outdir + '/generation' + str(generation)
    try:
        pedigree = pd.read_csv(basename + '_pedigree.csv', na_values='NA')
        return pedigree
    except:
        return None

def calculate_correlation_ancestry(pedigree):
    rho = pedigree['parent1_ancestry0'].corr(
            pedigree['parent2_ancestry0'])
    return rho

def calculate_delta_ancestry(pedigree):
    delta = np.mean(np.abs(pedigree['parent1_ancestry0'] -
        pedigree['parent2_ancestry0']))
    return delta

def count_migrant_pairs(pedigree, next_gen):
    if next_gen is None:
        return 0.0
    migrants = pedigree[np.isnan(pedigree['parent1'])]['id'].tolist()
    n_migrant_pairs = sum(next_gen['parent1'].isin(migrants) &
            next_gen['parent2'].isin(migrants))
    p_migrant_pairs = n_migrant_pairs / len(pedigree.index) * 100
    return p_migrant_pairs

def count_offspring(pedigree, next_gen):
    if next_gen is None:
        pedigree['offspring1'] = np.nan
        pedigree['offspring2'] = np.nan
        return pedigree
    for p in [1, 2]:
        offspring = next_gen['parent'+str(p)].value_counts().to_dict()
        pedigree['offspring'+str(p)] = pedigree['id'].map(offspring).fillna(0).astype('int')
    return pedigree

def calculate_psi_ratio(pedigree, ancestry0, args):
    def calculate_psi_max(args):
        if args['mate_choice'] == 'normal_distribution':
            psi_max = stats.norm.pdf(0, 0, args['sigma'])
        else:
            psi_max = 1
        return psi_max
    def calculate_psi_min(pedigree, ancestry0, args):
        delta = max(np.abs(pedigree['ancestry0'] - ancestry0))
        if args['mate_choice'] == 'exponential_decay':
            psi_min = np.exp(-args['alpha'] * delta)
        elif args['mate_choice'] == 'exponential_decay_normalized':
            sigma = pedigree['ancestry0'].std()
            psi_min = np.exp(-args['alpha']/sigma * delta)
        elif args['mate_choice'] == 'normal_distribution':
            psi_min = min(stats.norm.pdf(pedigree['ancestry0'], ancestry0, args['sigma']))
        else:
            psi_min = np.nan
        return psi_min
    psi_max = calculate_psi_max(args)
    psi_min = calculate_psi_min(pedigree, ancestry0, args)
    return psi_max/psi_min

def estimate_timing_from_variance(pedigree):
    mu = np.mean(pedigree['ancestry0'])
    sigma = np.var(pedigree['ancestry0'])
    estimated_timing = (np.log(mu * (1-mu))-np.log(sigma)) / np.log(2)
    return estimated_timing

def permute_mating_pairs(pedigree):
    def permute_column(pedigree, name, n=len(pedigree.index)):
        pedigree = pedigree.sample(n=n, replace=True).reset_index(drop=True)
        pedigree = pedigree.rename(columns={'id':name, 'ancestry0':name+'_ancestry0'})
        return pedigree
    pedigree = pedigree[['id', 'ancestry0']]
    permuted = pd.merge(permute_column(pedigree, 'parent1'), permute_column(pedigree, 'parent2'),
            left_index=True, right_index=True)
    selfing = permuted[permuted['parent1'] == permuted['parent2']].index
    for s in selfing:
        is_selfing = True
        while is_selfing:
            new_parent = permute_column(pedigree, 'parent2', n=1)
            is_selfing = new_parent.loc[0]['parent2'] == permuted.loc[s]['parent1']
        permuted.loc[s, 'parent2'] = new_parent.loc[0, 'parent2']
        permuted.loc[s, 'parent2_ancestry0'] = new_parent.loc[0, 'parent2_ancestry0']
    return permuted

def permute_pedigree(pedigree, n_permutations=1000):
    if pedigree is None:
        return None
    permuted = []
    for i in range(0, n_permutations):
        permuted.append(permute_mating_pairs(pedigree))
    return permuted

def calculate_expressed_bias(pedigree, permuted):
    if permuted is None:
        return np.nan, np.nan, []
    delta = calculate_delta_ancestry(pedigree)
    n_permutations = len(permuted)
    shuffled_delta = []
    for p in range(0, n_permutations):
        shuffled_delta.append(calculate_delta_ancestry(permuted[p]))
    B = 1 - (delta / np.mean(shuffled_delta))
    pval = (sum(shuffled_delta <= delta) + 1) / (n_permutations + 1)
    shuffled_delta = pd.DataFrame({'delta_perm': shuffled_delta})
    return B, pval, shuffled_delta

def calculate_correlation_bias(pedigree, permuted):
    if permuted is None:
        return np.nan, np.nan, []
    rho = calculate_correlation_ancestry(pedigree)
    n_permutations = len(permuted)
    shuffled_rho = []
    for p in range(0, n_permutations):
        shuffled_rho.append(calculate_correlation_ancestry(permuted[p]))
    pval = (sum(shuffled_rho >= rho) + 1) / (n_permutations + 1)
    shuffled_rho = pd.DataFrame({'rho_perm': shuffled_rho})
    return pval, shuffled_rho

def update_pedigree(pedigree, next_gen, args):
    pedigree = count_offspring(pedigree, next_gen)
    pedigree['psi_ratio'] = pedigree.apply(lambda x:
            calculate_psi_ratio(pedigree, x['ancestry0'], args), axis=1)
    return pedigree

def summarize_global_ancestry(pedigree, next_gen, permuted):
    summary,_ = calculate_summary_stats(pedigree['ancestry0'], 'ancestry0')
    rho = calculate_correlation_ancestry(pedigree)
    delta_obs = calculate_delta_ancestry(pedigree)
    p_migrant_pairs = count_migrant_pairs(pedigree, next_gen)
    B, pval_B, shuffled_delta = calculate_expressed_bias(pedigree, permuted)
    pval_rho, shuffled_rho = calculate_correlation_bias(pedigree, permuted)
    g_var = estimate_timing_from_variance(pedigree) 
    subpop_counts = pedigree.groupby('subpop').agg('count')['ancestry0'].values
    subpop_ancestry0 = pedigree.groupby('subpop').agg('mean')['ancestry0'].values
    summary_extended = pd.DataFrame({'parent_corr':rho,
        'count_0': subpop_counts[0], 'count_1':subpop_counts[1],
        'mean_0':subpop_ancestry0[0], 'mean_1':subpop_ancestry0[1],
        'delta_obs': delta_obs, 'B': B, 'pval_B': pval_B,
        'pval_rho': pval_rho, 'g_var': g_var, 'p_migrant_pairs': p_migrant_pairs},
        index = [0])
    summary = summary.join(summary_extended)
    return summary, shuffled_rho, shuffled_delta

def bin_global_ancestry(pedigree):
    groups = [[0], [1], [0, 1]]
    binned_ancestry = []
    for g in groups:
        p = pedigree[np.isin(pedigree['subpop'], g)]
        b = bin_data_continuous(p['ancestry0'], bin_size=0.02, max_value=1)
        b['proportion'] = (b['count'] / sum(b['count']))
        b['subpop'] = str(g[0]) if len(g) == 1 else 'all'
        b.drop('count', axis=1, inplace=True)
        binned_ancestry.append(b)
    binned_ancestry = pd.concat(binned_ancestry)
    return binned_ancestry

def bin_offspring(pedigree, next_gen):
    pedigree = count_offspring(pedigree, next_gen)
    bins = np.round(np.arange(0, 1.01, 0.01), decimals=2)
    bin_centers = [np.mean([x,y]) for x,y in zip(bins,bins[1:])]
    binned = pd.DataFrame({'bin': bin_centers})
    binned = binned.set_index('bin')
    for parent in [1, 2]:
        col = 'parent' + str(parent) + '_ancestry0'
        pedigree['bin'] = pd.cut(pedigree[col], bins, include_lowest=True, right=True)
        pedigree['bin'] = pedigree['bin'].apply(lambda x: x.mid)
        pedigree['bin'] = pedigree['bin'].astype('category')
        pedigree['bin'] = pedigree['bin'].cat.set_categories(bin_centers)
        col = 'offspring'+str(parent)
        counts = pedigree.groupby('bin').sum()[col].rename(col)
        binned[col] = counts / sum(counts)
    return binned

def calculate_binned_offspring_controls(binned, pedigree, permuted):
    b = []
    for p in range(0, len(permuted)):
        b.append(bin_offspring(pedigree, permuted[p]))
    b = pd.concat(b, axis=1)
    binned = pd.concat([binned, pd.DataFrame({'permuted_min': b.min(axis=1),
        'permuted_max': b.max(axis=1)})],axis=1)
    binned = binned.reset_index()
    binned = binned.rename(columns={'bin': 'ancestry0'})
    return binned

def bin_offspring_by_ancestry(pedigree, next_gen, permuted):
    if 'parent1_ancestry0' not in pedigree.columns:
        return None
    binned = bin_offspring(pedigree, next_gen)
    binned = calculate_binned_offspring_controls(binned, pedigree, permuted)
    return binned
