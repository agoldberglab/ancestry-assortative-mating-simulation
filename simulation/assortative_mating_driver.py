#! /usr/bin/python3

from argparse import ArgumentParser, ArgumentTypeError
from functools import partial
import numpy as np
import os
import pandas as pd
import subprocess
import sys
from utils import file_type, format_chroms, calculate_summary_stats, append_labels
from slim import run_slim
from global_ancestry import (bin_global_ancestry, bin_offspring_by_ancestry, load_pedigree,
        permute_pedigree, summarize_global_ancestry, update_pedigree)
from local_ancestry import (bin_tract_lengths, calculate_expected_length_dist,
        estimate_timing, load_tracts)

def calculate_mating_weight(mate_choice, theta, args_dict):

    def calculate_exponential_alpha(theta, normalized=False):
        alpha = np.log(theta)
        if normalized:
            alpha = alpha / 2
        return alpha

    def calculate_normal_dist_sigma(theta):
        if theta == 1:
            sigma = 1000
        else:
            denom = 2 * np.log(1/theta)
            sigma = np.sqrt(-1/denom)
        return sigma

    f = {'exponential_decay': partial(calculate_exponential_alpha,
                                      normalized=False),
         'exponential_decay_normalized':
                             partial(calculate_exponential_alpha,
                                                 normalized=True),
         'normal_distribution': calculate_normal_dist_sigma,
         'subpop' : lambda x: x
        }

    keys = {'exponential_decay': 'alpha',
            'exponential_decay_normalized': 'alpha',
            'normal_distribution': 'sigma',
            'subpop' : 'theta'
           }

    if args_dict['mate_choice'] not in keys:
        raise ValueError(("%s is not a valid mate-choice function."
                %(args_dict['mate_choice'])))
    args_dict.update({keys[mate_choice]: f[mate_choice](theta)})
    return args_dict

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('--source', required=True,
                        type=file_type('slim'))
    parser.add_argument('--seed', required=True, type=int)
    parser.add_argument('--theta', required=True, type=float)
    parser.add_argument('--mate_choice',  required=True)
    parser.add_argument('--recombination_rate_file',
                         required=True, type=file_type)
    parser.add_argument('--initial_proportion', default=0.5,
                        type=float)
    parser.add_argument('--migration_rate', default=0,
                        type=float)
    parser.add_argument('--population_size', default=10000,
                        type=int)
    parser.add_argument('--n_generations', default=20, type=int)
    parser.add_argument('--outdir', required=True)
    
    args = vars(parser.parse_args())

    args = calculate_mating_weight(args['mate_choice'],
                                   args['theta'], args)

    return args

def analyze_global_ancestry(args, generation):
    current_gen = load_pedigree(args['outdir'], generation)
    next_gen = load_pedigree(args['outdir'], generation+1)
    previous_gen = load_pedigree(args['outdir'], generation-1)
    current_gen = update_pedigree(current_gen, next_gen, args)
    binned_ancestry0 = bin_global_ancestry(current_gen)
    permuted = permute_pedigree(previous_gen)
    summary, shuffled_rho, shuffled_delta = summarize_global_ancestry(current_gen, next_gen, permuted)
    binned_offspring = bin_offspring_by_ancestry(previous_gen, current_gen, permuted)
    return current_gen, binned_ancestry0, summary, shuffled_rho, shuffled_delta, binned_offspring

def analyze_ancestry_tracts(args, generation, chroms):
    tracts = load_tracts(args['outdir'], generation)
    tract_summary, tract_outliers = calculate_summary_stats(tracts['length_cM'] * 100, 'tracts')
    tract_summary['tracts_q95'] = np.quantile(tracts[tracts['source'] == 0]['length_cM']*100, 0.95)
    estimated_timing, ld_decay = estimate_timing(tracts, chroms, args['seed'],
            generation, n=10000)
    q, exp_outliers = calculate_expected_length_dist(tracts,estimated_timing)
    tract_summary = pd.concat([tract_summary, q, estimated_timing], axis=1)
    return tract_summary, tract_outliers, exp_outliers, ld_decay

def analyze_generation(args, generation, chroms):
    (pedigree, binned_ancestry0, summary, shuffled_rho,
            shuffled_delta, binned_offspring) = analyze_global_ancestry(args, generation)
    tract_summary, tract_outliers, exp_outliers, ld_decay = analyze_ancestry_tracts(args,
            generation, chroms)
    summary = pd.concat([summary, tract_summary],axis=1)
    if generation != 20:
        tract_outliers = None
        exp_outliers = None
    output = [pedigree, binned_ancestry0, summary, shuffled_rho,
            shuffled_delta, binned_offspring, tract_outliers, exp_outliers, ld_decay]
    return output

def initialize_output(tables):
    output = {}
    for t in range(0, len(tables)):
        output[tables[t]] = {}
    return output

def append_output(output, data, tables, generation, args):
    for t in range(0, len(tables)):
        if data[t] is not None:
            output[tables[t]] = append_labels(data[t], generation, args)
    return output

def write_output(output, tables, generation, args):
    for t in range(0, len(tables)):
        table = tables[t]
        df = output[table]
        if not isinstance(df, pd.DataFrame):
            continue
        filename=args['outdir']+'/'+table+'.csv'
        header = False if os.path.isfile(filename) else True
        df.to_csv(args['outdir']+'/'+table+'.csv', mode='a',
                header=header, index=False, na_rep='NA')

def clean_up(outdir):
    cmd = 'rm ' + outdir + '/*.trees'
    subprocess.call(cmd, shell=True)
    cmd = 'rm ' + outdir + '/generation*pedigree.csv'
    subprocess.call(cmd, shell=True)
    cmd = 'rm ' + outdir + '/*local_ancestry_tracts.csv'
    subprocess.call(cmd, shell=True)
    cmd = 'rm ' + outdir + '/model.slim'
    subprocess.call(cmd, shell=True)

def main():
    args = parse_arguments()
    err = run_slim(args)
    if err != '':
        print(err)
        sys.exit(1)

    _,chroms,_ = format_chroms(args['recombination_rate_file'])
    tables = ['pedigree', 'binned_ancestry0', 'summary', 'shuffled_rho',
            'shuffled_delta', 'binned_offspring', 'tract_outliers', 'exp_outliers',
            'ld_decay']
    for generation in range(1, args['n_generations']+1):
        output = initialize_output(tables)
        data = analyze_generation(args, generation, chroms)
        output = append_output(output, data, tables, generation, args)
        write_output(output, tables, generation, args)

    clean_up(args['outdir'])

if __name__ == '__main__':
    main()
