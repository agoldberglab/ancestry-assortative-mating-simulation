#! /usr/bin/python3

import numpy as np
import os
import pandas as pd
import re
import tskit
from utils import file_type, format_chroms

def format_pedigree(ts):
    pedigree = pd.DataFrame(i for i in ts.individuals()).rename(
      columns={'id':'index'})
    pedigree['subpop'] = pd.DataFrame(
      pedigree['location'].to_list())[0].astype('int')
    pedigree[['genome1', 'genome2']] = pd.DataFrame(
      pedigree['nodes'].to_list())
    
    metadata = (pd.DataFrame(pedigree['metadata'].to_list(),
      columns=['pedigree_id', 'pedigree_p1', 'pedigree_p2', 'subpopulation']).
      rename(columns={'pedigree_id':'id', 'pedigree_p1':'parent1',
                      'pedigree_p2':'parent2'}).reset_index())
    
    pedigree = pedigree[['index', 'subpop', 'genome1', 'genome2']].merge(
      metadata, on='index')
    pedigree = pedigree[pedigree['subpopulation'] == 3].drop(
      ['index', 'subpopulation'], axis=1).copy()

    return pedigree


def record_parent_ancestry(trees, generation, pedigree):

    g = generation-1
    if g < 0:
        return pedigree

    previous = pd.read_csv(os.path.dirname(trees) + '/generation' +
                           str(g) + '_pedigree.csv', na_values='NA')
    pedigree = pedigree.merge(previous[['id', 'ancestry0']].set_index('id'),
                              left_on='parent1', right_on='id',
                              how='left', validate='many_to_one').rename(
                              columns={'ancestry0':'parent1_ancestry0'})
    pedigree = pedigree.merge(previous[['id', 'ancestry0']].set_index('id'),
                              left_on='parent2', right_on='id',
                              how='left', validate='many_to_one').rename(
                              columns={'ancestry0':'parent2_ancestry0'})
    index = np.isnan(pedigree['parent1_ancestry0'])
    pedigree.loc[index, 'parent1'] = np.nan
    pedigree.loc[index, 'parent2'] = np.nan
    return pedigree


def link_ancestors(ts, nodes, ancestors):
    edges = ts.tables.link_ancestors(nodes, ancestors)
    edges = [row for row in [edges.child, edges.parent,
       edges.left, edges.right]]
    edges = pd.DataFrame({'node': edges[0], 'source': edges[1],
                         'start' : edges[2], 'end': edges[3]}).astype('int')
    return edges

def format_local_ancestry_tracts(ts, batch, pedigree):
    tracts = link_ancestors(ts, batch, [0, 1, 2, 3])
    tracts.replace({'source': {0:0, 1:0, 2:1, 3:1}}, inplace=True)
    tracts.sort_values(['node','start'], inplace=True)
    ids = pedigree.apply(lambda x: {x.genome1:x.id, x.genome2:x.id}, axis=1)
    ids = {k:v for d in ids for k,v in d.items()}
    tracts['id'] = tracts['node'].map(ids)    
    phasing = pedigree.apply(lambda x: {x.genome1:0, x.genome2:1}, axis=1)
    phasing = {k:v for d in phasing for k,v in d.items()}
    tracts['phase'] = tracts['node'].map(phasing)
    tracts = tracts.astype('int')
    return tracts

def validate_chromosomes(tracts, chroms):
    labels = np.arange(1,len(chroms))
    chrom_left = pd.cut(tracts['start'], chroms, right=False, labels=labels)
    chrom_right = pd.cut(tracts['end'], chroms, right=True, labels=labels)
    assert(all(chrom_left == chrom_right))
    tracts['chrom'] = chrom_left.astype('int')
    return tracts

def split_by_chromosome(start, end, chroms):
    breaks = np.sort(np.unique([start]+[end]+chroms))
    breaks = breaks[(breaks >= start) & (breaks <= end)]
    tracts = [[x, y] for x,y in zip(breaks, breaks[1:])]
    return tracts

def merge_tracts_by_ancestry(tracts):
    tracts_shifted = tracts.shift().fillna(method='backfill')
    g = ((tracts['node'] == tracts_shifted['node']) &
         (tracts['source'] != tracts_shifted['source']))
    g = g.cumsum().rename('group')
    tracts = tracts.groupby(['id', 'phase', 'node', 'source', g],
            observed=True).agg({'start':'min',
                'end':'max'}).reset_index()
    tracts = tracts.drop('group', axis = 1)
    return tracts

def assign_chromosomes(tracts, chroms, chrom_sizes):
    tracts = validate_chromosomes(tracts, chroms)
    tracts['chrom_size'] = (tracts['chrom']).map(chrom_sizes)
    tracts['is_whole_chrom'] = tracts['length'] == tracts['chrom_size']
    tracts = tracts.drop('chrom_size', axis=1)
    return tracts

def call_tract_lengths(ts, b, pedigree, chroms, chrom_sizes):
    tracts = format_local_ancestry_tracts(ts, b, pedigree)
    tracts = merge_tracts_by_ancestry(tracts)
    tracts['coords'] = tracts.apply(
                lambda row: split_by_chromosome(row['start'],
                row['end'], chroms), axis=1)
    tracts = tracts.explode('coords')
    tracts[['start', 'end']] = tracts['coords'].tolist()
    tracts = tracts.drop('coords', axis=1)
    tracts['length'] = tracts['end'] - tracts['start']
    total_length = tracts[['node', 'length']].groupby('node').sum()
    assert(all(tracts.groupby('node').sum() == ts.sequence_length))
    tracts = assign_chromosomes(tracts, chroms, chrom_sizes)
    tracts = tracts.drop(['node', 'phase'], axis=1)
    return tracts

def call_local_ancestry_tracts(ts, pedigree, n):
    _, chroms, chrom_sizes = format_chroms(recombination_map)
    nodes = np.concatenate(pedigree[['genome1', 'genome2']].to_numpy())
    tracts = call_tract_lengths(ts, nodes[0:6], pedigree, chroms, chrom_sizes)
    batch_nodes = [nodes[i * n:(i + 1) * n]
                   for i in range((len(nodes) + n - 1) // n )]
    tracts = [call_tract_lengths(ts, b, pedigree, chroms, chrom_sizes) for b in batch_nodes]
    tracts = pd.concat(tracts).reset_index(drop = True)
    ancestry = calculate_ancestry_proportion(ts,tracts)
    return tracts, ancestry

def calculate_ancestry_proportion(ts, tracts):
    ancestry = tracts.groupby(['id', 'source']).sum()['length'].reset_index()
    ancestry['proportion'] = ancestry['length']/(2*ts.sequence_length)
    assert(all(ancestry.groupby('id').sum()['proportion'] == 1))
    ancestry = ancestry[ancestry['source'] == 0]
    ancestry = ancestry[['id', 'proportion']]
    return ancestry

def merge_ancestry_into_pedigree(pedigree, ancestry):
    pedigree = pedigree.merge(ancestry[['id', 'proportion']].set_index('id'),
                              on='id', how='left')
    pedigree = pedigree.rename(columns={'proportion':'ancestry0'})
    pedigree = pedigree.fillna({'ancestry0':0})
    pedigree = pedigree.drop(['genome1', 'genome2'], axis=1)
    pedigree = pedigree.sort_values('id').reset_index(drop=True)
    pedigree[['id', 'parent1', 'parent2']] = pedigree[['id', 'parent1', 'parent2']].astype('Int64')
    return pedigree

def main():
    tractfile = (os.path.dirname(trees) + '/' + basename +
               '_local_ancestry_tracts.csv')
    ts = tskit.load(trees).simplify()
    pedigree = format_pedigree(ts)
    pedigree = record_parent_ancestry(trees, generation, pedigree)
    tracts, ancestry = call_local_ancestry_tracts(ts, pedigree, batch_n)
    pedigree = merge_ancestry_into_pedigree(pedigree, ancestry)
    pedigree.to_csv((os.path.dirname(trees) + '/' + basename +
               '_pedigree.csv'), index=False, na_rep='NA')
    tracts.to_csv(tractfile, index=False) 

if __name__ == '__main__':
    from argparse import ArgumentParser, ArgumentTypeError

    parser = ArgumentParser()
    parser.add_argument('--trees', required=True, type=file_type('trees'))
    parser.add_argument('--recombination', required=True,
                         type=file_type) 
    parser.add_argument('--migration', required=True,
                         type=float)
    args = vars(parser.parse_args())

    trees = args['trees']
    basename = os.path.basename(trees).split('.')[0]
    generation = int(re.search('generation(.*)', basename).group(1)) 

    recombination_map = args['recombination']

    migration = args['migration']
    batch_n = 1000 if migration == 0 else 36

    main()
