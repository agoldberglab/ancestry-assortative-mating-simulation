#! /usr/bin/python3

from argparse import ArgumentParser, ArgumentTypeError
import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from scipy import stats

class file_type(str):
    def __init__(self, extension=None):
        self.extension = extension
    
    def __call__(self,name):
        if not os.path.isfile(name):
            raise ArgumentTypeError('File "%s" was not found.' %(name))
    
        if not self.extension==None:
            if not name.endswith(self.extension):
                raise ArgumentTypeError('"%s" is not a .%s file.'
                                         %(name,self.extension))
        return os.path.abspath(name)

def format_chroms(recombination_map):
    recombination_map = pd.read_csv(recombination_map,
                         names=['position', 'rate'])
    chroms = list(recombination_map[recombination_map['rate'] == 0.5]
            ['position'])
    chroms = [0] + chroms + [max(recombination_map['position'])+1]

    chrom_sizes = pd.DataFrame({'chrom':np.arange(1,len(chroms)),
      'length':np.diff(chroms)}).set_index('chrom')['length'].to_dict()

    return recombination_map, chroms, chrom_sizes

def calculate_quantiles(series):
    q = series.quantile([0, 0.25, 0.5, 0.75, 1]).values
    iqr = q[3]-q[1]
    lower = q[1]-(1.5*iqr)
    upper = q[3]+(1.5*iqr)
    filtered = series.clip(lower=lower, upper=upper)
    q[0] = np.min(filtered)
    q[4] = np.max(filtered)
    outliers = series[(series < lower) | (series > upper)].reset_index(drop=True) 
    return q, outliers

def calculate_summary_stats(series, name):
    q, outliers = calculate_quantiles(series)
    if series.empty:
        series = pd.Series(np.nan)
        outliers = None
    summary = pd.DataFrame({name+'_mean':np.mean(series), name+'_var':np.var(series),
        name+'_min':q[0], name+'_q25':q[1], name+'_q50':q[2],
        name+'_q75':q[3], name+'_max':q[4],
        name+'_skew':stats.skew(series), name+'_kurtosis':stats.kurtosis(series)},
        index=[0])
    return summary, outliers

def calculate_violin_density(df, col):

    def df_to_R(df):
        with (ro.default_converter + ro.pandas2ri.converter).context():
            df_R = ro.conversion.get_conversion().py2rpy(df)
        return df_R

    df_R = df_to_R(df)
    kde = ro.r['density'](df_R.rx2(col))
    kde = pd.DataFrame({col:kde.rx2('x'), 'density':kde.rx2('y')})
    kde = kde[(kde[col] >= min(df_R.rx2(col))) &
      (kde[col] <= max(df_R.rx2(col)))]
    return kde

def bin_data_continuous(series, bin_size, max_value):
    col = series.name
    bins = np.arange(0, max_value+bin_size, bin_size)
    binned_data = pd.cut(series, bins, include_lowest=True, right=True)
    binned_data = binned_data.value_counts().sort_index().reset_index()
    binned_data[col] = pd.IntervalIndex(binned_data[col]).mid.values
    return binned_data

def bin_data_discrete(series):
    binned_data = series.value_counts().sort_index().reset_index()
    return binned_data

def append_labels(df, generation, args):
    if args['theta'] == 1:
        model = "random"
    else:
        model = args['mate_choice']
    labels = pd.DataFrame({'generation': generation,
        'model': model, 'theta': args['theta'],
        'seed': args['seed'], 'migration': args['migration_rate'],
        'initial': args['initial_proportion']},
        index=df.index)
    df = pd.concat([df, labels], axis=1)
    return df
