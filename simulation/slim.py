#! /usr/bin/python3

import os
import subprocess
from numpy import isin

def construct_slim_command(args):
    slim_args = ['recombination_rate_file', 'initial_proportion',
                 'migration_rate', 'population_size', 'n_generations',
                 'outdir', 'theta', 'alpha', 'sigma']

    slim_args = {k:v for k,v in args.items()
                     if isin(k, slim_args)}

    cmd = ['slim -s ' + str(args['seed'])]
    cmd +=  ["-d " + k.upper() + "='" + str(v) + "'"
            if isinstance(v, str)
            else '-d ' + k.upper() + '=' + str(v)
            for k,v in slim_args.items()]
    cmd += [args['outdir'] + '/model.slim']
    cmd = ' '.join(cmd).split()
    return cmd

def make_output_dir(name):
    cmd = 'mkdir -p ' + name
    subprocess.call(cmd, shell=True)

def write_slim(parent_sim, mate_choice, output_dir):
    script_path = os.path.dirname(parent_sim)
    mate_choice_script = (script_path + '/mate_choice_' +
                          mate_choice + '.txt')
    if not os.path.isfile(mate_choice_script):
        raise ValueError("%s does not exist."%(mate_choice_script))
    subprocess.run([script_path + '/write_slim_file.sh',
      parent_sim, mate_choice_script, output_dir])

def run_slim(args):
    make_output_dir(args['outdir'])
    write_slim(args['source'], args['mate_choice'], args['outdir'])
    cmd = construct_slim_command(args)
    slim = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    err = slim.communicate()[1]
    return err
