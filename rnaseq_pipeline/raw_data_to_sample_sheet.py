#!/usr/bin/python3
# -*- coding: utf-8 -*

"""
    SAMPLESHEET BUILD
"""

import argparse
import glob
import pandas as pd
import yaml
import sys
from pathlib import Path

PARSER = argparse.ArgumentParser(description='Build sample_sheet parse')

PARSER.add_argument('-p', '--path', help='Relative Path to raw_data')
PARSER.add_argument('-r', '--ref', help='genome reference')
PARSER.add_argument('-k', '--kallisto', default="", help='kallisto params')

ARG = PARSER.parse_args()

class LenWrong(Exception):
    pass

def load_conf(in_file):
    cwd = Path.cwd()
    with cwd.joinpath(in_file).open(mode = 'r') as conf:
        conf_load = yaml.safe_load(conf)

    return conf_load

def write_conf(data, out_file, mode):
    cwd = Path.cwd()
    with cwd.joinpath(out_file).open(mode = mode) as conf:
        yaml.dump(data, conf, default_flow_style=False, allow_unicode=True, indent = 4)

def set_conf(sample_name, fq, conf, file):
    if sample_name+fq not in conf['rq']:
        conf['rq'].update({ sample_name+fq:[]})
    if file not in conf['rq'][sample_name+fq]:
        conf['rq'][sample_name+fq].append(file)

def build(path, ref, kallisto):

    build = {'sample_name':[], 'fq':[]}
    black = {}

    cwd = Path.cwd()
    path = cwd.joinpath(path)
    files = [file for file in path.iterdir() if '.fastq' in file.name]

    conf = load_conf('config/template/config.yaml')

    conf['resources']['ref']['transcriptome'] = ref
    conf['kallisto']["params"] = kallisto
    conf.update({'rq':{}})

    conf['path'] = str(Path.cwd())

    for file in files:

        filename = file.name.split('.')[0]
        sample_name = filename.split('_')[0]
        pair = filename.split('_')[1]

        print(f"Path : {sample_name} , pair : {pair}")
        if sample_name not in black:
            black.update({sample_name:{'R1':0, 'R2':0}})

        if pair == 'R1':
            black[sample_name]['R1'] += 1
            set_conf(sample_name, '_R1', conf, str(file))

        if pair == 'R2':
            black[sample_name]['R2'] += 1
            set_conf(sample_name, '_R2', conf, str(file))

    try:
        for u in black:
            if black[u]['R1'] != black[u]['R2'] and black[u]['R2'] != 0:
                raise LenWrong

    except LenWrong:
        print('WARNING : Multiple sample from different line had not the same method (PE or SE)')
        sys.exit(1)

    for u in black:
        build['sample_name'].append(u)
        build['fq'].append(None) if black[u]['R2'] == 0 else build['fq'].append('PE')

    write_conf(conf, 'config/config.yaml', 'w+')

    out_data_frame = pd.DataFrame.from_dict(build)
    out_data_frame.to_csv('config/units.tsv', sep='\t', index=False)

if __name__ == "__main__":
    build(ARG.path, ARG.ref, ARG.kallisto)
