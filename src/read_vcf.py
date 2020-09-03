'''
FILE: read_vcf.py
PURPOSE: reads vcf files
INPUT: path to file
OUTPUT: opened vcf file
CREATED BY: GitHub user dceoy
NOTES: repository can be found here: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
'''

#!/usr/bin/env python

import io
import os
import pandas as pd

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    