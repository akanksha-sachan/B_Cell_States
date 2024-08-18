#!/usr/bin/env python3
import argparse
import sys
from os.path import join as pjoin, exists as pexists
from os import listdir
from functools import reduce
from operator import add
from collections import Counter
import logging
import itertools

parser = argparse.ArgumentParser(description="Validates cell names in context-specific GRNs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir_data', type=str, default='data', help='Directory of input data folder.')
parser.add_argument('--dir_makefiles', type=str, default='makefiles', help='Directory of makefiles folder.')
parser.add_argument('-c', action='store_const', default=False, const=True, help='Continue validation as much as possible except for critical errors.')

args = parser.parse_args()
dirdata = args.dir_data
dirmakefiles = args.dir_makefiles
nerr = 0 if args.c else None

# Load cell names and validate context-specific GRNs
if not pexists(pjoin(dirmakefiles, 'static.mk')) or not pexists(pjoin(dirdata, 'subsets.txt')):
    logging.warning('Cannot find static.mk or subsets.txt. Skipping static network inference checks.')
else:
    # Load subsets
    with open(pjoin(dirdata, 'subsets.txt'), 'r') as f:
        names = f.readlines()
    names = [x.strip() for x in names]
    names = list(filter(lambda x: len(x) > 0, names))

    # Load cell names in each subset
    namec_srna = []
    namec_satac = []
    for xi in names:
        with open(pjoin(dirdata, 'subsets', xi, 'names_rna.txt'), 'r') as f:
            t1 = f.readlines()
        t1 = [x.strip() for x in t1]
        t1 = list(filter(lambda x: len(x) > 0, t1))
        namec_srna.append(t1)
        with open(pjoin(dirdata, 'subsets', xi, 'names_atac.txt'), 'r') as f:
            t1 = f.readlines()
        t1 = [x.strip() for x in t1]
        t1 = list(filter(lambda x: len(x) > 0, t1))
        namec_satac.append(t1)
    
    snamec_srna, snamec_satac = [[Counter(y) for y in x] for x in [namec_srna, namec_satac]] ###### counter object is unhashable

    # Validate RNA cell names
    t1 = reduce(add, snamec_srna)
    t2 = [x[0] for x in dict(t1).items() if x[1] > 1][:3]
    if len(t2) > 0:
        t2 = [(x, list(itertools.chain.from_iterable([[names[y]] * snamec_srna[y][x] for y in range(len(names))]))) for x in t2]
        t2 = '; '.join(['{}: {}'.format(x[0], ','.join(x[1])) for x in t2])
        logging.warning('Found RNA cells appearing multiple times in subsets at subsets/*/names_rna.txt. First three cells and their assigned subsets: ' + t2)

    # Extract keys from the Counter objects before converting to sets
    t1_keys = set(t1.keys())
    snamec_srna_keys = set(itertools.chain.from_iterable(counter.keys() for counter in snamec_srna))
    t2 = list(t1_keys - snamec_srna_keys)[:3]
    if len(t2) > 0:
        s = 'Subset RNA cells in subsets/*/names_rna.txt could not be found in expression.tsv.gz. First three missing cells: ' + ', '.join(t2)
        if nerr is None:
            raise ValueError(s)
        else:
            logging.error(s)
            nerr += 1

    # Validate ATAC cell names
    t1 = reduce(add, snamec_satac)
    t2 = [x[0] for x in dict(t1).items() if x[1] > 1][:3]
    if len(t2) > 0:
        t2 = [(x, list(itertools.chain.from_iterable([[names[y]] * snamec_satac[y][x] for y in range(len(names))]))) for x in t2]
        t2 = '; '.join(['{}: {}'.format(x[0], ','.join(x[1])) for x in t2])
        logging.warning('Found ATAC cells appearing multiple times in subsets at subsets/*/names_atac.txt. First three cells and their assigned subsets: ' + t2)

    # Extract keys (cell names) from the Counter objects before converting to sets
    t1_keys = set(t1.keys())
    snamec_satac_keys = set(itertools.chain.from_iterable(counter.keys() for counter in snamec_satac))

    # Check for missing ATAC cells
    t2 = list(t1_keys - snamec_satac_keys)[:3]
    if len(t2) > 0:
        s = 'Subset ATAC cells in subsets/*/names_atac.txt could not be found in bams folder. First three missing cells: ' + ', '.join(t2)
        if nerr is None:
            raise ValueError(s)
        else:
            logging.error(s)
            nerr += 1

if nerr is not None and nerr > 0:
    raise RuntimeError(f'Found {nerr} error(s) in total.')
