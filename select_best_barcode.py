#!/usr/bin/env python
"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgrimes@stanford.edu
:Creation date: 04/19/2021
:Description: 

This script selects the best barcode per read, using edit distance as primary criteria,
  with ties broken using cosine-similarity score.
    
Revisions: 

- mm/dd/yyyy	Description

"""
import argparse, sys, os, re, csv, gzip, string
import numpy as np, pandas as pd

import sc_barcodes as scb

script_name = os.path.basename(__file__)
print("Running ", script_name)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read program arguments and file(s)                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 
def parse_commandline():
  parser=argparse.ArgumentParser()
  parser.add_argument('--input', '-i', help='barcodes file (from soft-clips)', type=str, required=True)
 
  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args

args = parse_commandline()
out_fn = args.input.split('.')[0] + '.barcode_match.tsv'
df_softclip = pd.read_csv(args.input, sep='\t', header=0)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Functions for determining best barcode match                                #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
def find_best_barcode(sc_rows):
  best_rows = sc_rows.loc[sc_rows['dist'] == min(sc_rows['dist'])]
  if len(best_rows) > 1:
    best_cosine_sim = best_rows.loc[best_rows['score'] == max(best_rows['score'])]
    best_row = best_cosine_sim.iloc[0,]
  else:
    best_row = best_rows.iloc[0,]
  return best_row
    
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Main program logic                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
df_best_barcodes = df_softclip.groupby('rd_name').apply(find_best_barcode)
df_best_barcodes.loc[df_best_barcodes['dist'] < 5].to_csv(out_fn, sep='\t', index=False)
