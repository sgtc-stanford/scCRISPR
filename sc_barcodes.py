#!/usr/bin/env python
import sys, os, re, pysam, csv, gzip, string, distance
import numpy as np, pandas as pd

DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA") 
A1_3PRIME = 'CTGGCTCCTCTGTATGTTGGAGAAT'
SI_5PRIME = 'AAGCAGTGGTATCAACGCAGAGTAC'
FWD_10X_ADAPTER = 'TTTCTTATATGGG'
REV_10X_ADAPTER = 'CCCATATAAGAAA' 
ILLUMINA_R1_FWD = 'CTACACGACGCTCTTCCGATCT'
ILLUMINA_R1_REV = 'AGATCGGAAGAGCGTCGTGTAG'

ADAPTER_3PRIME = [A1_3PRIME, A1_3PRIME[::-1], A1_3PRIME[::-1].translate(DNA_COMPLEMENT)] 
ADAPTER_5PRIME = [SI_5PRIME, SI_5PRIME[::-1], SI_5PRIME[::-1].translate(DNA_COMPLEMENT)] 

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to be imported/called from long read barcode-matching scripts     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 

def extract_softclips(samrd,offset=0):
  soft_clip_fwd = 'NNNNNNNNN'; soft_clip_rev = 'NNNNNNNNN';
  sam_cigar = samrd.cigartuples
  if sam_cigar and sam_cigar[0][0] == 4:
    soft_clip_len = sam_cigar[0][1]
    soft_clip_fwd = samrd.query_sequence[:soft_clip_len+offset] 
  if sam_cigar and sam_cigar[-1][0] == 4:
    soft_clip_len = sam_cigar[-1][1]
    soft_clip_rev = samrd.query_sequence[(0-soft_clip_len-offset):]
  return {'fwd': soft_clip_fwd, 'rev': soft_clip_rev}

def reverse_complement(seq):
  return seq.translate(DNA_COMPLEMENT)[::-1] 
  
def calc_edit_distance(adapter, softclip_seq, offset):
  dist = []
  for i in range(0,len(softclip_seq)-offset):
    dist.append(distance.levenshtein(adapter, softclip_seq[i:i+offset]))
    if dist[-1] == 0: break
  edit_dist = min(dist)
  bc_pos = dist.index(edit_dist)
  return [edit_dist, bc_pos]
