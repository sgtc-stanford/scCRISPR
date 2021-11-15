#!/usr/bin/env python
"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgrimes@stanford.edu
:Creation date: 03/24/2021
:Description: 

This script extracts soft clipped bases at beginning (FWD strand) or end (REV strand)
of read.  These sequences will subsequently be searched for expected single cell barcode sequences.
    
Revisions: 

- 03/26/2021	Import reusable methods from sc_barcodes
- 04/28/2021	Add nbest command line argument

"""
import argparse, sys, os, re, pysam, csv, gzip, string, distance
import numpy as np, pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse import csr_matrix

import sc_barcodes as scb

script_name = os.path.basename(__file__)
print("Running ", script_name)

#Use SEARCH_OFFSET is desired to ignore part of 13bp 10X adapter, and 10bp UMI during matching processs 
SEARCH_OFFSET = 0
MAX_SEARCH_LEN = 55-SEARCH_OFFSET  
MIN_SCORE = 0.0

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define internal modules                                                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 
def parse_commandline():
  parser=argparse.ArgumentParser()
  parser.add_argument('--bam', '-b', help='input bam file', type=str, required=True)
  parser.add_argument('--barcodes', '-c', help='cellranger barcodes file', type=str, required=True)
  parser.add_argument('--strand', '-s', help='gene strand', type=str, required=True, choices=['plus','minus','both'])
  parser.add_argument('--exonrds', '-x', help='reads with exon skipping pattern identified', type=str, required=False)
  parser.add_argument('--kmer_len', '-k', help='k-mer length', type=int, required=False, default=8)
  parser.add_argument('--nbest', '-n', help='number of best matches to evaluate', type=int, required=False, default=5)

  args=parser.parse_args()
  print(args, file=sys.stderr)
  return args
  
def debug_samrd(samrd):
  strand = 'Rev' if samrd.is_reverse else 'Fwd'
  return [samrd.qname[-16:], strand, samrd.tid, samrd.pos]
  
def best_barcodes(string, barcodes, barcode_tfidf, vectorizer, nbest=5):
  best_barcodes = [['N',0]]
  barcode_seq_tfidf = vectorizer.transform([string])
  cos_sim = cosine_similarity(barcode_seq_tfidf, barcode_tfidf, dense_output=False)
  
  non_zero = [((i, j), cos_sim[i,j]) for i, j in zip(*cos_sim.nonzero())]
  nz_sorted = sorted(non_zero, key=lambda x: -x[1])
  idx_nbest = [x[0][1] for x in nz_sorted[0:nbest] if x[1] > MIN_SCORE]
  
  if len(idx_nbest) > 0:
    best_barcodes = zip([barcodes[i] for i in idx_nbest], [cos_sim[(0,i)] for i in idx_nbest])
  
  return best_barcodes

def format_bc_string(soft_clips, strand, bc_start):
  bc_end = pos+16
  if strand == 'fwd':  #positions will all be -ve offsets from end of sequence
    r1_adapter = soft_clips[max(bc_start-22, 0):bc_start]
    barcode = soft_clips[bc_start:bc_end]
    umi = soft_clips[bc_end:min(bc_end+10, len(soft_clips))]
    return '|'.join([r1_adapter, barcode, umi])
  else:
    umi = soft_clips[max(bc_start-10, 0):bc_start]
    barcode = soft_clips[bc_start:bc_end]
    r1_adapter = soft_clips[bc_end:min(bc_end+22, len(soft_clips))]
    return '|'.join([scb.reverse_complement(r1_adapter), scb.reverse_complement(barcode), scb.reverse_complement(umi)])

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
args = parse_commandline()

KMER_LEN = args.kmer_len
NBEST = args.nbest

sam_input = pysam.Samfile(args.bam,'rb') if args.bam[-3:] == 'bam' else pysam.Samfile(args.bam,'r')
sam_fname = os.path.basename(args.bam)

if args.barcodes.endswith('gz'):
  barcode_input = gzip.open(args.barcodes, 'r')
  barcodes = [line[:16].decode('ascii') for line in barcode_input]
else:
  barcode_input = open(args.barcodes, 'r')
  barcodes = [line[:16] for line in barcode_input]

xskip_fn = args.exonrds

try:
  out_fn = sam_fname[0:-4] + '.softclip.bestN.txt'
  out_file = open(out_fn, 'w')
  out_csv = csv.writer(out_file, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
except:
  print("Unable to open text file for output: ", out_fn)
  sys.exit(1)
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Vectorize barcodes                                                          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
def ngrams(string, kmer_len=KMER_LEN):
  ngrams = zip(*[string[i:] for i in range(kmer_len)])
  return [''.join(ngram) for ngram in ngrams]
  
#build kmer dictionary of all barcode seqs (in both forward and reverse orientation)
vectorizer_fwd = CountVectorizer(min_df=1, analyzer=ngrams)
fwd_barcode_tfidf = vectorizer_fwd.fit_transform(barcodes)

vectorizer_rev = CountVectorizer(min_df=1, analyzer=ngrams)
rev_barcode_tfidf = vectorizer_rev.fit_transform([scb.reverse_complement(barcode) for barcode in barcodes])

if args.exonrds: 
  xskip_rdnames = pd.read_csv(xskip_fn, sep='\t', header=None, index_col=0)
   
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read sam file and check for soft-clips at beginning or end of read          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Expected sequence (FWD), soft clips indicated by []
#   [Illumina R1 adapter][cell barcode][UMI][10X internal adapter] cDNA [polyA][adapter]
#
# Expected sequence (REV), soft clips indicated by []
#   [adapter][polyT] cDNA [10X internal adapter][UMI][barcode][Illumina R1 adapter]
i = 0; tot_rds = 0; adapter_flank_bp = 6; last_qname = 'none';

out_csv.writerow(['rd_name', 'exon_skip', 'strand', 'barcode', 'score', 'dist', 'pos', 'adapter|BC|UMI', 'search_len', 'align_start', 'align_end'])
sc_3or5 = '5prime' if args.strand == 'plus' else '3prime'

for samrd in sam_input.fetch(until_eof=True):
  i += 1
  
  if samrd.is_secondary:
    continue

  if args.exonrds and samrd.qname in xskip_rdnames.index:
    xskip_pattern = xskip_rdnames.loc[samrd.qname,1]
  elif args.exonrds:
    continue
  else:
    xskip_pattern = None

  align_strand = 'minus' if samrd.is_reverse else 'plus'

  tot_rds += 1
  soft_clips = scb.extract_softclips(samrd)
 
  barcodes_scores = [['N', 0]]

  if args.strand == 'plus':
    sc_5prime_len = len(soft_clips['fwd'])
    if sc_5prime_len > 16+SEARCH_OFFSET+1:    
      #Working backwards from end of soft clipped sequence (s/b 10X adapter, then UMI, then cell barcode) - determine start position for cell barcode search
      i_start = max(sc_5prime_len-MAX_SEARCH_LEN-SEARCH_OFFSET, 0)
      search_seq = soft_clips['fwd'][i_start:-SEARCH_OFFSET] if SEARCH_OFFSET > 0 else soft_clips['fwd'][i_start:]
      barcode_scores = best_barcodes(search_seq, barcodes, fwd_barcode_tfidf, vectorizer_fwd, NBEST)
    
  else: # args.strand == 'minus': 
    sc_3prime_len = len(soft_clips['rev'])
    if sc_3prime_len > 16+SEARCH_OFFSET+1:     
      i_end = min(MAX_SEARCH_LEN+SEARCH_OFFSET, sc_3prime_len)  
      search_seq = soft_clips['rev'][SEARCH_OFFSET:i_end]
      barcode_scores = best_barcodes(search_seq, barcodes, rev_barcode_tfidf, vectorizer_rev, NBEST)
  
  for bc_score in barcode_scores:
    if bc_score[0] != "N":
      if args.strand == 'plus':
        [dist, pos] = scb.calc_edit_distance(bc_score[0], search_seq, 16)
        barcode_with_flanking = format_bc_string(search_seq, 'fwd', pos)
        bc_pos = pos-len(search_seq)
      else:  #args.strand == 'minus'
        [dist, pos] = scb.calc_edit_distance(scb.reverse_complement(bc_score[0]), search_seq, 16)
        barcode_with_flanking = format_bc_string(search_seq, 'rev', pos)
        bc_pos = pos  
      out_csv.writerow([samrd.qname, xskip_pattern, align_strand, bc_score[0], bc_score[1], dist, bc_pos, barcode_with_flanking, len(search_seq),
	                    samrd.reference_start, samrd.reference_end])
      if dist < 2: break

print(i, "sam records read")
print("Evaluated", tot_rds, "primary (full transcript) alignments")

for fh in [sam_input, barcode_input, out_file]:
  fh.close()
