#!/usr/bin/env python
"""

:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgrimes@stanford.edu
:Creation date: 02/10/2020
:Description: 

This script counts exons for a given gene, in each long-read transcript.
	
Revisions:
	5/19/2020: Modify for situation where exons to compare do not start at 1
	8/24/2020: Optionally require only that reads align to last exon (vs first thru last)

"""
import sys, os, re, pysam, csv, numpy as np

script_name = os.path.basename(__file__)
print("Running ", script_name)

MIN_OVERLAP = 12  
CHKRD_NAMES = ['0c3ff184-b837-4671-a55a-24c7e5bc3df3', '57891cde-7423-48b7-a341-52ace43eba27']

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Files are assumed to be in current path (or fully qualified path is given)
if len(sys.argv) < 3:
  print("Usage: ", script_name, "<sam_file|bam_file> <gene_exon.bed> <all|last> [dtl]")
  sys.exit(1)

sam_file = sys.argv[1]
if os.path.isfile(sam_file) and os.access(sam_file, os.R_OK):
  sam_input = pysam.Samfile(sam_file,'rb') if sam_file[-3:] == 'bam' else pysam.Samfile(sam_file,'r')
  sam_fname = os.path.basename(sam_file)
else:
  print("Unable to open sam/bam file for input:", sam_file)
  sys.exit(1)
  
gene_exon_bed = sys.argv[2]
if os.path.isfile(gene_exon_bed) and os.access(gene_exon_bed, os.R_OK):
  gene_input = open(gene_exon_bed, 'r')
else:
  print("Unable to open gene_exon bed file for input:", gene_exon_bed)
  sys.exit(1)

try:
  out_file = sam_fname[0:-4] + '.exon_cts.txt'
  transcript_out = open(out_file, 'w')
except:
  print("Unable to transcript/exon counts file for output: ", out_file)
  sys.exit(1)

# Fix this later - create actual named parameters 
exon_span = 'all'; write_dtl = False;
if len(sys.argv) > 3:
  if sys.argv[3] == 'dtl':
    write_dtl = True
  else:
    exon_span = sys.argv[3]
    if len(sys.argv) > 4: write_dtl = True 

if write_dtl:
  out_detail = sam_fname[0:-4] + '.exon_reads.txt'
  reads_out = open(out_detail, 'w')
  reads_csv = csv.writer(reads_out, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define internal modules                                                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 
def createExonList(gene_input):
  exon_rds = list(csv.reader(gene_input, delimiter="\t"))
  for exon in exon_rds:
    exon[1] = int(exon[1])
    exon[2] = int(exon[2])
    exon[3] = int(exon[3][4:])
  return exon_rds 

def getTranscriptSpan(gene_exons, min_overlap, exon_span):
  gene_strand = 1 if gene_exons[-1][3] > gene_exons[0][3] else -1
  if exon_span == 'all':
    max_start = gene_exons[0][2] - min_overlap
    min_end   = gene_exons[-1][1] + min_overlap
  else:
    last_exon = gene_exons[-1] if gene_strand == 1 else gene_exons[0]
    max_start = last_exon[2] - min_overlap
    min_end   = last_exon[1] + min_overlap
  return [gene_strand, max_start, min_end]
  
def isSpanningRead(exon_rds, gene_strand, exon_span):
  if exon_span == 'all':
    spanning_read = (exon_rds[0] > 0 and exon_rds[-1] > 0)
  elif gene_strand == 1:
    spanning_read = exon_rds[-1] > 0
  else:
    spanning_read = exon_rds[0] > 1
  return spanning_read
 
# function to find first index >= x 
def lowerIndex(arr, n, x): 
  l = 0
  h = n-1
  while (l <= h): 
    mid = int((l + h)/2) 
    if (arr[mid] >= x): 
      h = mid - 1
    else: 
      l = mid + 1
  return l 
  
  
# function to find last index <= x 
def upperIndex(arr, n, x): 
  l = 0
  h = n-1
  while (l <= h): 
    mid = int((l + h)/2) 
    if (arr[mid] <= x): 
      l = mid + 1
    else: 
      h = mid - 1
  return h 
  
  
# function to count elements within given range 
def countInRange(arr, n, x, y): 
  # initialize result 
  count = 0; 
  count = upperIndex(arr, n, y) - lowerIndex(arr, n, x) + 1; 
  return count

def printDebug(i, qname, chkrdnames):
  return True if (i <20 or qname in chkrdnames) else False; 
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read sam file and check if spans complete genomic coordinate range (from .bed)
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
i = 0; tot_rds = 0;
transcript_exon_cts = []

# RACK1 example exon list:
# [['5', 181236937, 181237042, 8], ['5', 181237609, 181237719, 7], ['5', 181238099, 181238239, 6], 
#  ['5', 181239067, 181239177, 5], ['5', 181239487, 181239582, 4], ['5', 181241492, 181241639, 3], 
#  ['5', 181242174, 181242345, 2], ['5', 181243692, 181244209, 1]]
gene_exons = createExonList(gene_input)
print(gene_exons)

gene_chr = gene_exons[0][0]
transcript_info = getTranscriptSpan(gene_exons, MIN_OVERLAP, exon_span)
gene_strand = transcript_info[0]
transcript_max_start = transcript_info[1]
transcript_min_end   = transcript_info[2]

for samrd in sam_input.fetch(until_eof=True): 
  tot_rds += 1
  if samrd.is_unmapped:
    continue
  else:
    chrnm = sam_input.getrname(samrd.tid)
    if chrnm == gene_chr:
      #print("Gene transcript starting at: {0}, ending at: {1}".format(samrd.reference_start, samrd.reference_end))
      if samrd.qname in CHKRD_NAMES or (samrd.reference_start <= transcript_max_start and samrd.reference_end > transcript_min_end):
        i += 1
        if printDebug(i, samrd.qname, CHKRD_NAMES):
          print("{2} Full transcript start: {0}, end: {1}.  Alignment per exon below".format(samrd.reference_start, samrd.reference_end, samrd.qname))

        transcript_rds = samrd.get_reference_positions()
        exon_rds = np.zeros((len(gene_exons),), dtype=int)
        for j, exon in enumerate(gene_exons):
          nr_exon_bases = countInRange(transcript_rds, len(transcript_rds), exon[1], exon[2])
          exon_rds[j] = nr_exon_bases
          if printDebug(i, samrd.qname, CHKRD_NAMES):
            print("Exon {0}: {1} aligned bases".format(exon[3], nr_exon_bases))

        if printDebug(i, samrd.qname, CHKRD_NAMES):
          print(exon_rds)
        
        if isSpanningRead(exon_rds, gene_strand, exon_span):
          exon_rds_bin = np.where(exon_rds > 0, 1, 0)
          transcript_exon_cts.append(exon_rds_bin)
          if write_dtl:
            exon_skip_str = [str(exon[3]) for exon in gene_exons]
            for k, num in enumerate(exon_rds_bin):
              if num == 0: exon_skip_str[k] = '-'
            reads_csv.writerow([samrd.qname, ''.join(exon_skip_str)])
            
print("{0} total reads processed; {1} gene transcripts".format(tot_rds, i))
total_exon_cts = sum(transcript_exon_cts)
print(total_exon_cts.tolist())

if len(gene_exons[0]) > 4:
  exon_names = [exon[4] + '(' + str(exon[3]) + ')' for exon in gene_exons]
else:
  exon_names = ['Exon' + str(exon[3]) for exon in gene_exons]

out_csv = csv.writer(transcript_out, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
out_csv.writerow(exon_names)
out_csv.writerow(total_exon_cts.tolist())

sam_input.close()
gene_input.close()
transcript_out.close()
