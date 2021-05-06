1) [optional] Determine transcript isoform per long read, for a specific gene
   Requires .bed file with start/end positions for each exon of interest
   
   python ${CODE_DIR}/transcript_exon_cts.py <sample>.bam <gene>_exons.bed dtl

2) Determine up to top5 barcode matches per long read, using vectorization and cosine similarity metric
   Required input files are bam file, and barcode 'whitelist' from short-reads.
   If transcript isoform file is specified, also merge in the transcript isoform per read, and exclude
     any reads which do not span the exons of interest.
   Parameters available with --help, or see below.

   python softclip_bestN_barcodes.py -b <sample>.bam -c <sample>.barcodes.tsv.gz -s plus

   '--bam', '-b', help='input bam file', type=str, required=True
   '--barcodes', '-c', help='cellranger barcodes file', type=str, required=True
   '--strand', '-s', help='gene strand', type=str, required=True, choices=['plus','minus','both']
   '--exonrds', '-x', help='reads with exon skipping pattern identified', type=str, required=False
   '--kmer_len', '-k', help='k-mer length', type=int, required=False, default=8
   '--nbest', '-n', help='number of best matches to evaluate', type=int, required=False, default=5

3) Pick best of top5 barcode matches, using edit distance as primary criterion (discard reads with edit distance > 4)
   python select_best_barcode.py -i <sample>.softclip.bestN.txt 

4) Final filtering of barcode matches.  Keep all barcode matches with edit distance < 3, and keep edit distance 3,4
     only if cosine similarity metric > 0.15.  Reorder columns and fix column heading for later processing.

   awk -F'\t' '{if(NR == 1 || ($6 < 3 || $5 > 0.15)) print $4, $2, $5, $3, $6, $7, $8, $9}' \
     <sample>.barcode_match.tsv > <sample>.barcode_info.txt
   sed -i '1 s/^.*$/best_barcode exon_skip max_score strand edit_dist start_pos barcode+flanking search_len/' <sample>.barcode_info.txt

   