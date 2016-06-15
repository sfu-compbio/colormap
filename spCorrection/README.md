Shortest path correction
========================
For each long read, fixes the errors in each connected component of mapped short reads using the shortest path in that connected component.

USAGE: ./spCorrection -l LR.fasta -a SR2LR.sam > LR_corr.fasta
SR2LR.sam: 	alignment of SRs to LRs in sam format obtained from BWA-MEM; should be sorted using samtools sort 
 LR.fasta: 	long reads in FASTA format

to obtain mappings file:
bwa index LR.fasta
bwa mem -aY -t <threads> -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 16 -w 40 -r 1 -D 0 -y 20 -L 30,30 -T 2.5 LR.fasta SR.fastq > SR2LR.sam
blasr illumina.fastq pacbio.fasta -noSplitSubreads -m 5 -out sr2lr.m5 -nproc 8 -bestn 80 -nCandidates 80
to sort based on the ids of long reads:
samtools view -bS SR2LR.sam | samtools sort -@ <trheads> -o - SR2LR.sort.sam | samtools view -h -o SR2LR.sort.sam -
