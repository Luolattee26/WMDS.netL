# This code is used to convert gtf files into transcript fasta files using cufflinks on a Linux platform, 
# so that the fasta information of all involved lncRNAs can be extracted for protein coding capacity analysis 
# in the next step of the operation



# gtf to fasta by gffread(cufflinks)
# use original hg38(UCSC)
gffread gencode.v22.annotation.gtf -g hg38.fa -w transcripts.fa

# use latest hg38(UCSC)
gffread gencode.v22.annotation.gtf -g hg38_latest.fa -w transcripts_latest.fa