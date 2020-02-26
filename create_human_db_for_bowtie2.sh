#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>

PROJECT_DIR="~/HHV6"                                      # change as necessary
BOWTIE2_INDEX_DIR=$PROJECT_DIR/Bowtie2Index               # bowtie2 index dir
BOWTIE2_INDEX_BASE2=$BOWTIE2_INDEX_DIR/genome             # prefix for index 

mkdir -p $BOWTIE2_INDEX_DIR 
cd $BOWTIE2_INDEX_DIR 

# download relevant files
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_rna.fna.gz

# human 'contig' and 'clone' sequences
time zcat nt.gz | awk 'BEGIN{ ORS = ""; RS = ">"; FS="\n" }($1~"\\| Homo sapiens" || $1~"\\| Human") && ($1~"contig" || $1~"Contig" || $1~"clone" || $1~"Clone"){ print ">" $0 }' > genome.contig.clone.fa

# merge files
zcat GCF_000001405.31_GRCh38.p5_genomic.fna.gz GCF_000001405.31_GRCh38.p5_rna.fna.gz | cat - genome.contig.clone.fa > genome.fa

# build samtools index 
samtools faidx genome.fa
# build bowtie2 index
time bowtie2-build genome.fa $BOWTIE2_INDEX_BASE2
