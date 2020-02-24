#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Dependecies
# 1) tophat2
# 2) bamutils
# 3) cutadapt
# 4) prinseq-lite.pl
# 5) bowtie2
# 6) samtools
# 7) bowtie2 index (BOWTIE2_INDEX_BASE2 below): see 'create_human_db_for_bowtie2.sh'

# assuming single-end (SE125 in our case)

# An example input:
#~/HHV6/TopHat2/D701_D502/unmapped/unmapped.bam

# Example outputs:
#~/HHV6/TopHat2/D701_D502/unmapped/unmapped.prep.bt.fq.gz
#~/HHV6/TopHat2/D701_D502/unmapped/unmapped.prep.bt.local.fq.gz

# settings
NT=6                                                      # no. of threads
PROJECT_DIR="~/HHV6"                                      # change as necessary
BOWTIE2_INDEX_DIR=$PROJECT_DIR/Bowtie2Index               # bowtie2 index dir
BOWTIE2_INDEX_BASE2=$BOWTIE2_INDEX_DIR/genome             # prefix for index 
BARCODE="D701_D502"                                       # change as necessary (for this sample)
TOPHAT_OUT2=$PROJECT_DIR/TopHat2/$BARCODE/unmapped        # unmapped top directory
Unmapped_BAM=$TOPHAT_OUT2/unmapped.bam                    # input of preprocessing (unmapped BAM after 2-pass)
Unmapped_FastQ=$TOPHAT_OUT2/unmapped.fq.gz                # fastq format of above (converted from input)
Prep_Unmapped_FastQ=$TOPHAT_OUT2/unmapped.prep.fq.gz      # output of preprocessing
Unmapped_BT2_FastQ=$TOPHAT_OUT2/unmapped.prep.bt2.fq.gz 
Unmapped_BT2_local_FastQ=$TOPHAT_OUT2/unmapped.prep.bt2.local.fq.gz 

# unmapped bam file to fastq
if [ ! -s $Unmapped_FastQ ];then 
    if [ -s $Unmapped_BAM ]; then
        echo -e "\nbamutils tofastq $Unmapped_BAM | gzip -9 > $Unmapped_FastQ"
        time bamutils tofastq $Unmapped_BAM | gzip -9 > $Unmapped_FastQ
    else
        echo -e "\e[031m$Unmapped_BAM not found\e[0m\n"
        exit
    fi

fi

if [ ! -s $Unmapped_FastQ ]; then
    echo -e "\e[031m$Unmapped_FastQ not found\e[0m\n"
    exit
fi

##################################
## 1. Preprocess unmapped.fq.gz ##
## firstly by cutadapt (qc-trim)
## followed by prinseq (lc-trim)
## this setting is based on SE125
##################################
if [ ! -s $Prep_Unmapped_FastQ ];then 
    echo -e "\npreprocessing $Unmapped_FastQ"
    time cutadapt \
        -q 30 \
        -m 50 \
        --max-n 5 \
        --trim-n \
        $Unmapped_FastQ \
        2> ${Unmapped_FastQ%.fq.gz}.cutadapt.log \
        | perl ~/Install/prinseq-lite/prinseq-lite.pl \
        -fastq stdin \
        -min_len 50 \
        -ns_max_n 5 \
        -lc_method dust \
        -lc_threshold 7 \
        -noniupac \
        -trim_tail_left 5 \
        -trim_tail_right 5 \
        -out_good stdout \
        -out_bad null \
        2> ${Unmapped_FastQ%.fq.gz}.cutadapt.prinseq.log \
        | gzip -9 > $Prep_Unmapped_FastQ 
fi

##############################################################
## 2. Map the preprocessed unmapped reads against human ref ##
##############################################################
# via bowtie2 (end-to-end)
    if [ -s $Prep_Unmapped_FastQ ];then 
        echo -e "bowtie2 -U $Prep_Unmapped_FastQ | samtools view -f4 | awk | gzip > $Unmapped_BT2_FastQ"
        time bowtie2 \
            --threads $NT \
            -U $Prep_Unmapped_FastQ \
            -x $BOWTIE2_INDEX_BASE2 \
            2> $TOPHAT_OUT2/unmapped.prep.bt2.log \
            | samtools view -f4 - | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 > $Unmapped_BT2_FastQ
    else
        echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
        exit
    fi

# via bowtie2 --local
    if [ -s $Prep_Unmapped_FastQ ];then 
        echo -e "bowtie2 --local -U $Prep_Unmapped_FastQ | samtools view -f4 | awk | gzip > $Unmapped_BT2_local_FastQ"
        time bowtie2 \
            --local \
            --threads $NT \
            -U $Prep_Unmapped_FastQ \
            -x $BOWTIE2_INDEX_BASE2 \
            2> $TOPHAT_OUT2/unmapped.prep.bt2.local.log \
            | samtools view -f4 - | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 > $Unmapped_BT2_local_FastQ
    else
        echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
        exit
    fi
