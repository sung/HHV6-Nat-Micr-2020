## Introduction
This is a collection of custom scripts used for Gaccioli et al. Nat. Microbiol. 2020.

## Initial mapping of RNA-Seq data
The RNA-Seq reads were trimmed (using [cutadapt](https://github.com/marcelm/cutadapt) and [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) and mapped to the primary chromosomal assemblies of the [GRCh38.p3 version](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.29/) of the human reference genome using [TopHat2](https://github.com/infphilo/tophat), a splice-aware mapper built on top of [Bowtie2](https://github.com/BenLangmead/bowtie2) short-read aligner. For more details, read [Gong et al., Epigenetics, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5989156/) and [Gong et al., JCI Insight, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124516/).

## Preprocessing of unmapped RNA-Seq data 
The initial ‘unmapped’ reads were filtered out to remove poor quality reads, based on the following conditions: 
1. base-quality score (i.e. [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score)) < 30, 
2. read-length < 50bp, 
3. undetermined base (i.e. Ns) > 5bp, and 
4. poly A/T > 5bp, and 
5. low-complexity reads (defined by [the dust score](https://www.ncbi.nlm.nih.gov/pubmed/16796549)) > 7. 

Following shell script handles what's been described above: `preprocess_unmapped_reads.sh`

## A custom human sequence DB
In order to remove as many reads of human origin as possible, additional human reads were subtracted if they aligned to sequences present in the following databases: 
1. [GRCh38.p5](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.31/), 
2. [human RefSeq](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/), 
3. all human contigs and clone sequences from NCBI NT. 

See `create_human_db_for_bowtie2.sh` for details.

## Post processing
The remaining reads of each sample were mapped to a custom Kraken reference database, including the default bacterial and viral genomes and few additional eukaryotic genomes to remove residual unmapped human reads. Kraken (v0.10.6) was run using the metagm_run_kraken option. **Further details to follow**

## A list of dependent software and database
- [tophat2](https://github.com/infphilo/tophat)
- [bamutils](https://genome.sph.umich.edu/wiki/BamUtil)
- [cutadapt](https://github.com/marcelm/cutadapt)
- [prinseq-lite.pl](http://prinseq.sourceforge.net/)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [Kraken](https://ccb.jhu.edu/software/kraken/)


## Contacts
[Sung Gong](https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/)
[Marcus de Goffau](https://www.sanger.ac.uk/people/directory/de-goffau-marcus)
