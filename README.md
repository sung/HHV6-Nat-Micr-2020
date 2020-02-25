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
The remaining reads of each sample were mapped to a custom Kraken reference database, including the default bacterial and viral genomes and few additional eukaryotic genomes to remove residual unmapped human reads. Kraken (v0.10.6) was run using the metagm_run_kraken option.

### Step 1: Run Kraken 
Using the `-noclean` option on the paired filtered reads:
`metagm_run_kraken -t 4 -noclean scanner_DB Sample(x).report Sample(x)_kneaddata_paired_1.fastq Sample(x)_kneaddata_paired_2.fastq`

This will generate both a .report file and a .report.kraken_out file.

### Step 2: Open the Sample(x).report file. 
Write down all of the ID’s of the group of interest that you want to extract (penultimate column). This list can be very short or somewhat long, like in the case of E. coli. Make a column in a txt file of this list (Target_ids.tab).

### Step 3: Run the Perl script (shown below) on all of the Sample(x).report.kraken_out files.
```perl 
#!/usr/bin/perl
#use strict;

my $path = $ARGV[0];
if (scalar($path)==0) {
    $path = `pwd`;
    chomp $path;
}

undef @files;
opendir PATH, $path;
@files = grep /\.kraken_out/, readdir PATH;
closedir PATH;

my $ids = "";
open FILE, "$path\/Target\_ids\.tab";
while (<FILE>) {
    chomp;
    $ids .= "\|$_\|";
}
close FILE;

foreach (@files){
    my $f = $_;
    open OUT, ">$path\/$f\.selection";
    open FILE, "$path\/$f";
    while (<FILE>) {
        my $line = $_;
        my @line = split /\t/, $line;
        if ($ids =~ /\|$line[2]\|/) {
            print OUT "$line[1]\n";
        }
    }
    close FILE;
    close OUT;
}
```
 
### Step 4: Unix: Extract read pair headers of all reads of each sample: 
`less Sample(x)_kneaddata_paired_1.fastq | grep @HX1 > Sample(x).tab`

### Step 5: Run the following R code using the output from steps 3 and 4:
```r 
read_names <- scan("Sample(x).tab", what="\n")
selection <- scan("Sample(x).report.kraken_out.selection", what="\n")
selection <- gsub(selection, pattern="/1", replacement="")
selection <- as.numeric(selection)

read_selection <- read_names[selection]
read_selection <- gsub(read_selection, pattern="@", replacement="")

write(read_selection, file="Sample(x)Target.txt", ncolumns=1, sep="\t")
```

### Step 6: Unix: Extract all the target reads (forward and reverse).
```bash 
seqtk subseq Sample(x)_kneaddata_paired_1.fastq Sample(x)Target.txt > Sample(x)Target1.fq
seqtk subseq Sample(x)_kneaddata_paired_2.fastq Sample(x)Target.txt > Sample(x)Target2.fq
```
 
## A list of dependent software and database
- [tophat2](https://github.com/infphilo/tophat)
- [bamutils](https://genome.sph.umich.edu/wiki/BamUtil)
- [cutadapt](https://github.com/marcelm/cutadapt)
- [prinseq-lite.pl](http://prinseq.sourceforge.net/)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [Kraken](https://ccb.jhu.edu/software/kraken/)


## Contacts
+ [Sung Gong](https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/)
+ [Marcus de Goffau](https://www.sanger.ac.uk/people/directory/de-goffau-marcus)
