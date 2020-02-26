#!/usr/bin/env perl
use strict;
use warnings;
use utf8;

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
