#!/usr/bin/env perl
use warnings;
use strict;

my ($input, $dir_str) = @ARGV;
die "usage: $0 <list of all sim files> <dir string for sim files>" unless @ARGV == 2;

open(IN, "<", $input) or die "Could not open $input:$!\n";
open(OUT, ">", "$input.meta") or die "could not open output_file:$!\n";

while(<IN>) {
	my $line = $_;
	chomp $line;
	my ($attr_val, $rep_val) = $line =~ /attr_(.+)_rep_(\d+)_/;
	
	my $new_line = "$dir_str/" . $line . "\t$attr_val\t$rep_val\n";
	print OUT $new_line;
}

close IN;
close OUT;
