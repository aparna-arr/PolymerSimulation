#!/usr/bin/env perl
use warnings;
use strict;

my ($template, $outdir, $outfile_basename, @percents) = @ARGV;
die "usage: generate_domain_perc_dist_input_files.pl <template domain input file> <output directory> <output file basename> <perc1> <perc2> <...>\n

Note: the array of percentages is intended for ONE DOMAIN or multiple domains that all increase/decrease % in the SAME WAY at the SAME TIME. No matrix will be produced if you have multiple domains!! Read the code if you don't know what this means because it will mess up your output!\n" unless @ARGV >= 4;

open(TMPL, "<", $template) or die "could not open file $template: $!\n";
my @domainfiles;
my @filenames;

for (my $i = 0; $i < @percents; $i++) {
	push(@filenames, "$outdir/$outfile_basename\_perc_$percents[$i].txt");
	push(@domainfiles,"");
}


while(<TMPL>) {
	my $line = $_;
	chomp $line;

	if ($line !~ /MARKER_PERC/) {
		if ($line =~ /SAVE_FILENAME/) {
			my ($param, $value) = split(/\t/, $line);
			$value =~ s/\s/_/g;
			for (my $i = 0; $i < @percents; $i++) {
				$domainfiles[$i] .= "$param\t$value\_$percents[$i]\n";
			}
		}
		else {
			for (my $i = 0; $i < @percents; $i++) {
				$domainfiles[$i] .= "$line\n";
			}
		}
	}
	else {
		for (my $i = 0; $i < @percents; $i++) {
			$domainfiles[$i] .= "MARKER_PERC\t$percents[$i]\n";
		}

	}
}

close TMPL;

if (! -e $outdir) {
	`mkdir -p $outdir`;
}

for (my $i = 0; $i < @filenames; $i++) {
	open(OUT, ">", $filenames[$i]);
	print OUT $domainfiles[$i];
	close OUT;
}
