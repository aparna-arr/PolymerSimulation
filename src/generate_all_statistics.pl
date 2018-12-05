#!/usr/bin/env perl
use warnings;
use strict;

my ($statedir, $outdir, $outscrpt) = @ARGV;
die "usage: calculate_all_statistics.pl <dir of all state list files> <output dir> <output batch script name>\n" unless @ARGV == 3;

opendir(DIR, $statedir) or die "could not open dir $statedir: $!\n";

my @files = grep { /.*state_list.*/ } readdir DIR;

closedir(DIR);

open(OUT, ">", $outscrpt) or die "could not open $outscrpt: $!\n";

my $header = "#!/bin/bash
#SBATCH --job-name=generate_statistics
#SBATCH --output=generate_statistics.%j.out
#SBATCH --error=generate_statistics.%j.err
#SBATCH --time=8:00:00
#SBATCH -p normal
#SBATCH --mem=64000
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=arrajpur\@stanford.edu

";

print OUT $header;

foreach my $i (@files) {
	#print "On file [$i]\n";
	my ($outname) = $i =~ /\/*([^\/]+)_state_list.*/;
	#warn "outname is [$outname]";

	#warn "command is [calculate_polymer_statistics.py $statedir/$i $outdir/$outname]";
	print OUT "calculate_polymer_statistics.py $statedir/$i $outdir/$outname\_statistics.txt\n";
}
close OUT;
