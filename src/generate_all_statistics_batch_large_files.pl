#!/usr/bin/env perl
use warnings;
use strict;

my ($statedir, $outdir, $outscrpt, $polymer_len, @domains) = @ARGV;
die "usage: calculate_all_statistics.pl <dir of all state list files> <output dir for stats files> <output batch script name> <polymer len> <domain str>\n" unless @ARGV >= 4;

my $domainstr = "";
foreach my $domain (@domains) {
	$domainstr .= " $domain";
}

opendir(DIR, $statedir) or die "could not open dir $statedir: $!\n";

my @files = grep { /.*state_list.*/ } readdir DIR;

closedir(DIR);



foreach my $i (@files) {
	#print "On file [$i]\n";
	open(OUT, ">", $outscrpt . "_" . $i . "_batch.sh") or die "could not open $outscrpt: $!\n";

	my $header = "#!/bin/bash
#SBATCH --job-name=stats_$i
#SBATCH --output=stats_$i.%j.out
#SBATCH --error=stats_$i.%j.err
#SBATCH --time=8:00:00
#SBATCH -p normal
#SBATCH --mem=32000
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=arrajpur\@stanford.edu

";

	print OUT $header;
	my ($outname) = $i =~ /\/*([^\/]+)_state_list.*/;
	#warn "outname is [$outname]";

	#warn "command is [calculate_polymer_statistics.py $statedir/$i $outdir/$outname]";
	print OUT "calculate_polymer_statistics_large_files.py $statedir/$i $polymer_len $outdir/$outname $domainstr\n";
	close OUT;
}
