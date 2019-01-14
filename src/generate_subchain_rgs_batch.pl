#!/usr/bin/env perl
use warnings;
use strict;

my ($seed, $num_reps, $polymer_len, $list_of_subchain_lens, $list_of_state_files, $outname_script, $outfile_state) = @ARGV;
die "usage: generate_subchain_rgs_batch.pl <seed> <number of repetitions for random draws> <total_polymer_length> <list_of_subchain_lengths.txt> <list_of_state_files.txt> <outname for script> <path/to/outfile/outfile_base_name>" unless @ARGV == 7;

srand($seed);

open(SUB, "<", $list_of_subchain_lens) or die "could not open $list_of_subchain_lens: $!\n";

my @subchain_lens;

while(<SUB>) {
	my $line = $_;
	chomp $line;

	# this will also crash if we're handling a non-number
	if ($line >= $polymer_len) {
		warn "subchain length [$line] is >= total polymer length [$polymer_len]! Exiting.";
		close SUB;
		exit(1);
	}
	
	push(@subchain_lens, $line);
}

close SUB;


for (my $i = 0; $i < $num_reps; $i++) {
	my $header = "#!/bin/bash
#SBATCH --job-name=generate_subchain_stats_random_rep_$i
#SBATCH --output=generate_subchain_stats_random_rep_$i.%j.out
#SBATCH --error=generate_subchain_stats_random_rep_$i.%j.err
#SBATCH --time=6:00:00
#SBATCH -p normal
#SBATCH --mem=16000
#SBATCH --nodes=1
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=arrajpur\@stanford.edu

";

	open(OUT, ">", $outname_script . "_rep_$i.sh") or die "could not open $outname_script\_rep_$i.sh: $!\n";
	print OUT $header;


	my $domain_pairs_str = "";
	foreach my $sublen (@subchain_lens) {
		my $max_end = $polymer_len - $sublen;

		my $sub_start = int(rand($max_end));
		my $sub_end = $sub_start + $sublen;

		$domain_pairs_str .= " $sub_start,$sub_end";
	}
	print OUT "calculate_polymer_statistics_large_files.py $list_of_state_files $polymer_len $outfile_state\_random_subchains_rep_$i" . "$domain_pairs_str\n";

	close OUT;
}
