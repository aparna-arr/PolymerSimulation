#!/usr/bin/env perl
use warnings;
use strict;

my ($template_file, $xyz_contents, $output_dir, $output_prefix) = @ARGV;
die "usage: generate_index_files.pl <path/to/template parameter file> <path/to/xyz contents file> <output directory> <outfile prefix>\n" unless @ARGV == 4;

if (! -e $output_dir) {
	`mkdir -p $output_dir`;
}

`mkdir -p $output_dir/templates/`;
`mkdir -p $output_dir/inputs/`;
`mkdir -p $output_dir/batch/`;
`cp $template_file $output_dir/templates/`;
`cp $xyz_contents $output_dir/templates/`;

open(CONT, "<", $xyz_contents) or die "could not open $xyz_contents: $!\n";

my @xyz_filenames;
my @paramlines;
my @templatenames;

while(<CONT>) {
	my $line = $_;
	chomp $line;

	#warn "Line is [$line]";
	
	next if $line !~ /\.xyz$/;
	
	my ($basename) = $line =~ /\/*([^\/]+)\.xyz$/;

	#warn "FIRST basename is [$basename]";

	push(@xyz_filenames, $line);
	push(@paramlines, "");
	push(@templatenames, "$output_dir/inputs/$output_prefix\_$basename\.txt");
}

close CONT;

open(TEMPL, "<", $template_file) or die "could not open $template_file:$!\n";

while(<TEMPL>) {
	my $line = $_;
	chomp $line;

	if ($line !~ /^INPUT_POLYMER_FILE|^SAVE_PATH|^SAVE_FILENAME/) {
		for (my $i = 0; $i < @xyz_filenames; $i++) {
			$paramlines[$i] .= $line . "\n";
		}
	}
	else {
		if ($line =~ /INPUT_POLYMER_FILE/) {
			for (my $i = 0; $i < @xyz_filenames; $i++) {
				$paramlines[$i] .= "INPUT_POLYMER_FILE\t$xyz_filenames[$i]\n";
			}
		}
		elsif ($line =~ /SAVE_PATH/) {
			for (my $i = 0; $i < @xyz_filenames; $i++) {
				my ($basepath) = $output_dir . "/simulation";
				#warn "1 basepath is $basepath\n";

				if ($basepath eq "") {
					$basepath = ".";
				}

				#warn "2 basepath is $basepath\n";
				$paramlines[$i] .= "SAVE_PATH\t$basepath\n";
			}

		}
		elsif ($line =~ /SAVE_FILENAME/) {
			for (my $i = 0; $i < @xyz_filenames; $i++) {
				my ($basename) = $xyz_filenames[$i] =~ /\/*([^\/]+)\.xyz$/;
				$paramlines[$i] .= "SAVE_FILENAME\tsim_$basename\n";
			}

		}
	}
	
}

close TEMPL;

# now make batch file
my $batch_template = "#!/bin/bash
#SBATCH --job-name=$output_prefix
#SBATCH --output=$output_prefix.%j.out
#SBATCH --error=$output_prefix.%j.err
#SBATCH --time=5:00:00
#SBATCH -p gpu
#SBATCH --gres gpu:1
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=arrajpur\@stanford.edu

module load cuda

";

for (my $i = 0; $i < @templatenames; $i++) {
	open(TMP, ">", $templatenames[$i]);
	
	print TMP $paramlines[$i];
		
	close TMP;

	$batch_template .= "python simulate_polymer.py $templatenames[$i]\n";
}

open(BAT, ">", "$output_dir/batch/run_batch_$output_prefix\.sh");
print BAT $batch_template;
close BAT;



