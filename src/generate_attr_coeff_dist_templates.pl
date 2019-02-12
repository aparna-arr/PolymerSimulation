#!/usr/bin/env perl
use warnings;
use strict;

my ($template_file, $xyz_file, $output_dir, $output_prefix, $reps, @attr_coeffs) = @ARGV;
die "usage: generate_attr_coeff_dist_templates.pl <path/to/template parameter file> <path/to/xyz file> <output directory> <outfile prefix> <num_reps> <attr_coeff 1> <attr_coeff 2> <...>

Note the attraction coefficients will apply to ANY PARTICLE INTERACTION WITH ATTRACT_FORCE SET TO TRUE.\n" unless @ARGV >= 6;

if (! -e $output_dir) {
	`mkdir -p $output_dir`;
}

`mkdir -p $output_dir/templates/`;
`mkdir -p $output_dir/inputs/`;
`mkdir -p $output_dir/batch/`;
`mkdir -p $output_dir/statistics/`;
`cp $template_file $output_dir/templates/`;

my @templatenames;
my @paramlines;

my ($basename) = $xyz_file =~ /\/*([^\/]+)\.xyz$/;
foreach my $coeff (@attr_coeffs) {
	for (my $i = 0; $i < $reps; $i++) {
		push(@paramlines, "");
		push(@templatenames, "$output_dir/inputs/$output_prefix\_$basename\_attr_$coeff\_rep_$i\.txt");
		
	}
}

open(TEMPL, "<", $template_file) or die "could not open $template_file:$!\n";

while(<TEMPL>) {
	my $line = $_;
	chomp $line;

	if ($line !~ /^INPUT_POLYMER_FILE|^SAVE_PATH|^SAVE_FILENAME|^ATTR_E.+\$INIT/) {
	
		for (my $i = 0; $i < @templatenames; $i++) {
			$paramlines[$i] .= $line . "\n";
		}
	}
        else {
                if ($line =~ /INPUT_POLYMER_FILE/) {
                        for (my $i = 0; $i < @templatenames; $i++) {
                                $paramlines[$i] .= "INPUT_POLYMER_FILE\t$xyz_file\n";
                        }
                }
                elsif ($line =~ /SAVE_PATH/) {
                        for (my $i = 0; $i < @templatenames; $i++) {
                                my ($basepath) = $output_dir . "/simulation";
                                if ($basepath eq "") {
                                        $basepath = ".";
                                }

                                $paramlines[$i] .= "SAVE_PATH\t$basepath\n";
                        }
                }
                elsif ($line =~ /SAVE_FILENAME/) {
                        for (my $i = 0; $i < @templatenames; $i++) {
                                my ($basename) = $templatenames[$i] =~ /\/*([^\/]+)\.txt$/;
                                $paramlines[$i] .= "SAVE_FILENAME\tsim_$basename\n";
                        }

                }
		elsif ($line =~ /ATTR_E.+\$INIT/) {
			my $count = 0;
			for (my $i = 0; $i < @attr_coeffs; $i++) {
				my $attre = $line;
				$attre =~ s/\$INIT/$attr_coeffs[$i]/;
				for (my $j = 0; $j < $reps; $j++) {
					$paramlines[$count] .= $attre . "\n";
					$count++;
				}
			}
		}
        }
}

close TEMPL;


my $idx = 0;
for (my $i = 0; $i < @attr_coeffs; $i++) {

	my $batch_template = "#!/bin/bash
#SBATCH --job-name=$output_prefix\_attre_$attr_coeffs[$i]
#SBATCH --output=$output_prefix\_attre_$attr_coeffs[$i].%j.out
#SBATCH --error=$output_prefix\_attre_$attr_coeffs[$i].%j.err
#SBATCH --time=1-00:00:00
#SBATCH -p gpu
#SBATCH --gres gpu:1
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=arrajpur\@stanford.edu

module load cuda\n";

	open(BAT, ">", "$output_dir/batch/run_batch_coeff_$attr_coeffs[$i]\_$output_prefix\.sh");
	my $curr_batch_template = $batch_template;
	for (my $j = 0; $j < $reps ; $j++) {
		
        	open(TMP, ">", $templatenames[$idx]);
        	print TMP $paramlines[$idx];
        	close TMP;

        	$curr_batch_template .= "simulate_polymer.py $templatenames[$idx]\n";
		$idx += 1;
	}
	print BAT $curr_batch_template;
	close BAT;
}
