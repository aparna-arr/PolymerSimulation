#!/usr/bin/env perl
use warnings;
use strict;

my ($statefile, $xyzfile) = @ARGV;
die "usage: $0 <state file> <xyz file with domain labels>\n" unless @ARGV == 2;

open(XYZ, "<", $xyzfile) or die "couldn't open $xyzfile: $!\n";

my $xyz_state_str = "";
my @domain_labels;
while(<XYZ>) {
	my $line = $_;
	chomp $line;

	if ($line =~ /^\d+$/) {
		$xyz_state_str .= $line . "\n";
	}
	elsif ($line =~ /^\s*$/){
		$xyz_state_str .= "\n";
	}
	else {
		my (@lineAr) = split(/\t/, $line);
		my $label = $lineAr[0];
		push(@domain_labels,$label);
	}
}

close XYZ;

open(STATE, "<", $statefile) or die "couldn't open $statefile: $!\n";

my @xyzpos;

while(<STATE>) {
	my $line = $_;
	chomp $line;
	
	if ($line =~ /\s*<Position x.+/) {
		my ($x, $y, $z) = $line =~ /.*x="(.+)"\sy="(.+)"\sz="(.+)".*/;
		
		push(@xyzpos, "$x\t$y\t$z");		
	}	
}
close STATE;

open(XYZNEW, ">", $statefile . ".xyz");

print XYZNEW $xyz_state_str;

for (my $i = 0; $i < @domain_labels; $i++) {
	my $lab = $domain_labels[$i];
	my $xyz = $xyzpos[$i];

	print XYZNEW "$lab\t$xyz\n";
}

close XYZNEW;
