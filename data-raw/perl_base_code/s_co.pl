#!/usr/bin/perl

use strict;
use warnings;


##############################################################
#                         s_co.pl                            #
#  Returns selection coefficients for a table of scenarios   #
##############################################################



my $in = $ARGV[0];  #table containing allele frequency change scenarios
my $out = $ARGV[1]; #table with selection coefficients added


unless ($in and $out) {
	die "\n\nERROR: Not enough arguments.\n\n";
}

unless (-f $in) {
	die "\n\nERROR: Cannot find $in. Check file path.\n\n";
}


open IN, "$in";

my $header = <IN>;
chomp $header;

open OUT, ">$out";

$header = $header."\ts\n";

print OUT $header;

while (<IN>) {
	my $line = $_;
	chomp $line;
	my @l = split /\s/, $line;
	my $g = $l[0];
	my $q = $l[1];
	my $qx = $l[2];
	my $s = get_s($q, $qx, $g);
	print OUT "$g\t$q\t$qx\t$s\n";
}

close IN;
close OUT;




exit;


sub q_prime { # get frequency of q allele in next generation
	my $q = shift;
	my $s = shift;
	my $p = 1- $q;
	
	my $q1 = (($p * $q) + ($q ** 2) + ($s * ($q ** 2))) / (1 + ($s * ($q ** 2)));
	return $q1;


}

sub get_qx { #get final allele frequency given s and number of generations
	my $q = shift;
	my $g = shift;
	my $s = shift;
	
	my $p = 1 - $q;
	
	my $i = 0;
	
	until ($i == $g) {
		$q = q_prime($q, $s);
		++$i;
	}
	return $q;
	

}


sub get_s { #get the best s given initial q, final q and number of generations

	my $q = shift;
	my $qx = shift;
	my $g = shift;

	my @s = map{ $_/10000 } 1 .. 9999; #list of possible s to check
	

	my $final_q = $q;

	my $i = 0;
	my $s_test;
	until ($final_q == $qx or $i == scalar @s) {  #look for s that matches
		$s_test = $s[$i];
		my $testq = get_qx($q, $g, $s_test);
		$final_q = sprintf "%.2f", $testq;	
		++$i; 
	}

	return $s_test;
}