#!/usr/bin/perl

use strict;
use warnings;



#######################################################
#                  sgens.pl                           #
#  calculates the numbers of generations required     #
#  for frequency shift given s selection coefficient  #
#######################################################

my $q0 = $ARGV[0];
my $q1 = $ARGV[1];
my $s = $ARGV[2];

unless ($q0 and $q1 and $s) {
	die "\n\nERROR: Not enough arguments.\n\n";
}


my $g = 0; #number of generations

until ($q0 >= $q1) {
	$q0 = q_prime($q0, $s);
	++$g; #count generations
}

print "\n\n$g\n\n";





exit;

sub q_prime { # get frequency of q allele in next generation
	my $q = shift;
	my $s = shift;
	my $p = 1- $q;
	
	my $q1 = (($p * $q) + ($q ** 2) + ($s * ($q ** 2))) / (1 + ($s * ($q ** 2)));
	return $q1;
}

