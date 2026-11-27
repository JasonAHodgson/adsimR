#!/usr/bin/perl

#use strict;
use warnings;


############################################
#           fastadsim_cntrl.pl             #
#  runs for fastadsim.pl for all combos    #
#           of given conditions            #
############################################


my $conditions = $ARGV[0]; #file containing conditions to simulate
my $results = $ARGV[1]; #results directory. this directory has to exist
chomp $conditions;
if ($results) {chomp $results;}

unless ($conditions) {
	die "\n\nERROR: You must specify a conditions file.\n\n"; 
}

#unless ($s) {
#	$s = 1000; #default number of simulations
#}

#unless ($s =~ /^\d+$/) {
#	die "\n\nERROR: You have indicated $s simulations to run.  This must be a positive integer.\n\n";
#}

unless (-f $conditions) {
	die "\n\nERROR: Cannot find $conditions. Check file path.\n\n";
}




# -g number of generations of admixture
# -k population growth parameter
# -l number of generations prior to population growth
# -a desired final admixture proportion of population 1
# -m1 migration rate population 1
# -m2 migration rate population 2
# -n effective population size of generation 0
# -p1 allele 1 frequency of population 1
# -p2 allele 1 frequency of population 2
# -s number of simulations to run



open COND, "$conditions"; #open conditions files to get conditions

my @g; #number of generations of admixture
my @k; #population growth parameters
my @l; #number of generations prior to population growth
my @a; #admixture proportion of population 1
my @n; #effective population size of generation 0
my @m1; #migration rate population 1
my @m2; #migration rate population 2
my @p1; #allele frequency of population 1
my @p2; #allele frequency of population 2
my @s; #number of simulations to run for each parameter set


#######################################################
#                 Make a log file                     #
#######################################################

my $log = "fastadmigsim_cntrl.log"; #logfile name

open LOG, ">$log"; #open the logfile to write to

print LOG "fastadmigsim_cntrl.pl logfile\n";

my $datetime; #date and time program is run

my $time = localtime(time); #get the time data

#the following converts the time data into a useful readout

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my  @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
(my $second, my $minute, my $hour, my $dayOfMonth, my $month, my $yearOffset, my $dayOfWeek, my $dayOfYear, my $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

print LOG "$theTime\n\n"; #add date and time data to the log file




while (<COND>) { #read conditions file to get conditions
	my $line = $_; #save line
	chomp $line;
	my $flag; #indicated condition
	my $c; #conditions
	if ($line =~ /^(-\w{1,2}):\s+(.+)$/) { #look for flag and conditions
		$flag = $1;
		$c = $2;
		chomp $c;
		print "flag: $flag\n";
		print "conditions: $c\n";
		my @cond = split " ", $c; #split into conditions
		if ($flag eq "-g") { #if generations
			@g = @cond;
		} elsif ($flag eq "-k") {
			@k = @cond;
		} elsif ($flag eq "-l") {
			@l = @cond;
		} elsif ($flag eq "-a") {
			@a = @cond;
		} elsif ($flag eq  "-m1") {
			@m1 = @cond;
		} elsif ($flag eq "-m2") {
			@m2 = @cond;
		} elsif ($flag eq "-n") {
			@n = @cond;
		} elsif ($flag eq "-p1") {
			@p1 = @cond;
		} elsif ($flag eq "-p2") {
			@p2 = @cond;
		} elsif ($flag eq "-s") {
			@s = @cond;
		} else {
			die "\n\nERROR: $flag is not recognized.\n\n";
		}
	} 
}

close COND; #close conditions file

my @constant; #constant parameters
my @variable; #variable parameters

#separate into variable and constant parameters
if (scalar @g == 1) { 
	push @constant, "-g";
} elsif (scalar @g > 1) {
	push @variable, "-g";
} else {
	die "\n\nERROR: No -g specified.\n\n";
}

if (scalar @k == 1) { 
	push @constant, "-k";
} elsif (scalar @k > 1) {
	push @variable, "-k";
} else {
	die "\n\nERROR: No -k specified.\n\n";
}

if (scalar @l == 1) { 
	push @constant, "-l";
} elsif (scalar @l > 1) {
	push @variable, "-l";
} else {
	die "\n\nERROR: No -l specified.\n\n";
}

if (scalar @a == 1) { 
	push @constant, "-a";
} elsif (scalar @a > 1) {
	push @variable, "-a";
} else {
	die "\n\nERROR: No -a specified.\n\n";
}

if (scalar @m1 == 1) { 
	push @constant, "-m1";
} elsif (scalar @m1 > 1) {
	push @variable, "-m1";
} else {
	die "\n\nERROR: No -m1 specified.\n\n";
}

if (scalar @m2 == 1) { 
	push @constant, "-m2";
} elsif (scalar @m2 > 1) {
	push @variable, "-m2";
} else {
	die "\n\nERROR: No -m2 specified.\n\n";
}

if (scalar @n == 1) { 
	push @constant, "-n";
} elsif (scalar @n > 1) {
	push @variable, "-n";
} else {
	die "\n\nERROR: No -n specified.\n\n";
}

if (scalar @p1 == 1) { 
	push @constant, "-p1";
} elsif (scalar @p1 > 1) {
	push @variable, "-p1";
} else {
	die "\n\nERROR: No -p1 specified.\n\n";
}

if (scalar @p2 == 1) { 
	push @constant, "-p2";
} elsif (scalar @p2 > 1) {
	push @variable, "-p2";
} else {
	die "\n\nERROR: No -p2 specified.\n\n";
}

if (scalar @s == 1) { 
	push @constant, "-s";
} elsif (scalar @s > 1) {
	push @variable, "-s";
} else {
	die "\n\nERROR: No -s specified.\n\n";
}

print "\nThe following parameters are constant:\n"; #report to screen
print LOG "\nThe following parameters are constant:\n"; #report to screen
foreach (@constant) {
	print "$_\n";
	print LOG "$_\n";
}

print "\nThe following parameters are variable:\n"; #report to screen
print LOG "\nThe following parameters are variable:\n"; #report to screen


foreach (@variable) {
	print "$_\n";
	print LOG "$_\n";
}

my $combos = scalar @a * scalar @n * scalar @g * scalar @k * scalar @l * scalar @p1 * scalar @p2 * scalar @s * scalar @m1 * scalar @m2;

print "Number of parameter combinations to run: $combos\n";
print LOG "Number of parameter combinations to run: $combos\n";

my @all = ([@a], [@n], [@g], [@k], [@l], [@p1], [@p2], [@s], [@m1], [@m2]); #array of arrays

my @d; #all combinations of parameters in strings
push @d, "@{[ /.+/g ]}" for glob join ",", map { local $" = ","; "{@$_}" } @all; #no idea how this works, but it prints lines of all combinations of parameters


#my @args; #arguments to pass to fastadsim
foreach  (@d) {
	my @l = split ",", $_;
	
	#my $a10 = get_a10($l[1], $l[0], $l[3], $l[2], $l[8]);
	
	#my $a10 = $l[0];
	
	my $out = "a-$l[0]\_n-$l[1]_g-$l[2]_k-$l[3]_l-$l[4]_p1-$l[5]_p2-$l[6]_s-$l[7]_m1-$l[8]_m2-$l[9]"; #outfile and directory name
	my $argline = "-a $l[0] -n $l[1] -g $l[2] -k $l[3] -l $l[4] -p1 $l[5] -p2 $l[6] -s $l[7] -m1 $l[8] -m2 $l[9] -o $out";
	print "Running fastmadmigsim.pl with the following argument:\n$argline\n";
	print "Writing results to: $out\_results\n";
	print LOG "Running fastadmigsim.pl with the following argument:\n$argline\n";
	print LOG "Writing results to: $out\_results\n";
	mkdir "$out\_results";
	system "./scripts/fastadmigsim.pl $argline";
	system "./scripts/adsim2R.pl $out";
	system "./scripts/adsimPlot.R $out\.Rtable $out\.plot.pdf";
	system "mv $out $out\_results/";
	system "mv $out\.Rtable $out\_results/";
	system "mv $out\.plot.pdf $out\_results/";
	if (-d $results) {
		system "mv $out\_results $results";
	}
	#push @args, $argline;
}

print "Analyses finished.\n";
print LOG "Analyses finished.\n";



############################################################
#                      SUBROUTINES                         #
############################################################


sub get_a10 {
	
	my $ne0 = shift; # Ne at time 0
	my $a1g = shift; #  admixture proportion after g generations
	my $k = shift; #   Population growth parameter
	my $g = shift; #   number of generations
	my $m = shift; # migrants per generation
	
	my $a10 = $a1g - (($m * $g * (1 - $a1g)) / ($k**$g * $ne0)); #calculate value
	 
	return $a10; # return value
	
}



exit;