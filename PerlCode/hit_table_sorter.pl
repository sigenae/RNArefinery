#!/usr/bin/perl
use strict;
use warnings;

my $lastname = "";
my @records;

while (<STDIN>) {

	#
	# clean the artifacts and skip the empty lines
	chomp;
	if ( $_ eq "" ) {
		next;
	}

	# we expect nine columns in here
	#
	my @columns = split( '\t', $_ );

	if ( 9 != @columns ) {
		print "Unable to find 9 columns in the line \"" 
		  . $_
		  . " \".. found only "
		  . @columns . "\n";
	}

	if ( length($lastname) == 0 ) {    # we just entered the loop
		$lastname = $columns[0];
		push @records, [@columns];
	}
	else {

		if ( $lastname eq $columns[0] )
		{    # check if we still processing the same query
			    # filter this hit
			my $identity  = $columns[3];
			my $bits      = $columns[4];
			my $query_len = $columns[5];
			my $hit_len   = $columns[6];
			my $gaps      = $columns[7];
			my $evalue    = $columns[8];
			if ( $evalue > 0.0001 ) {
				push @records, [@columns];
			}
		}
		else {

			# if name is not the same sorting an array
			# using bit score
			#
			my @sorted =
			  sort { $b->[4] <=> $a->[4] } @records;

			#sort { $b->[4] <=> $a->[4] || $b->[7] <=> $a->[7] } @records;

			print "$sorted[0][0]\t$sorted[0][1]\t$sorted[0][2]\t$sorted[0][3]\t"
			  . "$sorted[0][4]\t$sorted[0][5]\t$sorted[0][6]\t$sorted[0][7]\t$sorted[0][8]\n";

			# and re-init the array
			@records = [@columns];

			$lastname = $columns[0];
		}
	}
}

#
# and make sure the last table was sorted too
my @sorted = sort { $b->[4] <=> $a->[4] } @records;
print "$sorted[0][0]\t$sorted[0][1]\t$sorted[0][2]\t$sorted[0][3]\t"
  . "$sorted[0][4]\t$sorted[0][5]\t$sorted[0][6]\t$sorted[0][7]\t$sorted[0][8]\n";
