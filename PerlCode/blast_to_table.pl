#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;
use Bio::SearchIO;
use Bio::Search::Tiling::MapTiling;

## configure input -- STDIN
#
my $input = Bio::SearchIO->new(
	-format => 'blast',
	-fh     => \*ARGV,
);

#my $input = Bio::SearchIO->new(
#	-format => 'blast',
#	-file   => 'blast.out.txt'
#);

# result is the next blast search result, i.e. Query
my $query_counter = 0;
while ( ( my $result = $input->next_result ) ) {

	# $hit is a Bio::Search::Hit::HitI compliant object, iterate over these
	while ( my $hit = $result->next_hit ) {

		#		my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
		#		my @contextsQ = $tiling->contexts('query');
		#		for my $contextQ (@contextsQ) {
		#			my ( $min, $max ) = $tiling->range( 'query', $contextQ );
		#			my $tLengthQ = $tiling->length( 'query', 'exact', $contextQ );
		#			print
		#			  "QUERY range in *$contextQ* => $min, $max (length=$tLengthQ)\n";
		#		}
		#		my @contextsS = $tiling->contexts('subject');
		#		for my $contextS (@contextsS) {
		#			my ( $min2, $max2 ) = $tiling->range( 'subject', $contextS );
		#			print "SUBJECT range in *$contextS* => $min2, $max2\n";
		#		}

		# $hsp is a Bio::Search::HSP::HSPI compliant object
		while ( my $hsp = $hit->next_hsp ) {

			#my $rounded_identity = floor( $hsp->percent_identity );
			print $result->query_name . "\t"     #1 query name
			  . $hit->name . "\t"                #2 hit name
			  . $hit->description . "\t"         #3 hit description
			  . $hsp->percent_identity . "\t"    #4 identity
			  . $hsp->bits . "\t"                #5 bits
			  . $result->query_length . "\t"     #6 query length
			  . $hsp->length('hit') . "\t"       #7 alignment length
			  . $hsp->gaps . "\t"                #8 gaps
			  . $hsp->evalue . "\t" . "\n";      #9 evalue

		}
	}

	$query_counter = $query_counter + 1;
}

print STDERR "processed "
  . $query_counter
  . " sequence alignmentes\n";
