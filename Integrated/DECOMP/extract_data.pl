

use strict;
use IO::Handle;
use POSIX;
use List::Util qw (min max);
use Getopt::Std;


# Inputs:
# -i <input data file>
# -t <datatype to extract>
# -n <take only this top number of edges of the datatype. Default: take all>
# -o <output PPI filename>
my %argopts;
if (! getopts('i:t:n:o:', \%argopts)) { die "invalid arguments"; }
my $infilename = $argopts{'i'};
my $outputfilename = $argopts{'o'};
my $datatype = $argopts{'t'};
my $num_top_edges = -1;
if (defined $argopts{'n'}) {
	$num_top_edges = $argopts{'n'}+0;
}

		
open (INFILE, $infilename) || die "cannot open input file";
open (OUTPUTFILE, ">$outputfilename") || die "cannot open output file";


my %edges = (); # $edges{$key} = $score
foreach my $line (<INFILE>) {
	chomp $line;
	my @toks = split(/\s/, $line);
	my $prot1 = $toks[0];
	my $prot2 = $toks[1];
	my $score = 1;
	if (scalar @toks == 4) {
		if ($toks[2] ne $datatype) { next; }
		$score = $toks[3]+0;
	}	
	elsif (scalar @toks == 3) {
		$score = $toks[2]+0;
	}
	if ($prot1 eq $prot2) { next; }
	my $key = $prot1 lt $prot2 ? "$prot1|$prot2" : "$prot2|$prot1";
	$edges{$key} = $score;
}
print "Num total $datatype edges = ".(scalar keys %edges)."\n";
	
if ($num_top_edges != -1) {
	my @sorted_edges = sort {$edges{$b} <=> $edges{$a}} keys %edges;
	for (my $idx = $num_top_edges; $idx < scalar @sorted_edges; $idx++) {
		delete $edges{$sorted_edges[$idx]};
	}
	print "Keeping top $num_top_edges $datatype edges, num $datatype edges = ".(scalar keys %edges)."\n";
}

# print resultant PPI network
foreach my $edge (sort {$edges{$b} <=> $edges{$a}} keys %edges) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	my $score = $edges{$edge};
	print OUTPUTFILE "$prot1\t$prot2\t$score\n";
}
	