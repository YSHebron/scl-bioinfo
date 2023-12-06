

use strict;
use IO::Handle;
use POSIX;
use List::Util qw (min max);
use Getopt::Std;


# Inputs:
# -i <input PPI file>
# -p <take only this top number of PPI edges. Default: take all>
# -n <Nhub: minimum number of neighbours a protein has to be considered a hub>
# -h <output hubs file>
# -o <output PPI filename>
my %argopts;
if (! getopts('i:p:n:h:o:', \%argopts)) { die "invalid arguments"; }
my $ppifilename = $argopts{'i'};
my $outputfilename = $argopts{'o'};
my $Nhub = $argopts{'n'} + 0;
my $outputhubsfile = $argopts{'h'};
my $num_top_ppis = -1;
if (defined $argopts{'p'}) {
	$num_top_ppis = $argopts{'p'}+0;
}

		
open (PPIFILE, $ppifilename) || die "cannot open PPI file";
open (OUTPUTFILE, ">$outputfilename") || die "cannot open output file";


my %ppi_edges = (); # $ppi_edges{$key} = $score
foreach my $line (<PPIFILE>) {
	chomp $line;
	my @toks = split(/\s/, $line);
	my $prot1 = $toks[0];
	my $prot2 = $toks[1];
	my $score = 1;
	if (scalar @toks == 4) {
		if ($toks[2] ne "PPIREL") { next; }
		$score = $toks[3]+0;
	}	
	elsif (scalar @toks == 3) {
		$score = $toks[2]+0;
	}
	if ($prot1 eq $prot2) { next; }
	my $key = $prot1 lt $prot2 ? "$prot1|$prot2" : "$prot2|$prot1";
	$ppi_edges{$key} = $score;
}
print "Num PPI edges = ".(scalar keys %ppi_edges)."\n";
	
if ($num_top_ppis != -1) {
	my @sorted_ppis = sort {$ppi_edges{$b} <=> $ppi_edges{$a}} keys %ppi_edges;
	for (my $idx = $num_top_ppis; $idx < scalar @sorted_ppis; $idx++) {
		delete $ppi_edges{$sorted_ppis[$idx]};
	}
	print "Keeping top $num_top_ppis PPIs, num PPI edges = ".(scalar keys %ppi_edges)."\n";
}
#print "PPIs:\n";
#foreach my $ppi (sort {$ppi_edges{$b} <=> $ppi_edges{$a}} keys %ppi_edges) {
#	print "$ppi\t$ppi_edges{$ppi}\n";
#}
	
# build %neighbours
my %neighbours = (); # $neighbours{$prot1}{$prot2} = 1
foreach my $ppi (keys %ppi_edges) {
	(my $prot1, my $prot2) = split(/\|/, $ppi);
	$neighbours{$prot1}{$prot2} = 1;
	$neighbours{$prot2}{$prot1} = 1;	
}
	
# identify hubs
my %hubs = (); 
foreach my $prot (keys %neighbours) {
	if (scalar keys %{$neighbours{$prot}} >= $Nhub) {
		$hubs{$prot} = scalar keys %{$neighbours{$prot}};
	}
}
print "Hub threshold = $Nhub, num hubs = ".(scalar keys %hubs)."\n";

# $degrees{$d} = num prots with degree $d
my %degrees = ();
foreach my $prot (keys %neighbours) {
	my $deg = scalar keys %{$neighbours{$prot}};
	$degrees{$deg}++;
}
print "Degree distribution:\n";
foreach my $deg (sort {$a <=> $b} keys %degrees) {
	print "$deg\t$degrees{$deg}\n";
}


# remove hubs from PPI network
my %hub_edges = ();
foreach my $edge (keys %ppi_edges) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	if (defined $hubs{$prot1} || defined $hubs{$prot2}) {
		$hub_edges{$edge} = $ppi_edges{$edge};
	}
}
foreach my $edge (keys %hub_edges) {
	delete $ppi_edges{$edge};
}
print "Num edges removed = ".(scalar keys %hub_edges).", num PPI edges remaining = ".(scalar keys %ppi_edges)."\n";


# print resultant PPI network
foreach my $edge (keys %ppi_edges) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	my $score = $ppi_edges{$edge};
	print OUTPUTFILE "$prot1\t$prot2\t$score\n";
}
	
# print hubs
if (defined $outputhubsfile) {
	open (OUTPUTHUBSFILE, ">$outputhubsfile") || die $!;
	print OUTPUTHUBSFILE "hubs\n";
	foreach my $prot (keys %hubs) {
		print OUTPUTHUBSFILE "$prot\t$hubs{$prot}\n";
	}
	print OUTPUTHUBSFILE "edges\n";
	foreach my $edge (keys %hub_edges) {
		print OUTPUTHUBSFILE "$edge\t$hub_edges{$edge}\n";
	}
	
	
	
}