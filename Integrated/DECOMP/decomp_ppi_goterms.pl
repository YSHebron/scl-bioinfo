#!/usr/bin/perl -w
use strict;
use IO::Handle;
use Getopt::Std;



# inputs:
# -i <input PPI file>
# -a <Input GO annotations file>
# -d <GO decomposition terms>
# -o <output basefile>
# perl decomp_ppi_goterms.pl -i "data_ppibiogrid.txt" -a "my_go_associations.sgd_propagated.txt" -d "decompGOterms30.txt" -o "data_ppibiogrid_godecomp30"

my %argopts;
if (! getopts('i:a:d:o:', \%argopts)) { die "invalid arguments"; }



# read decomposition terms
open (DECOMPFILE, "$argopts{'d'}") || die $!;
my %decomp_goterms = ();
foreach my $line (<DECOMPFILE>) {
	chomp $line;
	$decomp_goterms{$line} = 1;
}
	
	
# read GO annotations
my %annots = (); # $annots{$prot}{$goid} = 1. 
open(ANNOTFILE, "$argopts{'a'}") || die $!;
foreach my $line (<ANNOTFILE>) { 
	chomp($line);
  (my $prot, my $goid) = split(/\|/, $line);
	if (!defined $decomp_goterms{$goid}) { next; }
 	$annots{$prot}{$goid} = 1;
}
print "annotations read for ".(scalar keys %annots)." prots\n";


# read ppi data
my %ppi_edges = (); # $ppi_edges{$key} = $score
#my %neighbours = (); # $neighbours{$prot1}{$prot2} = 1
my %ppi_prots = ();
open (PPIFILE, "$argopts{'i'}") || die $!;
foreach my $line (<PPIFILE>) {
	chomp $line;
	my @toks = split(/\s/, $line);
	my $prot1 = $toks[0];
	my $prot2 = $toks[1];
	my $score = 1;
	if (scalar @toks == 4) {
		if ($toks[2] ne "PPI") { next; }
		$score = $toks[3]+0;
	}	
	elsif (scalar @toks == 3) {
		$score = $toks[2]+0;
	}
	if ($prot1 eq $prot2) { next; }
	my $key = $prot1 lt $prot2 ? "$prot1|$prot2" : "$prot2|$prot1";
#	$neighbours{$prot1}{$prot2} = 1;
#	$neighbours{$prot2}{$prot1} = 1;
	$ppi_edges{$key} = $score;
	$ppi_prots{$prot1} = 1;
	$ppi_prots{$prot2} = 1;
}
print "Num PPI edges = ".(scalar keys %ppi_edges).", num proteins = ".(scalar keys %ppi_prots)."\n";
	
	
my %prots_not_discarded = ();
my %edges_not_discarded = ();
my $curr_file_idx = 0;
foreach my $goid (keys %decomp_goterms) {
	my $outputfilename = "$argopts{'o'}"."_".$curr_file_idx.".txt";
	open (OUTPUTFILE, ">$outputfilename") || die $!;
	foreach my $edge (keys %ppi_edges) {
		(my $prot1, my $prot2) = split(/\|/, $edge);
		if (!defined $annots{$prot1}{$goid} || !defined $annots{$prot2}{$goid}) { next; }
		$prots_not_discarded{$prot1} = 1;
		$prots_not_discarded{$prot2} = 1;
		$edges_not_discarded{$edge} = 1;
		print OUTPUTFILE "$prot1\t$prot2\t$ppi_edges{$edge}\n";
	}
	$curr_file_idx++;
	close OUTPUTFILE;
	print "GO term $goid, outputfile $outputfilename\n";
}
print "Num edges used = ".(scalar keys %edges_not_discarded).", num proteins used = ".(scalar keys %prots_not_discarded)."\n";
print "Num edges discarded = ".((scalar keys %ppi_edges) - (scalar keys %edges_not_discarded)).", num proteins discarded = ".((scalar keys %ppi_prots) - (scalar keys %prots_not_discarded))."\n";


