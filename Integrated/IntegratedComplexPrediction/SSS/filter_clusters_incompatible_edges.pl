use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use POSIX;
use Storable qw(dclone);
use Getopt::Std;



# inputs:
# -i <input clusters file>
# -g <GO decomp removed edges file>
# -f <lower edge scores fuzzily (1), or just remove edges (0)? Default 1>
# -o <output clusters file>

my %argopts;
if (! getopts('i:g:f:o:', \%argopts)) { die "invalid arguments"; }


my $FUZZY_EDGE_REMOVE = 1;
if (defined $argopts{'f'}) {
	$FUZZY_EDGE_REMOVE = $argopts{'f'}+0;
}

my %remove_edges = ();
open (REMOVED_EDGES_FILE, $argopts{'g'}) || die $!;
foreach my $line (<REMOVED_EDGES_FILE>) {
	chomp $line;
	(my $edge, my $prob) = split(/\t/, $line);
	if (defined $prob) {
		$remove_edges{$edge} = $prob;
	}
	else {
		$remove_edges{$edge} = 0;
	}
}
print "Num edges to remove read = ".(scalar keys %remove_edges)."\n";



my %clusters = ();
ReadClusters($argopts{'i'}, \%clusters);
print "Num clusters read = ".(scalar keys %clusters)."\n";


## read complexes for debug
## $complexes{$cid}{$pid} = 1 if protein $pid is in complex $cid
#my %complexes = ();
#open (COMPLEXFILE, "complexes_cyc.txt") || die $!;
#foreach my $line (<COMPLEXFILE>) {
#	chomp($line);
#	(my $id, my $comp) = split(/\t/, $line);
#	$id = uc($id);
#	$complexes{$comp}{$id} = 1;
#}
#close COMPLEXES_FILE;
#print "num complexes read = ".(scalar keys %complexes)."\n";
#
#	
## read GO scheme, only for CC
#my %go_scheme = (); # $go_scheme{$goid}{"level"} = level, {"name"} = name, {"parents"}, {"children"}
#open (GOSCHEMEFILE, "go_scheme_all_mine.txt") || die $!;
#foreach my $line (<GOSCHEMEFILE>) {
#	chomp $line;
#	my @tokens = split(' ', $line, 2);
#	my @tokens2 = split(/\|/, $tokens[1]);
#	my $goclass = $tokens2[0];
#	if ($goclass ne "cellular_component") { next; }
#	$go_scheme{$tokens[0]}{"level"} = ($tokens2[1]+0);
#	$go_scheme{$tokens[0]}{"name"} = $tokens2[2];
#	if (scalar @tokens2 == 4) {
#		my @tokens3 = split(',', $tokens2[3]);
#		foreach (@tokens3) {
#			$go_scheme{$tokens[0]}{"parents"}{$_} = 1;
#			$go_scheme{$_}{"children"}{$tokens[0]} = 1;
#		}
#	}
#}
#print "read go scheme, num CC GO terms = ".(scalar keys %go_scheme)."\n";
#	
#	
## read GO annotations, only for CC
#my %annots = (); # $annots{$prot}{$goid} = 1. 
#open(ANNOTFILE, "my_go_associations.sgd_propagated.txt") || die $!;
#foreach my $line (<ANNOTFILE>) { 
#	chomp($line);
#  (my $prot, my $goid) = split(/\|/, $line);
#  if (!defined $go_scheme{$goid}) { next; }
# 	$annots{$prot}{$goid} = 1;
#}
#print "annotations read for ".(scalar keys %annots)." prots\n";
#
#
#
## get the set of most specific annots (those terms with no descendents) for each prot
#my %most_specific_annots = ();
#foreach my $prot (keys %annots) {
#	foreach my $goid (keys %{$annots{$prot}}) {
#		my $child_annotated = 0;
#		foreach my $child (keys %{$go_scheme{$goid}{children}}) {
#			if (defined $annots{$prot}{$child}) {
#				$child_annotated = 1;
#				last;
#			}
#		}
#		if ($child_annotated==0) {
#			$most_specific_annots{$prot}{$goid} = 1;
#		}
#	}
#	if (!defined $most_specific_annots{$prot}) { die; }
#}


#
#my %comp2 = ();
#my %comp3 = ();
#foreach my $comp (keys %complexes) {
#	if (scalar keys %{$complexes{$comp}} == 2) {
#		foreach my $prot (keys %{$complexes{$comp}}) {
#			$comp2{$comp}{$prot} = 1;
#		}		
#	}
#	elsif (scalar keys %{$complexes{$comp}} == 3) {
#		foreach my $prot (keys %{$complexes{$comp}}) {
#			$comp3{$comp}{$prot} = 1;
#		}		
#	}
#}

#my %comp2_removed = ();
#my %comp3_removed = ();
my $num_clus2_orig = 0;
my $num_clus3_orig = 0;
my $num_clus2_removed = 0;
my $num_clus3_removed = 0;
foreach my $clus (sort {$clusters{$b}{SCORE} <=> $clusters{$a}{SCORE}} keys %clusters) {
	my @prots = sort keys %{$clusters{$clus}{E}};
	if (scalar @prots == 2) {
		$num_clus2_orig++;
		my $edge = "$prots[0]|$prots[1]";
		if (defined $remove_edges{$edge}) {
#			my $iscomp = 0;
#			foreach my $comp (keys %comp2) {
#				if (defined $comp2{$comp}{$prots[0]} && defined $comp2{$comp}{$prots[1]}) {
#					print "Removed complex Clus2 $clus, score $clusters{$clus}{SCORE}, is a complex!!\n";
#					$comp2_removed{$clus}{E}{$prots[0]} = 1;
#					$comp2_removed{$clus}{E}{$prots[1]} = 1;
#					$comp2_removed{$clus}{SCORE} = $clusters{$clus}{SCORE};
#					$iscomp = 1;
#				}
#			}
#			if ($iscomp==0) {
#				print "Removed non-complex clus2 $clus, score $clusters{$clus}{SCORE}\n";
#			}
			if ($FUZZY_EDGE_REMOVE==0) {
				delete $clusters{$clus};
			}
			else {
				$clusters{$clus}{SCORE} *= $remove_edges{$edge};
			}			
			$num_clus2_removed++;				
		}
	}
	elsif (scalar @prots == 3) {
		$num_clus3_orig++;
		my $edge1 = "$prots[0]|$prots[1]";
		my $edge2 = "$prots[0]|$prots[2]";
		my $edge3 = "$prots[1]|$prots[2]";
		if (defined $remove_edges{$edge1} || defined $remove_edges{$edge2} || defined $remove_edges{$edge3}) {
#			my $iscomp = 0;
#			foreach my $comp (keys %comp3) {
#				if (defined $comp3{$comp}{$prots[0]} && defined $comp3{$comp}{$prots[1]} && defined $comp3{$comp}{$prots[2]}) {
#					print "Removed complex Clus3 $clus, score $clusters{$clus}{SCORE}, is a complex!!\n";
#					$comp3_removed{$clus}{E}{$prots[0]} = 1;
#					$comp3_removed{$clus}{E}{$prots[1]} = 1;
#					$comp3_removed{$clus}{E}{$prots[2]} = 1;
#					$comp3_removed{$clus}{SCORE} = $clusters{$clus}{SCORE};
#					$iscomp = 1;
#				}
#			}
#			if ($iscomp==0) {
#				print "Removed non-complex clus3 $clus, score $clusters{$clus}{SCORE}\n";
#			}
			if ($FUZZY_EDGE_REMOVE==0) {
				delete $clusters{$clus};
			}
			else {
#				my $min_remove_prob = 1;
#				if (defined $remove_edges{$edge1} && $remove_edges{$edge1} < $min_remove_prob) {
#					$min_remove_prob = $remove_edges{$edge1};
#				}
#				if (defined $remove_edges{$edge2} && $remove_edges{$edge2} < $min_remove_prob) {
#					$min_remove_prob = $remove_edges{$edge2};
#				}
#				if (defined $remove_edges{$edge3} && $remove_edges{$edge3} < $min_remove_prob) {
#					$min_remove_prob = $remove_edges{$edge3};
#				}
#				$clusters{$clus}{SCORE} *= $min_remove_prob;
				if (defined $remove_edges{$edge1}) { $clusters{$clus}{SCORE} = $clusters{$clus}{SCORE} * 2/3 + $clusters{$clus}{SCORE} * $remove_edges{$edge1} / 3; }
				if (defined $remove_edges{$edge2}) { $clusters{$clus}{SCORE} = $clusters{$clus}{SCORE} * 2/3 + $clusters{$clus}{SCORE} * $remove_edges{$edge2} / 3; }
				if (defined $remove_edges{$edge3}) { $clusters{$clus}{SCORE} = $clusters{$clus}{SCORE} * 2/3 + $clusters{$clus}{SCORE} * $remove_edges{$edge3} / 3; }
			}
			$num_clus3_removed++;
		}
	}
	else { die ("clus size ".(scalar @prots)."\n"); }
}
print "Num size-2 clusters originally = $num_clus2_orig, num size-3 clusters originally = $num_clus3_orig\n";
print "Num size-2 clusters removed = $num_clus2_removed, num size-3 clusters removed = $num_clus3_removed\n";
print "Num clusters remaining = ".(scalar keys %clusters)."\n";


# print
open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $clus (sort {$clusters{$b}{SCORE} <=> $clusters{$a}{SCORE}} keys %clusters) {
	my @prots = sort keys %{$clusters{$clus}{E}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
	print OUTPUTFILE "C$clus(".$size."_$clusters{$clus}{SCORE}): $protsstring\n";
}	

#
#
#
#print "----------------Comp-2 removed:\n";
#foreach my $comp (sort {$comp2_removed{$b}{SCORE} <=> $comp2_removed{$a}{SCORE}} keys %comp2_removed) {
#	print "$comp\t$comp2_removed{$comp}{SCORE}\n";
#	foreach my $prot (keys %{$comp2_removed{$comp}{E}}) {
#		foreach my $annot (keys %{$most_specific_annots{$prot}}) {
#			print "$prot\t$annot\t$go_scheme{$annot}{level}\t$go_scheme{$annot}{name}\n";
#		}
#	}
#}
#
#print "\n---------------- Comp-3 removed:\n";
#foreach my $comp (sort {$comp3_removed{$b}{SCORE} <=> $comp3_removed{$a}{SCORE}} keys %comp3_removed) {
#	print "$comp\t$comp3_removed{$comp}{SCORE}\n";
#	foreach my $prot (keys %{$comp3_removed{$comp}{E}}) {
#		foreach my $annot (keys %{$most_specific_annots{$prot}}) {
#			print "$prot\t$annot\t$go_scheme{$annot}{level}\t$go_scheme{$annot}{name}\n";
#		}
#	}
#}


# $clusters{$c}{SCORE} = score, $clusters{$c}{E}{$pid} = 1
sub ReadClusters ($$$) {
	my $filename = $_[0];
	my $clusters_ref = $_[1];
	
	open (CLUSTERS_FILE, $filename) || die $!;
	foreach my $line (<CLUSTERS_FILE>) {
		chomp($line);
		my @toks = split(' ', $line);
		(my $clus, my $rest) = split(/\(/, $toks[0]);
		(my $rest2, my $score) = split("_", $rest);
		$score = substr($score, 0, length($score)-2);
		$$clusters_ref{$clus}{SCORE} = $score + 0;
		shift @toks;
		foreach (@toks) {
			$$clusters_ref{$clus}{E}{$_} = 1;
		}
	}		
	close CLUSTERS_FILE;	
}


