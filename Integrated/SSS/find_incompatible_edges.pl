#!/usr/bin/perl -w
use strict;
use IO::Handle;
use Getopt::Std;


# When decomposing, if a protein's most specific annotation is ancestor of a decomposition term, then also include that protein in the subnet
# inputs:
# -i <input PPI file>
# -a <Input GO annotations file>
# -s <GO scheme file>
# -d <GO decomposition terms>
# -o <output basefile>
# perl decomp_ppi_goterms.pl -i "data_ppibiogrid.txt" -a "my_go_associations.sgd_propagated.txt" -s "go_scheme_all_mine.txt" -d "decompGOterms30.txt" -o "data_ppibiogrid_godecomp30"

my %argopts;
if (! getopts('i:a:d:s:t:o:', \%argopts)) { die "invalid arguments"; }

# threshold to add a prot with non-specific annotation to a subnet. 0 = add if the prot has any GO term with a subnet-GO descendent, 1 = add only if the prot has the subnet's goterm (same as original decomp)
# Leave it at 1 to output all removed edges with their removal probabilities, so that we can use the probabilities to reduce the scores later on
my $SUBNET_THRESHOLD = 1;


# read decomposition terms
open (DECOMPFILE, "$argopts{'d'}") || die $!;
my %decomp_goterms = ();
foreach my $line (<DECOMPFILE>) {
	chomp $line;
	$decomp_goterms{$line} = 1;
}
	
	
# read GO scheme, only for CC
my %go_scheme = (); # $go_scheme{$goid}{"level"} = level, {"name"} = name, {"parents"}, {"children"}
open (GOSCHEMEFILE, "$argopts{'s'}") || die $!;
foreach my $line (<GOSCHEMEFILE>) {
	chomp $line;
	my @tokens = split(' ', $line, 2);
	my @tokens2 = split(/\|/, $tokens[1]);
	my $goclass = $tokens2[0];
	if ($goclass ne "cellular_component") { next; }
	$go_scheme{$tokens[0]}{"level"} = ($tokens2[1]+0);
	$go_scheme{$tokens[0]}{"name"} = $tokens2[2];
	if (scalar @tokens2 == 4) {
		my @tokens3 = split(',', $tokens2[3]);
		foreach (@tokens3) {
			$go_scheme{$tokens[0]}{"parents"}{$_} = 1;
			$go_scheme{$_}{"children"}{$tokens[0]} = 1;
		}
	}
}
print "read go scheme, num CC GO terms = ".(scalar keys %go_scheme)."\n";
	
	
# read GO annotations
my %annots = (); # $annots{$prot}{$goid} = 1. 
open(ANNOTFILE, "$argopts{'a'}") || die $!;
foreach my $line (<ANNOTFILE>) { 
	chomp($line);
  (my $prot, my $goid) = split(/\|/, $line);
	if (!defined $go_scheme{$goid}) { next; }
 	$annots{$prot}{$goid} = 1;
}
print "annotations read for ".(scalar keys %annots)." prots\n";


my %goterms_count = (); # $goterms_count{$goid} = num prots annotated to it
foreach my $goid (keys %go_scheme) {
	$goterms_count{$goid} = 0;
}
foreach my $prot (keys %annots) {
	foreach my $goid (keys %{$annots{$prot}}) {
		$goterms_count{$goid}++;
	}
}
#foreach my $goterm (keys %goterms_count) {
#	print "$goterm\t$goterms_count{$goterm}\t$go_scheme{$goterm}{name}\n";
#}



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
	
foreach my $prot (keys %ppi_prots) {
	if (!defined $annots{$prot}) {
		$annots{$prot}{"0005575"} = 1;
#		print "PPI protein $prot has no CC annotations\n";
	}
}



# get the set of most specific annots (those terms with no descendents) for each prot
my %most_specific_annots = ();
foreach my $prot (keys %annots) {
	foreach my $goid (keys %{$annots{$prot}}) {
		my $child_annotated = 0;
		foreach my $child (keys %{$go_scheme{$goid}{children}}) {
			if (defined $annots{$prot}{$child}) {
				$child_annotated = 1;
				last;
			}
		}
		if ($child_annotated==0) {
			$most_specific_annots{$prot}{$goid} = 1;
		}
	}
	if (!defined $most_specific_annots{$prot}) { die; }
}


# for each decomp term g, for each of its ancestors a, calculate n_g / n_a
# $decomp_terms_ancestors_probs{$goterm}{$ancestor} = num annot to $goterm / num annot to $ancestor
my %decomp_terms_ancestors_probs = (); 
foreach my $goterm (keys %decomp_goterms) {
	my %ancestors = ();
	GetAncestors($goterm, \%ancestors);
	foreach my $anc (keys %ancestors) {
		$decomp_terms_ancestors_probs{$goterm}{$anc} = $goterms_count{$goterm} / $goterms_count{$anc};
	}
}
print "Calculated decomp terms ancestors probs, for ".(scalar keys %decomp_terms_ancestors_probs)." decomp terms\n";
foreach my $goterm (keys %decomp_terms_ancestors_probs) {
	print "$goterm\t$go_scheme{$goterm}{name}\n";
	foreach my $anc (sort {$decomp_terms_ancestors_probs{$goterm}{$b} <=> $decomp_terms_ancestors_probs{$goterm}{$a}} keys %{$decomp_terms_ancestors_probs{$goterm}}) {
		print "$anc\t$go_scheme{$anc}{name}\t$decomp_terms_ancestors_probs{$goterm}{$anc}\n";
	}
}



my %edges_discarded = ();
foreach my $edge (keys %ppi_edges) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	my $remove = 1;
	my $highest_prob = 0;
	foreach my $goid (keys %decomp_goterms) {
		my $prot1insubnet_prob = PutInSubnetFuzzy($prot1, $goid);
		my $prot2insubnet_prob = PutInSubnetFuzzy($prot2, $goid);
		my $this_goid_lowest_prob = $prot1insubnet_prob < $prot2insubnet_prob ? $prot1insubnet_prob : $prot2insubnet_prob;
		if ($this_goid_lowest_prob > $highest_prob) {
			$highest_prob = $this_goid_lowest_prob;
		}
		if ($prot1insubnet_prob >= $SUBNET_THRESHOLD && $prot2insubnet_prob >= $SUBNET_THRESHOLD) { 
			$remove = 0;
			last;
		}
	}
	if ($remove==1) {
		$edges_discarded{$edge} = $highest_prob;
	}
}


open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $edge (keys %edges_discarded) {
	print OUTPUTFILE "$edge\t$edges_discarded{$edge}\n";
}


##### debug #####
#
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
#close COMPLEXFILE;
#print "num complexes read = ".(scalar keys %complexes)."\n";
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
#
#foreach my $edge (keys %edges_discarded) {
#	(my $prot1, my $prot2) = split(/\|/, $edge);
#	foreach my $comp (keys %comp2) {
#		if (defined $comp2{$comp}{$prot1} && defined $comp2{$comp}{$prot2}) {
#			print "$edge discarded, from size-2 complex $comp\n";
#			foreach my $goterm (keys %decomp_goterms) {
#				if (defined $annots{$prot1}{$goterm}) {
#					print "$prot1 has decomp term $goterm\n";
#				}
#				else {
#					my $highest_prob = -1;
#					my $highest_prob_annot = "";
#					foreach my $spec_annot (keys %{$most_specific_annots{$prot1}}) {
#						if (!defined $decomp_terms_ancestors_probs{$goterm}{$spec_annot}) { next; }
#						if ($decomp_terms_ancestors_probs{$goterm}{$spec_annot} > $highest_prob) {
#							$highest_prob = $decomp_terms_ancestors_probs{$goterm}{$spec_annot};
#							$highest_prob_annot = $spec_annot;
#						}
#					}
#					if ($highest_prob>-1) { print "$prot1 has annot $highest_prob_annot which is ancestor of $goterm, prob = $highest_prob\n"; }
#						
#				}
#				if (defined $annots{$prot2}{$goterm}) {
#					print "$prot1 has decomp term $goterm\n";
#				}
#				else {
#					my $highest_prob = -1;
#					my $highest_prob_annot = "";
#					foreach my $spec_annot (keys %{$most_specific_annots{$prot2}}) {
#						if (!defined $decomp_terms_ancestors_probs{$goterm}{$spec_annot}) { next; }
#						if ($decomp_terms_ancestors_probs{$goterm}{$spec_annot} > $highest_prob) {
#							$highest_prob = $decomp_terms_ancestors_probs{$goterm}{$spec_annot};
#							$highest_prob_annot = $spec_annot;
#						}
#					}
#					if ($highest_prob>-1) { print "$prot2 has annot $highest_prob_annot which is ancestor of $goterm, prob = $highest_prob\n";}
#				}
#			}
#		}
#	}
#	foreach my $comp (keys %comp3) {
#		if (defined $comp3{$comp}{$prot1} && defined $comp3{$comp}{$prot2}) {
#			print "$edge discarded, from size-3 complex $comp\n";
#			foreach my $goterm (keys %decomp_goterms) {
#				if (defined $annots{$prot1}{$goterm}) {
#					print "$prot1 has decomp term $goterm\n";
#				}
#				if (defined $annots{$prot2}{$goterm}) {
#					print "$prot1 has decomp term $goterm\n";
#				}
#			}
#		}
#	}
#}
#			



###############################################################################################

#sub PutInSubnet {
#	my $prot = $_[0];
#	my $goterm = $_[1];
#		
#	# protein is annotated to $goterm
#	if (defined $annots{$prot}{$goterm}) {
#		return 1;
#	}
#	
#	# protein has no annotations at all
#	if (!defined $most_specific_annots{$prot}) { 
#		die "Prot $prot has no most_specifc_annots\n";
#	}
#	
#	# at least one of protein's most specific annotation is ancestor of $goterm
#	foreach my $spec_annot (keys %{$most_specific_annots{$prot}}) {
#		if (IsAncestor($spec_annot, $goterm)==1) {
#	#		print "Prot $prot most specific annot is $most_specific_annots{$prot}, which is ancestor of $goterm, so put it in $goterm subnet\n";
#			return 1;
#		}
#	}
#	return 0;	
#}



sub PutInSubnetFuzzy {
	my $prot = $_[0];
	my $goterm = $_[1];
		
	# protein is annotated to $goterm
	if (defined $annots{$prot}{$goterm}) {
		return 1;
	}
	
	# protein has no annotations at all
	if (!defined $most_specific_annots{$prot}) { 
		die "Prot $prot has no most_specifc_annots\n";
	}
	
	# $goterm should be a decomposition term
	if (!defined $decomp_terms_ancestors_probs{$goterm}) { die; }

	# find highest probability that $prot is annotated to $goterm, based on whther $prot is annotated to $goterm's ancestors
	my $highest_prob = 0;
	my $highest_prob_annot = "";
	foreach my $spec_annot (keys %{$most_specific_annots{$prot}}) {
		if (!defined $decomp_terms_ancestors_probs{$goterm}{$spec_annot}) { next; }
		if ($decomp_terms_ancestors_probs{$goterm}{$spec_annot} > $highest_prob) {
			$highest_prob = $decomp_terms_ancestors_probs{$goterm}{$spec_annot};
			$highest_prob_annot = $spec_annot;
		}
	}
	return $highest_prob;
	
}


# quirk: $goterm is considered its own ancestor
sub IsAncestor {
	my $ancestor = $_[0];
	my $goterm = $_[1];
	
	if ($goterm eq $ancestor) {
		return 1;
	}
	
	foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
		if (IsAncestor($ancestor, $parent)==1) {
			return 1;
		}
	}
	return 0;	
}

# quirk: $goterm is considered its own ancestor
sub GetAncestors {
	my $goterm = $_[0];
	my $ancestors_ref = $_[1];
	
	$$ancestors_ref{$goterm} = 1;
	foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
		GetAncestors($parent, $ancestors_ref);
	}
}