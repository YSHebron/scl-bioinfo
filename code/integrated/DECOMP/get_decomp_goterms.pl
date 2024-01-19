#!/usr/bin/perl -w
use strict;
use IO::Handle;
use Getopt::Std;



# inputs:
# -a <Input GO annotations file>
# -s <GO scheme file>
# -n <ngo parameter>
# -m <how to fill Missing GO terms. 0 = do nothing (default). 1 = just add the most general (highest-level) missing terms>
# -t <minimum number of proteins annotated to each selected term. Default 100>
# -o <outputfile>
# perl get_decomp_goterms.pl -a "my_go_associations.sgd_propagated.txt" -s "go_scheme_all_mine.txt" -n 300 -o "decompGOterms300_2.txt"

my %argopts;
if (! getopts('a:s:n:m:t:o:', \%argopts)) { die "invalid arguments"; }

my $ngo = $argopts{'n'}+0;

my $FILL_MISSING_MODE = 0;
if (defined $argopts{'m'}) {
	$FILL_MISSING_MODE = $argopts{'m'}+0;
}

my $TERM_MINPROT = 100;
if (defined $argopts{'t'}) {
	$TERM_MINPROT = $argopts{'t'}+0;
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
	
	
# read GO annotations, only for CC
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
#print "Num prots annotated to GO terms:\n";
#foreach my $goid (sort keys %goterms_count) {
#	print "$goid\t$goterms_count{$goid}\n";
#}
#print "\n";

my %decomp_goterms = ();
CheckGOTerm("0005575"); # start check on cellular component

# Check if any of the chosen terms have ancestors that are also chosen
foreach my $goterm (keys %decomp_goterms) {
	foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
		if (HasSelectedAncestor($parent, \%decomp_goterms)==1) {
			print "Selected term $goterm has selected ancestor! Deleting it\n";
			delete $decomp_goterms{$goterm};
		}
	}
}

print "Num decomposition terms = ".(scalar keys %decomp_goterms)."\n";
foreach my $goid (sort  {$goterms_count{$b} <=> $goterms_count{$a}} keys %decomp_goterms) {
	print "$goid\t$goterms_count{$goid}\n";
}
print "\n";



# GO terms that are more specific than ngo (annotated to fewer than ngo), yet have no ancestors that are informative
my %missedGOTerms = ();
CheckMissedGOTerms("0005575"); # start check on cellular component
print "Num missed GO terms (excluding their descendents) = ".(scalar keys %missedGOTerms)."\n";
#foreach my $goterm (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys %missedGOTerms) {
#	print "$goterm\t$goterms_count{$goterm}\n";
#	foreach my $parent (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys  %{$go_scheme{$goterm}{"parents"}}) {
#		print "\t$parent\t$goterms_count{$parent}\n";
#	}
#}



############ Try to fill in the missed terms #################


# just add the missing terms (those with at least TERM_MINPROT prots)
if ($FILL_MISSING_MODE==1) {
	while (1) {
		my @missed_terms_sorted = sort {$goterms_count{$b} <=> $goterms_count{$a}} keys %missedGOTerms;
		if (scalar @missed_terms_sorted == 0 || $goterms_count{$missed_terms_sorted[0]} < $TERM_MINPROT) {
			last;
		}
		my $add_candidate = $missed_terms_sorted[0];
		$decomp_goterms{$add_candidate} = 1;
#		print "Adding $add_candidate\n";
			
		# Check if any of the chosen terms have ancestors/descendents that are also chosen
		foreach my $goterm (keys %decomp_goterms) {
			foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
				if (HasSelectedAncestor($parent, \%decomp_goterms)==1) {
#					print "Selected term $goterm has selected ancestor! Deleting it\n";
					delete $decomp_goterms{$goterm};
				}
			}
		}
		
		print "Num decomposition terms = ".(scalar keys %decomp_goterms)."\n";
#		foreach my $goid (sort  {$goterms_count{$b} <=> $goterms_count{$a}} keys %decomp_goterms) {
#			print "$goid\t$goterms_count{$goid}\n";
#		}
#		print "\n";
					
		%missedGOTerms = ();
		CheckMissedGOTerms("0005575"); # start check on cellular component
#		print "Num missed GO terms (excluding their descendents) = ".(scalar keys %missedGOTerms)."\n";
#		foreach my $goterm (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys %missedGOTerms) {
#			print "$goterm\t$goterms_count{$goterm}\n";
#			foreach my $parent (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys  %{$go_scheme{$goterm}{"parents"}}) {
#				print "\t$parent\t$goterms_count{$parent}\n";
#			}
#		}
			
	} # end while(1)
	
} # end FILL_MISSING_MODE==1





	
# check if any selected term has fewer than TERM_MINPROT proteins
foreach my $goid (sort keys %decomp_goterms) {
	if ($goterms_count{$goid} < $TERM_MINPROT) {
		print "**** WARNING! Term $goid has $goterms_count{$goid} proteins < TERM_MINPROT $TERM_MINPROT, so deleting it ****\n";
		delete $decomp_goterms{$goid};
	}
}


print "\n\n================ Final selected terms =================\n";
print "Num decomposition terms = ".(scalar keys %decomp_goterms)."\n";
foreach my $goid (sort  {$goterms_count{$b} <=> $goterms_count{$a}} keys %decomp_goterms) {
	print "$goid\t$goterms_count{$goid}\t$go_scheme{$goid}{name}\n";
}
print "\n";

%missedGOTerms = ();
CheckMissedGOTerms("0005575"); # start check on cellular component
print "Num missed GO terms (excluding their descendents) = ".(scalar keys %missedGOTerms)."\n";
#foreach my $goterm (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys %missedGOTerms) {
#	print "$goterm\t$goterms_count{$goterm}\n";
#	foreach my $parent (sort {$goterms_count{$b} <=> $goterms_count{$a}} keys  %{$go_scheme{$goterm}{"parents"}}) {
#		print "\t$parent\t$goterms_count{$parent}\n";
#	}
#}


open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $goid (sort keys %decomp_goterms) {
	print OUTPUTFILE "$goid\n";
}

###################


# Find GO terms g, s.t. g has >= ngo prots, but none of its children has >= ngo prots
sub CheckGOTerm {
	my $goterm = $_[0];
	if (!defined $go_scheme{$goterm}) { die "$goterm"; }
	if ($goterms_count{$goterm} < $ngo) {
		return;
	}
	if (!defined $go_scheme{$goterm}{"children"} || scalar keys %{$go_scheme{$goterm}{"children"}}==0) {
		$decomp_goterms{$goterm} = 1;
		return;
	}
	my $addthisterm = 1;
	foreach my $child (keys %{$go_scheme{$goterm}{"children"}}) {
		if ($goterms_count{$child} >= $ngo) {
			$addthisterm = 0;
			CheckGOTerm($child);
		}
	}
	if ($addthisterm==1) {
		$decomp_goterms{$goterm} = 1;
		return;
	}
	
}



# Find GO terms g, s.t. g annotated to >= ngo prots, and either (at least one child annotated to less than ngo prots, OR g has no children)
sub CheckGOTermComprehensive {
	my $goterm = $_[0];
	if (!defined $go_scheme{$goterm}) { die "$goterm"; }
	if ($goterms_count{$goterm} < $ngo) {
		return;
	}
	if (!defined $go_scheme{$goterm}{"children"} || scalar keys %{$go_scheme{$goterm}{"children"}}==0) {
		$decomp_goterms{$goterm} = 1;
		return;
	}
	foreach my $child (keys %{$go_scheme{$goterm}{"children"}}) {
		if ($goterms_count{$child} < $ngo) {
			$decomp_goterms{$goterm} = 1;			
		}
		else {
			CheckGOTermComprehensive($child);
		}
	}	
	
}




# Find GO terms G, s.t. g in G annotated to <= ngo prots, and all its parents have > ngo prots
sub CheckGOTermComprehensive2 {
	my $goterm = $_[0];
	if (!defined $go_scheme{$goterm}) { die "$goterm"; }
	if (defined $decomp_goterms{$goterm}) { return;}
	if ($goterms_count{$goterm} > $ngo) {
		foreach my $child (keys %{$go_scheme{$goterm}{"children"}}) {
			CheckGOTermComprehensive2($child);
		}
		return;
	}
	if ($goterms_count{$goterm} <= $ngo) {
		my $add_current = 1;
		foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
			if ($goterms_count{$parent} <= $ngo) {
				$add_current = 0;
			}
		}
		if ($add_current==1) {
			$decomp_goterms{$goterm} = 1;
		}
	}
	
}



# search for GO terms that are more specific than ngo (annotated to fewer than ngo), yet have no ancestors that are selected. Ignore those with 0 prots
# Put in %missedGOTerms
sub CheckMissedGOTerms {
	my $goterm = $_[0];
	if (!defined $go_scheme{$goterm}) { die "$goterm"; }
	if ($goterms_count{$goterm}==0) { return; }
	if ($goterms_count{$goterm} <= $ngo) {
		my $all_parents_lt_ngo = 1;
		my $has_inf_parent = 0;
		foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
			if ($goterms_count{$parent} >= $ngo) {
				$all_parents_lt_ngo = 0;
			}
		}
		if ($all_parents_lt_ngo==1) {
			die "GO term $goterm, num prots $goterms_count{$goterm} < $ngo, but all parents < $ngo\n";
		}
		if (HasSelectedAncestor($goterm, \%decomp_goterms)==0) {
			$missedGOTerms{$goterm} = 1;
		}
		return;
	}
	foreach my $child (keys %{$go_scheme{$goterm}{"children"}}) {
		CheckMissedGOTerms($child);
	}
	
}


sub HasSelectedAncestor {
	my $goterm = $_[0];
	my $selected_ref = $_[1];
	
	if (defined $$selected_ref{$goterm}) {
		return 1;
	}
	
	foreach my $parent (keys %{$go_scheme{$goterm}{"parents"}}) {
		if (HasSelectedAncestor($parent, $selected_ref)) {
			return 1;
		}
	}
	return 0;
	
}