use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use POSIX;
use Storable qw(dclone);
use Getopt::Std;

# inputs:
# -1 <SWC clusters file. Small clusters will be removed from this file>
# -2 <DECOMP clusters file. Small clusters will be removed from this file>
# -3 <SSS clusters file>
# -x <first clusters score reduction factor>
# -y <second clusters score reduction factor>
# -z <third clusters score reduction factor>
# -t <Match threshold for duplicate clusters. Small clusters are NOT checked with other small or large complexes>
# -c <Combine score mode. 0 = just take score of highest source. 1 = sum the scores (default)>
# -o <Output filename>


my %argopts;
if (! getopts('1:2:3:t:c:o:x:y:z:', \%argopts)) { die "invalid arguments"; }

my $firstfilename = "";
my $secondfilename = "";
my $thirdfilename = "";
if (defined $argopts{'1'}) {
	$firstfilename = $argopts{'1'};
}
if (defined $argopts{'2'}) {
	$secondfilename = $argopts{'2'};
}
if (defined $argopts{'3'}) {
	$thirdfilename = $argopts{'3'};
}
my $matchscore_thr = 1;
if (defined $argopts{'t'}) {
	$matchscore_thr = $argopts{'t'} + 0;
}
my $outputfilename = $argopts{'o'};
my $COMBINE_SCORE_MODE = 1;
if (defined $argopts{'c'}) {
	$COMBINE_SCORE_MODE = $argopts{'c'} + 0;
}

my $first_scorefac = 1;
my $second_scorefac = 1;
my $third_scorefac = 1;
if (defined $argopts{'x'}) {
	$first_scorefac = $argopts{'x'}+0;
}
if (defined $argopts{'y'}) {
	$second_scorefac = $argopts{'y'}+0;
}
if (defined $argopts{'z'}) {
	$third_scorefac = $argopts{'z'}+0;
}


my %firstclusters; # $clusters{$c}{SCORE} = score, $clusters{$c}{ELEMENTS}{$pid} = 1
my %secondclusters;
my %thirdclusters;
my %combinedclusters;
if (defined $argopts{'1'}) {
	ReadClusters($firstfilename, \%firstclusters, "1", $first_scorefac);
} 		
if (defined $argopts{'2'}) {
	ReadClusters($secondfilename, \%secondclusters, "2", $second_scorefac); 		
}
if (defined $argopts{'3'}) {
	ReadClusters($thirdfilename, \%thirdclusters, "3", $third_scorefac);	
}

# remove small clusters from %firstclusters and %secondclusters
foreach my $clus (keys %firstclusters) {
	if (scalar keys %{$firstclusters{$clus}{ELEMENTS}} < 4) {
		delete $firstclusters{$clus};
	}
}
foreach my $clus (keys %secondclusters) {
	if (scalar keys %{$secondclusters{$clus}{ELEMENTS}} < 4) {
		delete $secondclusters{$clus};
	}
}
print "Removed small clusters from firstclusters and secondclusters\n";

AnalyzeClusters(\%firstclusters, "1");
AnalyzeClusters(\%secondclusters, "2");
AnalyzeClusters(\%thirdclusters, "3");


if ($matchscore_thr < 1) {
	print "Matchscore threshold = $matchscore_thr < 1, so removing duplicates from individual sources first\n";
	print "Removing from first clusters...\n";
	RemoveDuplicateClusters(\%firstclusters, $matchscore_thr);
	print "Removing from second clusters...\n";
	RemoveDuplicateClusters(\%secondclusters, $matchscore_thr);
	print "Removing from third clusters...\n";
	RemoveDuplicateClusters(\%thirdclusters, $matchscore_thr);
	AnalyzeClusters(\%firstclusters, "1");
	AnalyzeClusters(\%secondclusters, "2");
	AnalyzeClusters(\%thirdclusters, "3");
}
print "\n";


my %clusters_overlap_info = (); # $clusters_overlap_info{3} = # clusters from all 3 sources, {2} = from 2 sources, etc. {1}{1} = clusters uniquely from source 1, {1}{2}, clusters uniquely from source 2, etc.
CombineClusters (\%firstclusters, \%secondclusters, \%thirdclusters, \%combinedclusters);
print "Combined clusters: ".(scalar keys %firstclusters)." first clusters, ".(scalar keys %secondclusters)." second clusters, ".(scalar keys %thirdclusters)." third clusters, ".(scalar keys %combinedclusters)." combined clusters\n";
if ($matchscore_thr < 1) {
	RemoveDuplicateClustersAnalyze(\%combinedclusters, $matchscore_thr, \%clusters_overlap_info);
	print "Removed duplicates, ".(scalar keys %combinedclusters)." combined clusters\n";
	print "Overlapping clusters info:\n";
	foreach my $src (sort keys %clusters_overlap_info) {
		if ($src==1) {
			foreach my $uniqsrc (sort keys %{$clusters_overlap_info{$src}}) {
				print "$uniqsrc only\t$clusters_overlap_info{$src}{$uniqsrc}\n";
			}
		}
		else {
			print "$src\t$clusters_overlap_info{$src}\n";
		}
	}
}
print "\n";
	
open (OUTPUTFILE, ">$outputfilename") || die $!;
foreach my $clus (sort {$combinedclusters{$b}{SCORE} <=> $combinedclusters{$a}{SCORE}} keys %combinedclusters) {
	my @prots = sort keys %{$combinedclusters{$clus}{ELEMENTS}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
	my @srcs = split(/\|/, $combinedclusters{$clus}{SRC});
	my $num_srcs = scalar @srcs;
	my $srcs_str = "|$num_srcs|".$combinedclusters{$clus}{SRC}."|";
#	if ($size < 4) { next; }
	print OUTPUTFILE "C$clus$srcs_str(".$size."_$combinedclusters{$clus}{SCORE}): $protsstring\n";
}	
close OUTPUTFILE;



# Analyze and print characteristics of combined clusters
print "-------- Combined clusters characteristics: --------\n";
AnalyzeClusters(\%combinedclusters, "COMBINED");
print "\n";


#my %clusters_overlap_info_alliters = (); # $clusters_overlap_info_alliters{5}{MEAN} = mean # clusters generated by all 5 methods across all iters, etc., {STD} = std
#for (my $iter = 0; $iter < $numiters; $iter++) {
#	foreach my $src (keys %{$clusters_overlap_info{$iter}}) {
#		$clusters_overlap_info_alliters{$src}{MEAN} += $clusters_overlap_info{$iter}{$src};
#		$clusters_overlap_info_alliters{$src}{STD} += ($clusters_overlap_info{$iter}{$src}) ** 2;
#	}
#}
#foreach my $src (keys %clusters_overlap_info_alliters) {
#	$clusters_overlap_info_alliters{$src}{MEAN} /= $numiters;
#	$clusters_overlap_info_alliters{$src}{STD} /= $numiters;
#	$clusters_overlap_info_alliters{$src}{STD} = $clusters_overlap_info_alliters{$src}{STD} - $clusters_overlap_info_alliters{$src}{MEAN} ** 2;
#	$clusters_overlap_info_alliters{$src}{STD} = $clusters_overlap_info_alliters{$src}{STD} ** .5;
#}
#	
#print "NumSrcs\tMean\tStd\n";
#foreach my $src (sort keys %clusters_overlap_info_alliters) {
#	print "$src\t$clusters_overlap_info_alliters{$src}{MEAN}\t$clusters_overlap_info_alliters{$src}{STD}\n";
#}




# ------------------------------------------------------------



# match cluster with a lower bound. If match score is definitely below lower bound, then just return 0, no need to compute exactly
sub MatchClustersLB ($$$) {
	my $clus1_ref = $_[0];
	my $clus2_ref = $_[1];
	my $lowerbound = $_[2];
	
	if (scalar keys %{$clus1_ref} < scalar keys %{$clus2_ref}) {
		my $tmp = $clus1_ref;
		$clus1_ref = $clus2_ref;
		$clus2_ref = $tmp;
	}
	
	if ((scalar keys %{$clus2_ref}) / (scalar keys %{$clus1_ref}) < $lowerbound) { 
		return 0;
	}
	
	# for speed, clus1 should be the bigger cluster
	my $num_intersect = 0;
	foreach my $elem (keys %$clus2_ref) {
		if (defined $$clus1_ref{$elem}) {
			$num_intersect++;
		}
	}
	my $num_union = (scalar keys %{$clus1_ref}) + (scalar keys %{$clus2_ref}) - $num_intersect;
	return ($num_intersect/$num_union);
}




sub AnalyzeClusters {
	my $clusters_ref = $_[0];
	my $clustersource = $_[1];
	
	if (scalar keys %{$clusters_ref}==0) { return; }
	
	my $numclus = 0;
	my $sizemean = 0;
	my $sizestd = 0;
	foreach my $clus (keys %{$clusters_ref}) {
		$numclus++;
		$sizemean += scalar keys %{$$clusters_ref{$clus}{ELEMENTS}};
		$sizestd += (scalar keys %{$$clusters_ref{$clus}{ELEMENTS}})**2;				
	}
	if ($numclus==0) { return; }
	$sizemean /= $numclus;
	$sizestd /= $numclus;
	$sizestd = $sizestd - $sizemean**2;
	$sizestd = $sizestd ** .5;
	print "$clustersource: num clusters = $numclus, size mean = $sizemean, size std = $sizestd\n";
}

 
# $clusters{$c}{SCORE} = score, $clusters{$c}{SRC} = "1" or "2" etc, $clusters{$c}{ELEMENTS}{$pid} = 1
sub ReadClusters ($$$) {
	my $filename = $_[0];
	my $clusters_ref = $_[1];
	my $clustersource = $_[2];
	my $scorefac = $_[3];
	
	if ($filename eq "") {
		next;
	}
	open (CLUSTERS_FILE, "$filename") || die "$filename not found";
	foreach my $line (<CLUSTERS_FILE>) {
		chomp($line);
		my @toks = split(' ', $line);
		(my $clus, my $rest) = split(/\(/, $toks[0]);
		(my $rest2, my $score) = split("_", $rest);
		$score = substr($score, 0, length($score)-2);
		if ($score < .0001) { next; }
		$clus = ($clustersource."_".$clus);
		$$clusters_ref{$clus}{SCORE} = $score + 0;
		$$clusters_ref{$clus}{SCORE} *= $scorefac;
		$$clusters_ref{$clus}{SRC} = $clustersource;
		shift @toks;
		foreach (@toks) {
			$$clusters_ref{$clus}{ELEMENTS}{$_} = 1;
		}
		close CLUSTERS_FILE;
	}
}


sub RemoveDuplicateClusters ($$) {
	my $clusters_ref = $_[0];
	my $match_thres = $_[1];
	
	my @clusters = sort {$$clusters_ref{$b}{SCORE} <=> $$clusters_ref{$a}{SCORE}} keys %$clusters_ref;
	my %clusters_to_delete = ();
	
	for (my $idx1 = 0; $idx1 < scalar @clusters; $idx1++) {
		for (my $idx2 = $idx1+1; $idx2 < scalar @clusters; $idx2++) {
			my $clus1 = $clusters[$idx1];
			my $clus2 = $clusters[$idx2];
			if (defined $clusters_to_delete{$clus1} && $clusters_to_delete{$clus1}==1) { next; }
			if (defined $clusters_to_delete{$clus2} && $clusters_to_delete{$clus2}==1) { next; }
			if (scalar keys %{$$clusters_ref{$clus1}{ELEMENTS}} <= 3 || scalar keys %{$$clusters_ref{$clus2}{ELEMENTS}} <= 3) { next; }
			my $matchscore = MatchClustersLB($$clusters_ref{$clus1}{ELEMENTS}, $$clusters_ref{$clus2}{ELEMENTS}, $match_thres);
			if ($matchscore >= $match_thres) {
				$clusters_to_delete{$clus2} = 1;
#				$$clusters_ref{$clus1}{SCORE} += $$clusters_ref{$clus2}{SCORE};
			}
		}
	}
		
	print "In RemoveDuplicateClusters, num to delete = ".(scalar keys %clusters_to_delete)."\n";
	
	foreach my $clus (keys %clusters_to_delete) {
		delete $$clusters_ref{$clus};
	}	
}



sub RemoveDuplicateClustersAnalyze ($$) {
	my $clusters_ref = $_[0];
	my $match_thres = $_[1];
	my $clusters_overlap_info_ref = $_[2];
	
	my @clusters = sort {$$clusters_ref{$b}{SCORE} <=> $$clusters_ref{$a}{SCORE}} keys %$clusters_ref;
	my %clusters_to_delete = (); # $clusters_to_delete{$clus1}{$clus2} = 1, if $clus1 is deleted because it matches with $clus2 and $clus2 is higher score
	my %clusters_representatives = (); # $clusters_representatives{$clus2}{$clus1} = 1, if $clus2 is a cluster representative and $clus1 is deleted because of $clus2
	
	for (my $idx1 = 0; $idx1 < scalar @clusters; $idx1++) {
		for (my $idx2 = $idx1+1; $idx2 < scalar @clusters; $idx2++) {
			my $clus1 = $clusters[$idx1];
			my $clus2 = $clusters[$idx2];
			if (defined $clusters_to_delete{$clus1} && scalar keys %{$clusters_to_delete{$clus1}}>0) { next; }
			if (defined $clusters_to_delete{$clus2} && scalar keys %{$clusters_to_delete{$clus2}}>0) { next; }
			if (scalar keys %{$$clusters_ref{$clus1}{ELEMENTS}} <= 3 || scalar keys %{$$clusters_ref{$clus2}{ELEMENTS}} <= 3) { next; }
			my $matchscore = MatchClustersLB($$clusters_ref{$clus1}{ELEMENTS}, $$clusters_ref{$clus2}{ELEMENTS}, $match_thres);
			if ($matchscore >= $match_thres) {
				$clusters_to_delete{$clus2}{$clus1} = 1;
				$clusters_representatives{$clus1}{$clus2} = 1;
				if ($COMBINE_SCORE_MODE==0) {
					
				}
				elsif ($COMBINE_SCORE_MODE==1) {
					my $contribute_score = 1;
					foreach my $clus (keys %{$clusters_representatives{$clus1}}) {
						if ($clus eq $clus2) { next; }
						if ($$clusters_ref{$clus}{SRC} eq $$clusters_ref{$clus2}{SRC}) {
							$contribute_score = 0;
							last;
						}
					}
					if ($contribute_score==1) {
						$$clusters_ref{$clus1}{SCORE} += $$clusters_ref{$clus2}{SCORE};
					}	
				}
			}
		}
	}
	
	# add clusters that are not deleted, and do not represent other clusters, to %clusters_representatives
	foreach my $clus (@clusters) {
		if (!defined $clusters_to_delete{$clus} && !defined $clusters_representatives{$clus}) {
			$clusters_representatives{$clus} = ();
		}
	}
		
	# do some checks
	my %clusters_check = ();
	foreach my $clus (keys %clusters_representatives) {
		$clusters_check{$clus} = 1;
		foreach my $clus_del (keys %{$clusters_representatives{$clus}}) {
			if (!defined $clusters_to_delete{$clus_del}) { die; }
			if (!defined $clusters_to_delete{$clus_del}{$clus}) { die; }
			$clusters_check{$clus_del} = 1;
		}
	}
	if (scalar keys %clusters_check != scalar @clusters) { die; }
	
	# set the clusters' sources and print some info	
	print "In RemoveDuplicateClustersAnalyze, num to delete = ".(scalar keys %clusters_to_delete).", num representatives = ".(scalar keys %clusters_representatives)."\n";
#	print "Cluster representatives:\n";
#	foreach my $clus_rep (sort {scalar keys %{$clusters_representatives{$b}} <=> scalar keys %{$clusters_representatives{$a}}} keys %clusters_representatives) {
#		print "$clus_rep\t".(scalar keys %{$clusters_representatives{$clus_rep}});
#		foreach my $clus_del (keys %{$clusters_representatives{$clus_rep}}) {
#			print "\t$clus_del";
#		}
#		print "\n";
#	}
	# $clusters_overlap_info{3} = # clusters generated by all 3 methods, {1}{1} = num clusters uniquely from 1, etc, {1|2} = num overlap between 1 and 2
	foreach my $clus_rep (keys %clusters_representatives) {
		if (scalar keys %{$clusters_representatives{$clus_rep}} == 0) {
			my $src = $$clusters_ref{$clus_rep}{SRC};
			$$clusters_overlap_info_ref{1}{$src}++;
		}
		else {
			my %clus_srcs = ();
			foreach my $clus (keys %{$clusters_representatives{$clus_rep}}) {
				$clus_srcs{$$clusters_ref{$clus}{SRC}} = 1;
			}
			$clus_srcs{$$clusters_ref{$clus_rep}{SRC}} = 1;
			my $new_srcs = join('|', sort keys %clus_srcs);
			$$clusters_ref{$clus_rep}{SRC} = $new_srcs;
			my $num_srcs = scalar keys %clus_srcs;
			$$clusters_overlap_info_ref{$num_srcs}++;
		}
	}
	
	
	foreach my $clus (keys %clusters_to_delete) {
		delete $$clusters_ref{$clus};
	}	
}




sub CombineClusters($$$$$) {
	my $firstclusters_ref = $_[0];
	my $secondclusters_ref = $_[1];
	my $thirdclusters_ref = $_[2];
	my $combinedclusters_ref = $_[3];
#	print "In CombinemClusters, num cmc = ".(scalar keys %{$cmcclusters_ref}).", num ipca = ".(scalar keys %{$ipcaclusters_ref})."\n";
	foreach my $clus (keys %{$firstclusters_ref}) {
		$$combinedclusters_ref{"$clus"}{SCORE} = $$firstclusters_ref{$clus}{SCORE};
		$$combinedclusters_ref{"$clus"}{SRC} = $$firstclusters_ref{$clus}{SRC};
		foreach my $prot (keys %{$$firstclusters_ref{$clus}{ELEMENTS}}) {
			$$combinedclusters_ref{"$clus"}{ELEMENTS}{$prot} = 1;;
		}
	}
	foreach my $clus (keys %{$secondclusters_ref}) {
		$$combinedclusters_ref{"$clus"}{SCORE} = $$secondclusters_ref{$clus}{SCORE};
		$$combinedclusters_ref{"$clus"}{SRC} = $$secondclusters_ref{$clus}{SRC};
		foreach my $prot (keys %{$$secondclusters_ref{$clus}{ELEMENTS}}) {
			$$combinedclusters_ref{"$clus"}{ELEMENTS}{$prot} = 1;;
		}
	}
	foreach my $clus (keys %{$thirdclusters_ref}) {
		$$combinedclusters_ref{"$clus"}{SCORE} = $$thirdclusters_ref{$clus}{SCORE};
		$$combinedclusters_ref{"$clus"}{SRC} = $$thirdclusters_ref{$clus}{SRC};
		foreach my $prot (keys %{$$thirdclusters_ref{$clus}{ELEMENTS}}) {
			$$combinedclusters_ref{"$clus"}{ELEMENTS}{$prot} = 1;;
		}
	}
#	print "num combined = ".(scalar keys %{$combinedclusters_ref})."\n";
}


#sub FindIntersection {
#	my $clusters1ref = $_[0];
#	my $clusters2ref = $_[1];
#	my $matchscore_thr = $_[2];
#	my $intersectionref = $_[3];
#	
#	my $num_score_exceeded = 0;
#	foreach my $clus1 (keys %{$clusters1ref}) {
#		foreach my $clus2 (keys %{$clusters2ref}) {
#			my $score = MatchClusters($$clusters1ref{$clus1}{ELEMENTS}, $$clusters2ref{$clus2}{ELEMENTS});
#			
#			if ($score >= $matchscore_thr) {
#				$$intersectionref{INTS}{$clus1}{SCORE} = $$clusters1ref{$clus1}{SCORE};
#				foreach my $prot (keys %{$$clusters1ref{$clus1}{ELEMENTS}}) {
#					$$intersectionref{INTS}{$clus1}{ELEMENTS}{$prot} = 1;
#				}
#				$$intersectionref{INTS}{$clus2}{SCORE} = $$clusters2ref{$clus2}{SCORE};
#				foreach my $prot (keys %{$$clusters2ref{$clus2}{ELEMENTS}}) {
#					$$intersectionref{INTS}{$clus2}{ELEMENTS}{$prot} = 1;
#				}
#			}
#		}
#	}
#	$$intersectionref{JAC} = (scalar keys %{$$intersectionref{INTS}}) / ((scalar keys %{$clusters1ref}) + (scalar keys %{$clusters2ref}));
#}


