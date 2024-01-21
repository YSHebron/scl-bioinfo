# Imports and Compiler Flags
use strict;					# safety pragma
use IO::Handle;				# supply object methods for I/O handles
use List::Util 'shuffle';	# imports shuffle, randomize order of elements in a list
use List::Util 'max';		# imports max, returns maximum value in a list
use POSIX;					# provides access to OS-specific functions
use Storable qw(dclone);	# deep clone function, extract from string to list with delimiter
use Getopt::Std;			# import getopts, easier parsing of cmd options and arguments


# Inputs
# -i <input clusters file, OR basefilename (see -x)>
# -c <complexes file>
# -n <Number of iterations. If defined, will iterate over the input file name, eg if input is "clusters_cmc swc5k", will process "clusters_cmc swc5k iter0.txt", "clusters_cmc swc5k iter1.txt", etc.>
# -x <Xval file, only when -n option is used>
# -m <small complexes matching requirement is always 1? default 0>
# -l <which complexes to consider for test: 0 = all, 1 = only size 4 and above, 23 = sizes 2 and 3. Note that clusters that match unconsidered complexes will NOT be false positives> default 1
# -k <remove small clusters? 0 = no (default). If this is 1, but -l option still considers small complexes, then recall will be lowered because small complexes cannot be matched. But precision may be higher because small clusters are noisy>
## Usage:
#   To keep all clusters, eval on all complexes: -l 0. Small complexes match requirement 1: -l 0 -m 1
# 	To keep only large clusters, eval on large complexes: -l 1
# 	To keep only small clusters, eval on small complexes: -l 23. Small complexes match requirement 1: -l 23 -m 1
# 	To keep only large clusters, eval on ALL complexes: -k 1 -l 0. Small complexes match requirement 1 (so that the large clusters cannot match small comp): -k 1 -l 0 -m 1
# -s <score>. only consider clusters with score over this, default 0
# -a <number of clusters to keep per iteration. Default 9999999999>
# -u <filter clusters to keep only unique ones (matchthres = 0.5): 1 or 0 (default)>

# declare argopts dictionary
my %argopts;
if (! getopts('i:c:n:x:l:k:s:u:m:a:', \%argopts)) { die "invalid arguments"; }

my $SMALLCOMP_MATCH1 = 0;
if (defined $argopts{'m'}) {
	$SMALLCOMP_MATCH1 = $argopts{'m'}+0;
}
my $REMOVE_SMALLCLUS = 0;
if (defined $argopts{'k'}) {
	$REMOVE_SMALLCLUS = $argopts{'k'}+0;
}

my $which_comps_consider = 1; 
if (defined $argopts{'l'}) {
	$which_comps_consider = $argopts{'l'};
}
my $score_threshold = 0;
if (defined $argopts{'s'}) {
	$score_threshold = $argopts{'s'} + 0;
}
my $filter_unique = 0;
if (defined $argopts{'u'}) {
	$filter_unique = $argopts{'u'} + 0;
}

my $num_iters = 0;
if (defined $argopts{'n'}) {
	$num_iters = $argopts{'n'} + 0;
	if ($num_iters == 0) { die; }
}

my $NUM_CLUS_KEEP = 999999999;
if (defined $argopts{'a'}) {
	$NUM_CLUS_KEEP = $argopts{'a'}+0;
}

my $xvalfilename;
if ($num_iters > 0) {
	$xvalfilename = $argopts{'x'};
}


# read complexes
# $complexes{$cid}{$pid} = 1 if protein $pid is in complex $cid
my %complexes = ();		# create hash %complexes, equivalent to C union T
open (COMPLEXFILE, $argopts{'c'}) || die $!;	# open $argopts{'c'} (ref cplxs) as file handle COMPLEXFILE
# for each line in COMPLEXFILE, chomp newline at the end,
# then split assign line content to $id and $comp (prot ID and cplx ID) at tab delimiter
# then convert all IDs to uppercase
# then set $complexes{$cid}{$pid} = 1 if protein $pid is in complex $cid
foreach my $line (<COMPLEXFILE>) {
	chomp($line);
	(my $id, my $comp) = split(/\t/, $line);
	$id = uc($id);
	$complexes{$comp}{$id} = 1;
}
close COMPLEXFILE;
print "num complexes read = ".(scalar keys %complexes)."\n";

# The hash %complexes would look something like this:
# {'1' => {'YLR342W' => 1,'YPR165W' => 1},
#  '2' => {'YGR032W' => 1,'YPR054W' => 1}}
# where 1 is a complex ID, and YLR342W and YPR165W are its co-complexed protein IDs
# where 2 is a complex ID, and YGR032W and YPR054W are its co-complexed protein IDs

# ==================== no iters ===============================
if ($num_iters==0) {

}

# ==================== iters ===============================
else {
		
	# $test_complexes{$iter}{$complexid} = 1, if $complexid is a TEST complex in $iter
	# We need an xval logfile such as xval_yeast.txt because each iteration uses
	# different training set and testing set from the reference complexes, and the xval
	# logfile logs which complexes (given by their complex id at complexes_yeast.txt) are being used as test complexes
	# (important for the precision metric which filters out predicted complexes that match training complexes)
	my %test_complexes = ();	# equivalent to C
	my $xval_iters = ReadXValData($xvalfilename, \%test_complexes);	
		
	# count total number of test complexes
	my $num_test_complexes = 0;
	# Filter test complexes
	foreach my $iter (keys %test_complexes) {
		# consider all complexes
		if ($which_comps_consider==0) {
			# scalar keys %hash gives the number of test complexes in that hash iteration
			# we are adding up $num_test_complexes in each iter to get total num of test_complexes for verification purposes
			$num_test_complexes += scalar keys %{$test_complexes{$iter}};
		}
		# consider only large complexes (size >= 4)
		elsif ($which_comps_consider==1) {
			# 
			foreach my $comp (keys %{$test_complexes{$iter}}) {
				if (scalar keys %{$complexes{$comp}} <= 3) {
					# by using info from the hash complexes, delete small test_complexes
					delete $test_complexes{$iter}{$comp};
				}
				else {
					$num_test_complexes++;
				}
			}
		}
		# consider only small complexes (size 2 or 3)
		elsif ($which_comps_consider==23) {
			foreach my $comp (keys %{$test_complexes{$iter}}) {
				if (scalar keys %{$complexes{$comp}} != 2 && scalar keys %{$complexes{$comp}} != 3) {
					# by using info from the hash complexes, delete large test_complexes
					delete $test_complexes{$iter}{$comp};
				}
				else {
					$num_test_complexes++;
				}
			}
			
		} # end which_comps_consider==23
	}
	
	# print some statistics
	print "num test complexes = $num_test_complexes\n";
	foreach my $iter (sort keys %test_complexes) {
		print "Iteration $iter, num test complexes = ".(scalar keys %{$test_complexes{$iter}})."\n";
	}
	
	# Read the clusters i.e. predicted complexes P
	# $clusters{$iter}{$c}{SCORE} = score, $clusters{$iter}{$c}{ELEMENTS}{$pid} = 1
	# Data structure:
	# clusters = {iter => c => {SCORE, ELEMENTS => pid}}
	my $inputfilename = $argopts{'i'};
	my %clusters_orig = ();		# equivalent to P

	# Recall that $num_iters comes from -n, and if it is set, we iterate through the file name
	# which if it looks like "clusters integrated.txt" will iterate through "clusters integrated iter{0..num_iters}"
	ReadClustersIters($inputfilename, \%clusters_orig, $num_iters);
	# ReadClustersIters populates %cluster_orig with key-value pairs like this:
	# {1 => {C1_CCL1_C125|4|CL1|CMC|COACH|IPCA||2|1|2| => {SCORE => 1.50999445763554, ELEMENTS => {YER157W => 1, YGL005C => 1, ...}}, ...}}
	for (my $iter=0; $iter<$num_iters; $iter++) {
		# Filter P using -l -s and -u parameters
		FilterClusters($clusters_orig{$iter}, $which_comps_consider, $score_threshold, $filter_unique);
	}

	# P, C, and also T (for precision caveat) has already been generated here for each iteration. Time to evaluate.
	my %clusters = ();
	
	# ------------ Get Precision vs. Recall for matchscore = 0.5 -------------
	# make a working copy of clusters
	%clusters = ();
	for (my $iter=0; $iter<$num_iters; $iter++) {
		foreach my $clus (keys %{$clusters_orig{$iter}}) {
			$clusters{$iter}{$clus}{SCORE} = $clusters_orig{$iter}{$clus}{SCORE};
			foreach my $prot (keys %{$clusters_orig{$iter}{$clus}{ELEMENTS}}) {
				$clusters{$iter}{$clus}{ELEMENTS}{$prot} = $clusters_orig{$iter}{$clus}{ELEMENTS}{$prot};
			}
		}
	}
	my $avg_auc = 0;
	my $se_auc = 0;
	my $avg_rec = 0;
	my $se_rec = 0;	
	# $test_complexes_matched{$comp}{N} = num iters that $comp is tested in, {M} = num iters it is matched
	my %test_complexes_matched = (); 
	foreach my $iter (keys %test_complexes) {
		foreach my $comp (keys %{$test_complexes{$iter}}) {
			$test_complexes_matched{$comp}{"N"}++;
		}
	}
	print "*************** Match threshold = 0.5 *************\n";
	if ($SMALLCOMP_MATCH1==1) {
		print "*************** Match threshold = 1 for small comps *************\n";
	}
	for (my $iter=0; $iter<$num_iters; $iter++) {
		print "Iter $iter\n";
		my $rec = 0;
		my %matched_complexes = ();
		my $auc = CalcPrecRecCompPred(0.5, \%{$clusters{$iter}}, \%{$test_complexes{$iter}}, \%complexes, \$rec, \%matched_complexes);
		foreach my $comp (keys %matched_complexes) {
			$test_complexes_matched{$comp}{"M"}++;
		}
		$avg_auc += $auc;
		$se_auc += $auc ** 2;
		$avg_rec += $rec;
		$se_rec += $rec ** 2;
	}
	$avg_auc /= $num_iters;
	$se_auc /= $num_iters;
	$se_auc = $se_auc - $avg_auc**2;
	$se_auc = ($se_auc * $num_iters / ($num_iters-1)) ** .5;
	$se_auc = $se_auc / ($num_iters ** .5);
	
	$avg_rec /= $num_iters;
	$se_rec /= $num_iters;
	$se_rec = $se_rec - $avg_rec**2;
	$se_rec = ($se_rec * $num_iters / ($num_iters-1)) ** .5;
	$se_rec = $se_rec / ($num_iters ** .5);
	
	print "\nAUCmean\t$avg_auc\n";
	print "AUCse\t$se_auc\n";
	print "RECALLmean\t$avg_rec\n";
	print "RECALLse\t$se_rec\n";
	print "\n";

	foreach my $comp (keys %test_complexes_matched) {
		$test_complexes_matched{$comp}{"P"} = $test_complexes_matched{$comp}{"M"} / $test_complexes_matched{$comp}{"N"};
	}
	print "Test complexes matched:\n";
	print "Comp\tPercentMatched\tNumTimesMatched\tNumTimesTested\n";
	foreach my $comp (sort {$test_complexes_matched{$a}{"P"} <=>$test_complexes_matched{$b}{"P"}} keys %test_complexes_matched) {
		print "$comp\t$test_complexes_matched{$comp}{P}\t$test_complexes_matched{$comp}{M}\t$test_complexes_matched{$comp}{N}\n";
	}
	
	# print overall precision-recall for all iterations
	%clusters = ();
	for (my $iter=0; $iter<$num_iters; $iter++) {
		foreach my $clus (keys %{$clusters_orig{$iter}}) {
			$clusters{$iter}{$clus}{SCORE} = $clusters_orig{$iter}{$clus}{SCORE};
			foreach my $prot (keys %{$clusters_orig{$iter}{$clus}{ELEMENTS}}) {
				$clusters{$iter}{$clus}{ELEMENTS}{$prot} = $clusters_orig{$iter}{$clus}{ELEMENTS}{$prot};
			}
		}
	}
	CalcPrecRecCompPredAllIters (0.5, \%clusters, \%test_complexes, \%complexes, $num_iters);
	print "\n";
		
		
		
		
	# ------------ Get Precision vs. Recall for matchscore = 0.75 -------------
	my %clusters = ();
	for (my $iter=0; $iter<$num_iters; $iter++) {
		foreach my $clus (keys %{$clusters_orig{$iter}}) {
			$clusters{$iter}{$clus}{SCORE} = $clusters_orig{$iter}{$clus}{SCORE};
			foreach my $prot (keys %{$clusters_orig{$iter}{$clus}{ELEMENTS}}) {
				$clusters{$iter}{$clus}{ELEMENTS}{$prot} = $clusters_orig{$iter}{$clus}{ELEMENTS}{$prot};
			}
		}
	}
	my $avg_auc = 0;
	my $se_auc = 0;
	my $avg_rec = 0;
	my $se_rec = 0;
	# $test_complexes_matched{$comp}{N} = num iters that $comp is tested in, {M} = num iters it is matched
	my %test_complexes_matched = (); 
	foreach my $iter (keys %test_complexes) {
		foreach my $comp (keys %{$test_complexes{$iter}}) {
			$test_complexes_matched{$comp}{"N"}++;
		}
	}
	print "*************** Matchscore 0.75 *************\n";
	if ($SMALLCOMP_MATCH1==1) {
		print "*************** Match threshold = 1 for small comps *************\n";
	}
	for (my $iter=0; $iter<$num_iters; $iter++) {
		print "Iter $iter\n";
		my $rec = 0;
		my %matched_complexes = ();
		my $auc = CalcPrecRecCompPred(0.75, \%{$clusters{$iter}}, \%{$test_complexes{$iter}}, \%complexes, \$rec, \%matched_complexes);
		foreach my $comp (keys %matched_complexes) {
			$test_complexes_matched{$comp}{"M"}++;
		}
		$avg_auc += $auc;
		$se_auc += $auc ** 2;
		$avg_rec += $rec;
		$se_rec += $rec ** 2;
	}
	$avg_auc /= $num_iters;
	$se_auc /= $num_iters;
	$se_auc = $se_auc - $avg_auc**2;
	$se_auc = ($se_auc * $num_iters / ($num_iters-1)) ** .5;
	$se_auc = $se_auc / ($num_iters ** .5);
	$avg_rec /= $num_iters;
	$se_rec /= $num_iters;
	$se_rec = $se_rec - $avg_rec**2;
	$se_rec = ($se_rec * $num_iters / ($num_iters-1)) ** .5;
	$se_rec = $se_rec / ($num_iters ** .5);
	print "\nAUCmean\t$avg_auc\n";
	print "AUCse\t$se_auc\n";
	print "RECALLmean\t$avg_rec\n";
	print "RECALLse\t$se_rec\n";
	print "\n";

	foreach my $comp (keys %test_complexes_matched) {
		$test_complexes_matched{$comp}{"P"} = $test_complexes_matched{$comp}{"M"} / $test_complexes_matched{$comp}{"N"};
	}
	print "Test complexes matched:\n";
	print "Comp\tPercentMatched\tNumTimesMatched\tNumTimesTested\n";
	foreach my $comp (sort {$test_complexes_matched{$a}{"P"} <=>$test_complexes_matched{$b}{"P"}} keys %test_complexes_matched) {
		print "$comp\t$test_complexes_matched{$comp}{P}\t$test_complexes_matched{$comp}{M}\t$test_complexes_matched{$comp}{N}\n";
	}

	# print overall precision-recall for all iterations
	%clusters = ();
	for (my $iter=0; $iter<$num_iters; $iter++) {
		foreach my $clus (keys %{$clusters_orig{$iter}}) {
			$clusters{$iter}{$clus}{SCORE} = $clusters_orig{$iter}{$clus}{SCORE};
			foreach my $prot (keys %{$clusters_orig{$iter}{$clus}{ELEMENTS}}) {
				$clusters{$iter}{$clus}{ELEMENTS}{$prot} = $clusters_orig{$iter}{$clus}{ELEMENTS}{$prot};
			}
		}
	}
	CalcPrecRecCompPredAllIters (0.75, \%clusters, \%test_complexes, \%complexes, $num_iters);
	print "\n";
	
	
} # end for $num_iters > 0


 
 

sub ReadXValData ($$) {
	my $filename = $_[0];
	my $xvaldataref = $_[1]; # $xval_date{$iter}{$complexid} = 1, if $complexid is a TEST complex in $iter
	
	open (XVALFILE, "$filename") || die ("Cannot open xval sampling file $filename");
	my $curriter;
	foreach my $line (<XVALFILE>) {
		chomp $line;
		my @toks = split(/\t/, $line);
		# Get iteration from xval file
		if ($toks[0] eq "iter") {
			$curriter = $toks[1] + 0;
		}
		# Get complex ID from xval file, 
		elsif (defined $toks[0] && $toks[0] ne "") {
			my $compid = $toks[0] + 0;
			$$xvaldataref{$curriter}{$compid} = 1;
		}
	}
	close XVALFILE;
	
	# print some stats
	print "Read xval data:\n";
	foreach my $iter (sort {$a <=> $b} keys %{$xvaldataref}) {
		print "Iter $iter, ".(scalar keys %{$$xvaldataref{$iter}})." test complexes\n";
	}
	print "\n";
	
	return ($curriter+1);
}




# match cluster with a lower bound. If match score is definitely below lower bound, then just return 0, no need to compute exactly
sub MatchClustersLB ($$$) {
	my $clus1_ref = $_[0];
	my $clus2_ref = $_[1];
	my $lowerbound = $_[2];
	
	# require exact match for small comps
	if ($SMALLCOMP_MATCH1==1) {
		# both are small comps
		if (scalar keys %{$clus1_ref} <= 3 && scalar keys %{$clus2_ref} <= 3) {
			if (scalar keys %{$clus1_ref} != scalar keys %{$clus2_ref}) {
				return 0;
			}
			my $num_intersect = 0;
			foreach my $elem (keys %$clus2_ref) {
				if (defined $$clus1_ref{$elem}) {
					$num_intersect++;
				}
			}
			if ($num_intersect == scalar keys %{$clus1_ref}) {
				return 1;
			}
			else {
				return 0;
			}
		}
		# only one is small comp, return 0
		elsif (scalar keys %{$clus1_ref} <= 3) {
			return 0;
		}
		elsif (scalar keys %{$clus2_ref} <= 3) {
			return 0;
		}		
	}
	
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

 

# $clusters{$iter}{$c}{SCORE} = score, $clusters{$iter}{$c}{ELEMENTS}{$pid} = 1
sub ReadClustersIters ($$$) {
	my $filename = $_[0];
	my $clusters_ref = $_[1];
	my $numiters = $_[2];
	
	for (my $iter=0; $iter<$numiters; $iter++) {		
		open (CLUSTERS_FILE, "$filename iter".$iter.".txt") || die ("$filename iter".$iter.".txt");
		# Sample line: C1_CCL1_C125|4|CL1|CMC|COACH|IPCA||2|1|2|(8_1.50999445763554): YER157W YGL005C YGL223C YGR120C YML071C YNL041C YNL051W YPR105C
		foreach my $line (<CLUSTERS_FILE>) {
			chomp($line);
			my @toks = split(' ', $line);
			# toks = [C1_CCL1_C125|4|CL1|CMC|COACH|IPCA||2|1|2|(8_1.50999445763554):,YER157W,YGL005C,YGL223C,YGR120C,YML071C,YNL041C,YNL051W,YPR105C]
			(my $clus, my $rest) = split(/\(/, $toks[0]);
			# $clus = C1_CCL1_C125|4|CL1|CMC|COACH|IPCA||2|1|2|
			# $rest = 8_1.50999445763554):
			(my $rest2, my $score) = split("_", $rest);
			# $rest2 = 8
			# $score = 1.50999445763554)
			$score = substr($score, 0, length($score)-2);
			# $score = 1.50999445763554 (then typecast to number)
			$$clusters_ref{$iter}{$clus}{SCORE} = $score + 0;
			shift @toks;
			# $toks = [YER157W,YGL005C,YGL223C,YGR120C,YML071C,YNL041C,YNL051W,YPR105C]
			foreach (@toks) {
				$$clusters_ref{$iter}{$clus}{ELEMENTS}{$_} = 1;
			}
		}		
		close CLUSTERS_FILE;
	}	
#	foreach my $clus (keys %$clusters_ref) {
#		print "$clus\t$$clusters_ref{$clus}{SCORE}\n";
#		foreach my $elem (keys %{$$clusters_ref{$clus}{ELEMENTS}}) {
#			print "\t$elem\n";
#		}
#	}	
}


# Not used
# $clusters{$c}{SCORE} = score, $clusters{$c}{ELEMENTS}{$pid} = 1
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
			$$clusters_ref{$clus}{ELEMENTS}{$_} = 1;
		}
	}		
	close CLUSTERS_FILE;	
}






sub FilterClusters {
	my $clustersref = $_[0];
	my $which_comps_consider = $_[1];
	my $score_threshold = $_[2];
	my $filter_unique = $_[3];
	print "FilterClusters, num clusters originally = ".(scalar keys %{$clustersref})."\n";
	# Filter clusters
	if ($which_comps_consider == 1) {		
		my %clusters_to_remove = ();
		foreach my $clus (keys %{$clustersref}) {
			if (scalar keys %{$$clustersref{$clus}{ELEMENTS}} <= 3) {
				$clusters_to_remove{$clus} = 1;
			}
		}
		foreach my $clus (keys %clusters_to_remove) {
			delete $$clustersref{$clus};
		}
		print "keeping only clusters of size > 3, num clusters = ".(scalar keys %{$clustersref})."\n";
	}
	elsif ($which_comps_consider == 23) {	
		foreach my $clus (keys %{$clustersref}) {
			if (scalar keys %{$$clustersref{$clus}{ELEMENTS}} > 3) {
				delete $$clustersref{$clus};
			}
		}	
		print "keeping only clusters of size <= 3, num clusters = ".(scalar keys %{$clustersref})."\n";
	}
	if ($REMOVE_SMALLCLUS==1) {
		foreach my $clus (keys %{$clustersref}) {
			if (scalar keys %{$$clustersref{$clus}{ELEMENTS}} <= 3) {
				delete $$clustersref{$clus};
			}
		}
		print "Discarding small clusters, num clusters = ".(scalar keys %{$clustersref})."\n";
	}
	
	if ($score_threshold > 0) {
		my %clusters_to_remove = ();
		foreach my $clus (keys %{$clustersref}) {
			if ($$clustersref{$clus}{SCORE} < $score_threshold) {
				$clusters_to_remove{$clus} = 1;
			}
		}
		foreach my $clus (keys %clusters_to_remove) {
			delete $$clustersref{$clus};
		}
		print "keeping only clusters of score > $score_threshold, num clusters = ".(scalar keys %{$clustersref})."\n";
	}
	
	if ($filter_unique==1) {
		RemoveDuplicateClusters($clustersref, 0.5);
		print "keeping only unique clusters (match_thre=0.5), num clusters = ".(scalar keys %{$clustersref})."\n";
	}
	
	if ($NUM_CLUS_KEEP < scalar keys %{$clustersref}) {
		my @clus_sorted = sort {$$clustersref{$b}{SCORE} <=> $$clustersref{$a}{SCORE}} keys %{$clustersref};
		for (my $idx=$NUM_CLUS_KEEP; $idx < scalar @clus_sorted; $idx++) {
			delete $$clustersref{$clus_sorted[$idx]};
		}
		print "Keeping only top $NUM_CLUS_KEEP clusters, num clusters = ".(scalar keys %{$clustersref})."\n";
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
			my $matchscore = MatchClustersLB($$clusters_ref{$clus1}{ELEMENTS}, $$clusters_ref{$clus2}{ELEMENTS}, $match_thres);
			if ($matchscore >= $match_thres) {
				$clusters_to_delete{$clus2} = 1;
			}
		}
	}
		
	print "In RemoveDuplicateClusters, num to delete = ".(scalar keys %clusters_to_delete)."\n";
	
	foreach my $clus (keys %clusters_to_delete) {
		delete $$clusters_ref{$clus};
	}
	
}






sub CalcPrecRecCompPred {
	my $matchscore_thr = $_[0];
	my $clusters_ref = $_[1];   # $$clusters_ref{$clus}{SCORE} = $score, {ELEMENTS}{$prot} = 1
	my $test_comps_ref = $_[2]; # $$test_comps_ref{$comp} = 1  iff $comp is test complex in iter $iter
	my $all_comps_ref = $_[3];  # $$all_comps_ref{$comp}{$prot} = 1
	my $recall_ref = $_[4];
	my $matched_complexes_ref = $_[5];
	
	# $results[i][0] = score
	# $results[i][1] = number of correct matches
	# $results[i][2] = cluster id
	# $results[i][3]{complex_id} = 1 if matched to complex_id
	my @results;
	print "Num predicted clusters = ".(scalar keys %{$clusters_ref})."\n";
	
	# 1. For each cluster in %clusters:
	#    - Check if it matches a test complex. If so, put in @results as correct
	#    - If it has matched any test complexes, delete from %clusters
	# 2. For each remaining cluster in %clusters:
	#    - Check if it matches a nontest complex. If so, delete from %clusters
	# 3. For each remaining cluster in %clusters:
	#    - Put in @results as incorrect
	my %correct_clusters = ();
	my %correct_smallclusters = ();
	foreach my $clus (keys %{$clusters_ref}) {
		my @cluster_matches; 	# since each cluster may be matched to multiple complexes, this structure keeps track of them
		$cluster_matches[0] = $$clusters_ref{$clus}{SCORE};
		$cluster_matches[1] = 0;
		$cluster_matches[2] = $clus;
		$cluster_matches[3] = ();
		foreach my $comp (keys %{$test_comps_ref}) {
			my $matchscore = MatchClustersLB($$clusters_ref{$clus}{ELEMENTS}, $$all_comps_ref{$comp}, $matchscore_thr);
			if ($matchscore >= $matchscore_thr) {
				$correct_clusters{$clus} = 1;
				if (scalar keys %{$$clusters_ref{$clus}{ELEMENTS}} <= 3) {
					$correct_smallclusters{$clus} = 1;
				}
				$cluster_matches[1]++;
				$cluster_matches[3]{$comp} = 1;
				$$matched_complexes_ref{$comp} = 1;
			}
		}
		if ($cluster_matches[1]>0) {
			push (@results, \@cluster_matches);
		}
	}
	foreach my $clus (keys %correct_clusters) {
		delete $$clusters_ref{$clus};
	}
	print "\tNum correct clusters = ".(scalar keys %correct_clusters)."\n";
	print "\tNum correct small clusters = ".(scalar keys %correct_smallclusters)."\n";
	print "\tNum test complexes matched = ".(scalar keys %{$matched_complexes_ref})."\n";
#		foreach my $comp (keys %{$matched_complexes_ref}) {
#			print "\t\t$comp\n";
#		}
	print "\tNum clusters not correct (matched to nontest complex, or wrong) = ".(scalar keys %{$clusters_ref})."\n";
	my %nontest_correct_clusters = ();
	foreach my $clus (keys %{$clusters_ref}) {
		foreach my $comp (keys %{$all_comps_ref}) {
			if (defined $$test_comps_ref{$comp}) { next; }
			my $matchscore = MatchClustersLB($$clusters_ref{$clus}{ELEMENTS}, $$all_comps_ref{$comp}, $matchscore_thr);
			if ($matchscore >= $matchscore_thr) {
				$nontest_correct_clusters{$clus} = 1;
				last;
			}
		}
	}
	foreach my $clus (keys %nontest_correct_clusters) {
		delete $$clusters_ref{$clus};
	}
	print "\tNum clusters matched to nontest partial complex = ".(scalar keys %nontest_correct_clusters)."\n";
	print "\tNum wrong clusters = ".(scalar keys %{$clusters_ref})."\n";
	foreach my $clus (keys %{$clusters_ref}) {
		my @tmparray;
		$tmparray[0] = $$clusters_ref{$clus}{SCORE};
		$tmparray[1] = 0;
		$tmparray[2] = $clus;
		$tmparray[3] = ();
		push (@results, \@tmparray);
	}	
	
	
	my $numtestcomps = scalar keys %{$test_comps_ref};
	print "Num test comps = $numtestcomps\n";
	
#	print "Clusters (excluding those that match nontest complexes:\n";
#	print "Score\tNumMatches\tID\tNumMatches\tIter\n";
#	foreach (sort {$$b[0] <=> $$a[0]} @results) {
#		print "$$_[0]\t$$_[1]\t$$_[2]\t".(scalar keys %{$$_[3]})."\t$$_[4]\n";
#	}
	
	print "Precision-recall for match_threshold = $matchscore_thr:\n";	
	print "Score threshold\tPredictions\tRecall\tPrecision\n";
	# precision = $corrects / $predicts
	# recall = $matched / $complex_count
	my $predicts = 0;		
	my $corrects= 0;
	my $matched = 0;
	my %matched_complexes = ();
	my $lastvalue = -1;
	my $threshold = 0.01;
	my $recall = 0;
	my $auc = 0;
	my $auc_prevrecall = 0;
	foreach (sort {$$b[0] <=> $$a[0]} @results) {
		if ($lastvalue != $$_[0]) {
	   	$recall = $matched/$numtestcomps;
	   	if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)    {
	    	print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
	    	while($threshold < $recall) {
	     		$threshold += 0.01;
	    	}
			  # accumulate AUC
			  $auc += ($recall - $auc_prevrecall) * ($corrects/$predicts);
			  $auc_prevrecall = $recall;
	   	}
	   	$lastvalue = $$_[0];
	  }
	  if ($$_[1] > 0){
	  	$corrects++;
	  	foreach my $comp (keys %{$$_[3]}) {
	  		$matched_complexes{$comp} = 1;
	  	}
	  	$matched = scalar keys %matched_complexes;
	  }
	  $predicts++;
	}
	if ($predicts > 0) {
	  $recall = $matched/$numtestcomps;
		print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
	}
	else {
	  $recall = $matched/$numtestcomps;
		print "$lastvalue\t$predicts\t$recall\t0\n";
	}
	# accumulate AUC
	$auc += ($recall - $auc_prevrecall) * ($corrects/$predicts);
	print "\nAUC\t$auc\n\n";
	$$recall_ref = $recall;
	return $auc;
}
	 





sub CalcPrecRecCompPredAllIters {
	my $matchscore_thr = $_[0];
	my $clusters_ref = $_[1];   # $$clusters_ref{$iter}{$clus}{SCORE} = $score, {ELEMENTS}{$prot} = 1
	my $test_comps_ref = $_[2]; # $$test_comps_ref{$iter}{$comp} = 1  iff $comp is test complex in iter $iter
	my $all_comps_ref = $_[3];  # $$all_comps_ref{$comp}{$prot} = 1
	my $num_iters = $_[4];
	
	# $results[i][0] = score
	# $results[i][1] = number of correct matches
	# $results[i][2] = cluster id
	# $results[i][3]{complex_id} = 1 if matched to complex_id
	my @results;
	
	
	for (my $iter = 0; $iter<$num_iters; $iter++) {	
#		print "Iteration $iter, num predicted clusters = ".(scalar keys %{$$clusters_ref{$iter}})."\n";
		
		# 1. For each cluster in %clusters:
		#    - Check if it matches a test complex. If so, put in @results as correct
		#    - If it has matched any test complexes, delete from %clusters
		# 2. For each remaining cluster in %clusters:
		#    - Check if it matches a nontest complex. If so, delete from %clusters
		# 3. For each remaining cluster in %clusters:
		#    - Put in @results as incorrect
		my %correct_clusters = ();
		my %matched_complexes = (); # just for debugging
		foreach my $clus (keys %{$$clusters_ref{$iter}}) {
			my @cluster_matches; 	# since each cluster may be matched to multiple complexes, this structure keeps track of them
			$cluster_matches[0] = $$clusters_ref{$iter}{$clus}{SCORE};
			$cluster_matches[1] = 0;
			$cluster_matches[2] = $clus;
			$cluster_matches[3] = ();
			$cluster_matches[4] = $iter;
			foreach my $comp (keys %{$$test_comps_ref{$iter}}) {
				my $matchscore = MatchClustersLB($$clusters_ref{$iter}{$clus}{ELEMENTS}, $$all_comps_ref{$comp}, $matchscore_thr);
				if ($matchscore >= $matchscore_thr) {
					$correct_clusters{$clus} = 1;
					$cluster_matches[1]++;
					$cluster_matches[3]{$comp} = 1;
					$matched_complexes{$comp} = 1;
				}
			}
			if ($cluster_matches[1]>0) {
				push (@results, \@cluster_matches);
			}
		}
		foreach my $clus (keys %correct_clusters) {
			delete $$clusters_ref{$iter}{$clus};
		}
#		print "\tNum correct clusters = ".(scalar keys %correct_clusters)."\n";
#		print "\tNum test complexes matched = ".(scalar keys %matched_complexes)."\n";
#		print "\tNum clusters not correct (matched to nontest complex, or wrong) = ".(scalar keys %{$$clusters_ref{$iter}})."\n";
		my %nontest_correct_clusters = ();
		foreach my $clus (keys %{$$clusters_ref{$iter}}) {
			foreach my $comp (keys %{$all_comps_ref}) {
				if (defined $$test_comps_ref{$iter}{$comp}) { next; }
				my $matchscore = MatchClustersLB($$clusters_ref{$iter}{$clus}{ELEMENTS}, $$all_comps_ref{$comp}, $matchscore_thr);
				if ($matchscore >= $matchscore_thr) {
					$nontest_correct_clusters{$clus} = 1;
					last;
				}
			}
		}
		foreach my $clus (keys %nontest_correct_clusters) {
			delete $$clusters_ref{$iter}{$clus};
		}
#		print "\tNum clusters matched to nontest complex = ".(scalar keys %nontest_correct_clusters)."\n";
#		print "\tNum wrong clusters = ".(scalar keys %{$$clusters_ref{$iter}})."\n";
		foreach my $clus (keys %{$$clusters_ref{$iter}}) {
			my @tmparray;
			$tmparray[0] = $$clusters_ref{$iter}{$clus}{SCORE};
			$tmparray[1] = 0;
			$tmparray[2] = $clus;
			$tmparray[3] = ();
			$tmparray[4] = $iter;
			push (@results, \@tmparray);
		}	
	}
	
	my $numtestcomps = 0;
	for (my $iter = 0; $iter<$num_iters; $iter++) {	
		$numtestcomps += scalar keys %{$$test_comps_ref{$iter}};
	}
	
	print "Precision-recall across all iters for match_threshold = $matchscore_thr:\n";	
	if ($SMALLCOMP_MATCH1==1) {
		print "*************** Match threshold = 1 for small comps *************\n";
	}
	print "Score threshold\tPredictions\tRecall\tPrecision\n";
	# precision = $corrects / $predicts
	# recall = $matched / $complex_count
	my $predicts = 0;		
	my $corrects= 0;
	my $matched = 0;
	my %matched_complexes = ();
	my $lastvalue = -1;
	my $threshold = 0.01;
	my $recall = 0;
	my $auc = 0;
	my $auc_prevrecall = 0;
	foreach (sort {$$b[0] <=> $$a[0]} @results) {
		if ($lastvalue != $$_[0]) {
	   	$recall = $matched/$numtestcomps;
	   	if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)    {
	    	print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
	    	while($threshold < $recall) {
	     		$threshold += 0.01;
	    	}
			  # accumulate AUC
			  $auc += ($recall - $auc_prevrecall) * ($corrects/$predicts);
			  $auc_prevrecall = $recall;
	   	}
	   	$lastvalue = $$_[0];
	  }
	  if ($$_[1] > 0){
	  	$corrects++;
	  	foreach my $comp (keys %{$$_[3]}) {
	  		my $iter = $$_[4];
	  		$matched_complexes{$iter}{$comp} = 1;
	  	}
	  	$matched = 0;
	  	foreach my $iter (keys %matched_complexes) {
	  		$matched += scalar keys %{$matched_complexes{$iter}};
	  	}
	  }
	  $predicts++;
	}
	if ($predicts > 0) {
	  $recall = $matched/$numtestcomps;
		print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
	}
	else {
	  $recall = $matched/$numtestcomps;
		print "$lastvalue\t$predicts\t$recall\t0\n";
	}
	# accumulate AUC
	$auc += ($recall - $auc_prevrecall) * ($corrects/$predicts);
	print "\nAUC\t$auc\n\n";
}
	 

