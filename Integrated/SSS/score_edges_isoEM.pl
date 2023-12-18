###############################################################################
#
# I use a data structure, "data", to contain instances for learning (ie. edges).
# Data contains only instances (edges) that have at least 1 data source.
# For each edge, if a data source is missing, then its value may be set to 0, or undefined.
#
# data
#  |- SUBTYPES
#  |   |- Pubmed
#  |   |   |- TYPE
#  |   |       |- real
#  |   |- L12
#  |       |- TYPE
#  |           |- real
#  |
#  |- NUMPROTEINS
#  |   |- 6150
#  |
#  |- NUMCOMPLEXEDGES
#  |   |- 11255
#  |
#  |- INSTANCES
#      |- SGD0000053|SGD0019021
#      |   |- DATA
#      |   |   |- Pubmed
#      |   |   |   |- 0.51
#      |   |   |- L12
#      |   |       |- 0.85
# 		 |	 |
# 		 |	 | - CLASS
#      |   | 	 |- ISCOMPLEX
#      |   |   |   |- 1
#			 |	 |	 |
#			 |	 |	 | - ISBORDER
#			 |	 |	 |   |- 1
#			 |	 |	 |
#      |   |   | - COMPLEXES
#      |   |       |- 78
#      |   |           |- 1
#      |   |       |- 103
#      |   |           |- 1
#      |   |
#      |   |- SCORE
#      |   |   |- 0.656503
#      |   |
#      |   |- ITER  (not used anymore)
#      |       |- 3
#      |           |- 1
#      |    
#      |- SGD0000053|SGD0019021
#   etc.


use strict;
use IO::Handle;
use POSIX;
use Getopt::Std;

	
use constant INSTANCES => 1; 
use constant SUBTYPES => 2; 
use constant DATA => 3;
use constant TYPE => 4;
use constant ITER => 5;
use constant SCORE => 6;
use constant LOG2 => log(2);


# inputs:
# -i <input data file>
# -c <complexes file>
# -m <mode: x (xval, do all iters), s (split xval, do from iters -y to -z), e (evaluate prec-rec, read files from -r), t (trainall)>
# -y <in mode s (split xval), first iteration to do xval>
# -z <in mode s (split xval), last iteration to do xval>
# -r <in mode e, base filename for xval results>
# -x <cross-validation file>
# -a <approach (sss (size-sensitive swc), sss_isoem, s (swc), t (topo), r (string), n (nowei)). Default = sss_isoem>
# -p <Only use PPI edges with reliability greater than this. Default 0.>
# -s <for SSS_ISOEM approach, how many EM iters, default = 2>
# -v <for SSS_ISOEM approach, use only neighbours above this score. Default = 1e-7>
# -h <for SSS/swc approach, smoothen likelihoods so that they are non-decreasing, default 0>
# -l <for SSS/SWC approach, MDL discretization factor, higher = fewer splits, default 1>
# -0 <for SSS approach, use this value for prior probability of class 0>
# -2 <for SSS approach, use this value for prior probability of class 2>
# -4 <for SSS approach, use this value for prior probability of class 4>
# -d <for SWC approach, when calculating posterior probability, the likelihood products are raised to this power. Default 1.>
# -t <for SWC approach, use posterior probability (0, default) or posterior ratio (1)>
# -e <1: count zero-edges (use zero edges as implicit edges), 0: only count zeros if it is in a data edge, 2 (NOT DONE!): ignore all zeros, 3: use zero edges as implicit edges, but bin 0 likelihoods always 1. default 0>
# -u <1: calculate likelihoods using cumulative probability, 0: calculate likelihoods at each feature value (default)>
# -o <output basefilename>
# perl score_edges.pl -i data_PPI,L2,COEXP,STR,LIT.txt -c complexes_CYC.txt -m x -x xval_.8_.2_25.txt -a s -o "swc xval_.8_.2"

my %argopts;
if (! getopts('i:b:c:m:y:z:r:x:h:l:a:s:v:p:t:d:e:o:u:0:2:4:', \%argopts)) { die "invalid arguments"; }

my $inputfilename = $argopts{'i'};
my $compfilename = $argopts{'c'};
my $xvalfilename;
my $xvalresultsbasefilename;
my $start_iter;
my $end_iter;
my $mode = "t";
if (defined $argopts{'m'}) {
	$mode = $argopts{'m'};
	if ($mode eq "x") {
		$xvalfilename = $argopts{'x'};
	}
	elsif ($mode eq "s") {
		$xvalfilename = $argopts{'x'};
		$start_iter = $argopts{'y'} + 0;
		$end_iter = $argopts{'z'} + 0;
	}
	elsif ($mode eq "e") {
		$xvalfilename = $argopts{'x'};
		$xvalresultsbasefilename = $argopts{'r'};
	}
}
my $approach = "sss_isoem";
if (defined $argopts{'a'}) {
	$approach = $argopts{'a'};
}
my $outputbasefilename = $argopts{'o'};
my $COUNT_ZERO_EDGES = 0;
if (defined $argopts{'e'}) {
	$COUNT_ZERO_EDGES = $argopts{'e'} + 0;
}

my $DEPENDENCE_FACTOR = 1;
if (defined $argopts{'d'}) {
	$DEPENDENCE_FACTOR = $argopts{'d'}+0;
}

my $USE_POSTERIOR_RATIO = 0;
if (defined $argopts{'t'}) {
	$USE_POSTERIOR_RATIO = $argopts{'t'} + 0;
}

my $ONLY_PPI_THR = 0;
if (defined $argopts{'p'}) {
	$ONLY_PPI_THR = $argopts{'p'} + 0;
}

my $MDL_ALPHA = 1;
if (defined $argopts{'l'}) {
	$MDL_ALPHA = $argopts{'l'} + 0;
}


my $SMOOTHEN_LIKELIHOODS = 0;
if (defined $argopts{'h'}) {
	$SMOOTHEN_LIKELIHOODS = $argopts{'h'} + 0;
}

my $CUMULATIVE_PROB_LIKELIHOODS = 0;
if (defined $argopts{'u'}) {
	$CUMULATIVE_PROB_LIKELIHOODS = $argopts{'u'} + 0;
}

my $NUM_ISOEMITER = 2;
if (defined $argopts{'s'}) {
	$NUM_ISOEMITER = $argopts{'s'}+0;
}
my $ISO_THR = 1e-7;
if (defined $argopts{'v'}) {
	$ISO_THR = $argopts{'v'}+0;
}

my %SSS_PRIORS = ();
$SSS_PRIORS{0} = -1;
$SSS_PRIORS{2} = -1;
$SSS_PRIORS{4} = -1;
if (defined $argopts{'0'}) {
	$SSS_PRIORS{0} = $argopts{'0'}+0;
}
if (defined $argopts{'2'}) {
	$SSS_PRIORS{2} = $argopts{'2'}+0;
}
if (defined $argopts{'4'}) {
	$SSS_PRIORS{4} = $argopts{'4'}+0;
}

# =========================================================================
# Feature selection , discretization, learn parameters using entire data set. 
if ($mode eq "t") {
	
	# $inputdata{SUBTYPES}{$subtype}{TYPE} = "real"
	# $inputdata{INSTANCES}{$key}{DATA}{$type} = $score, $inputdata{INSTANCES}{$key}{CLASS}{"COMPLEXES"}{$complexid} = 1/0, $inputdata{INSTANCES{$key}{CLASS}{"ISCOMPLEX"} = 1/0, {COCOMP2}, {COCOMP3}, {COCOMP4}
	# $inputdata{INSTANCES}{$key}{SCORE} = $score (predicted score)
	my %inputdata = ();
	my %complexes = ();
	
	print "Input data: $inputfilename\nComplex data: $compfilename\nApproach: $approach\nMode: $mode\nXval: $xvalfilename\n\n";
	ReadData(\%inputdata, \%complexes, $inputfilename, $compfilename);
	FilterDataByApproach(\%inputdata, $approach, \%complexes);
	
	# $likelihoods{$iter}{$subtype}{$bin}{"COMPLEX"} = complex edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"NONCOMP"} = noncomp edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"ALL"} = all edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"LIKERATIO"} = likelihood ratio of complex edge vs all edge
	# $likelihoods{$iter}{"PRIORCOMP"} = prior prob of complex edge vs all edge
	# if $iter is "ALL", this gives the likelihoods for entire data
	my %likelihoods = ();
			
	# $disc_cutpoints{$iter}{$subtype} is an array of cutpoints for this subtype
	# if $iter is "ALL", this gives cutpoints when using entire dataset
	my %disc_cutpoints = ();

		
	# Use Isolatedness with EM
	if ($approach eq "sss_isoem") {
		
		# feature selection 
		my %finalselectedfeats = ();
		my $done = 0;
		while ($done==0) {
			my %selectedfeats = ();
			LearnFeatureSelection(\%inputdata, \%selectedfeats, "sss");
			SelectFeatures(\%inputdata, \%selectedfeats, "sss", \%complexes);
			if (scalar keys %selectedfeats == scalar keys %finalselectedfeats) {
				$done = 1;
			}
			else {
				%finalselectedfeats = %selectedfeats;
			}
		}
		print "After feature selection, num instances = ".(scalar keys %{$inputdata{INSTANCES}}).", NUMCOMPLEXEDGES = $inputdata{NUMCOMPLEXEDGES}, NUMPROTEINS = $inputdata{NUMPROTEINS}, NUM_NONCOMP_BORDER = $inputdata{NUM_NONCOMP_BORDER}\n";
			
		# ----- Do discretization using training set -----
		LearnDiscretizationRanges(\%inputdata, \%{$disc_cutpoints{ALL}}, "sss");
		DiscretizeData(\%inputdata, \%{$disc_cutpoints{ALL}});
			
		# calculate likelihoods
		CalcLikelihoodsSSS(\%inputdata, \%{$likelihoods{ALL}});
		print "Calculated likelihoods\n\n";
		ScoreEdgesPosteriorProbSSS(\%inputdata, \%{$likelihoods{ALL}});
			
		# print scores
#		PrintScoresSSS(("$outputbasefilename isoemiter1 scored_edges.txt"),\%inputdata);
			
		# calculate isolatedness_2
		for (my $curr_isoemiter=2; $curr_isoemiter <= $NUM_ISOEMITER; $curr_isoemiter++) {
			
			my %neighbours = (); # $neighbours{$prot1}{$prot2} = 1, iff $prot1 and $prot2 are neighbours
			CalcNeighbours(\%inputdata, \%neighbours, $ISO_THR);
			my %triangles = (); # $ppitriangles{$prot1|$prot2|$prot3} = 1, iff they form a PPI triangle
			CalcTriangles(\%inputdata, \%triangles, $ISO_THR);
				
			CalcIso23combined(\%inputdata, \%neighbours, \%triangles, $ISO_THR);
			LearnDiscretizationRanges(\%inputdata, \%{$disc_cutpoints{ALL}}, "sss", "ISO2");
			DiscretizeData(\%inputdata, \%{$disc_cutpoints{ALL}}, "ISO2");
			CalcLikelihoodsSSS(\%inputdata, \%{$likelihoods{ALL}}, "ISO2");
			print "ISO EM Iter $curr_isoemiter, calculated likelihoods\n\n";
			ScoreEdgesPosteriorProbSSS(\%inputdata, \%{$likelihoods{ALL}});
			PrintScoresSSS(("$outputbasefilename isoemiter$curr_isoemiter scored_edges.txt"),\%inputdata);
		}
	}
		
	else {
	}
}

# =========================================================================
# print data
if ($mode eq "p") {
	
	# $inputdata{SUBTYPES}{$subtype}{TYPE} = "real"
	# $inputdata{INSTANCES}{$key}{DATA}{$type} = $score, $inputdata{INSTANCES}{$key}{CLASS}{"COMPLEXES"}{$complexid} = 1/0, $inputdata{INSTANCES{$key}{CLASS}{"ISCOMPLEX"} = 1/0
	# $inputdata{INSTANCES}{$key}{SCORE} = $score (predicted score)
	my %inputdata = ();
	my %complexes = ();
	my %num_comp_testedges = (); # $num_comp_testedges{$iter} = number of complex test edges in iteration $iter, including non-data edges
	
	print "Input data: $inputfilename\nComplex data: $compfilename\nApproach: $approach\nMode: $mode\nXval: $xvalfilename\n\n";
	ReadData(\%inputdata, \%complexes, $inputfilename, $compfilename);
#	PrintArff(("$outputbasefilename.arff"), \%inputdata);
	PrintTabDelimited("$outputbasefilename data.txt", \%inputdata);
#	FilterDataByApproach(\%inputdata, $approach, \%complexes);
	
}

# =========================================================================
# Do xvalidation
elsif ($mode eq "x") {
		# $inputdata{SUBTYPES}{$subtype}{TYPE} = "real"
	# $inputdata{INSTANCES}{$key}{DATA}{$type} = $score, $inputdata{INSTANCES}{$key}{CLASS}{"COMPLEXES"}{$complexid} = 1/0, $inputdata{INSTANCES{$key}{CLASS}{"ISCOMPLEX"} = 1/0
	# $inputdata{INSTANCES}{$key}{SCORE} = $score (predicted score)
	my %inputdata = ();
	my %complexes = (); # $complexes{$comp}{$prot} = 1
	
	print "Input data: $inputfilename\nComplex data: $compfilename\nApproach: $approach\nMode: $mode\nXval: $xvalfilename\n\n";
	ReadData(\%inputdata, \%complexes, $inputfilename, $compfilename);
#	my $total_edges_space = $inputdata{NUMPROTEINS} * ($inputdata{NUMPROTEINS}-1) * .5; # total number of edges regardless of approach used, for calculating ROC
	FilterDataByApproach(\%inputdata, $approach, \%complexes);
	
	# $likelihoods{$iter}{$subtype}{$bin}{"COMPLEX"} = complex edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"NONCOMP"} = noncomp edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"ALL"} = all edges likelihood for this subtype and bin
	# $likelihoods{$iter}{$subtype}{$bin}{"LIKERATIO"} = likelihood ratio of complex edge vs all edge
	# $likelihoods{$iter}{"PRIORCOMP"} = prior prob of complex edge vs all edge
	# if $iter is "ALL", this gives the likelihoods for entire data
	my %likelihoods = ();
			
	# $disc_cutpoints{$iter}{$subtype} is an array of cutpoints for this subtype
	# if $iter is "ALL", this gives cutpoints when using entire dataset
	my %disc_cutpoints = ();
		
	# $xval_data{$iter}{$complexid} = 1, if $complexid is a TEST complex in $iter
	my %xval_data = (); 
	my $xval_iters = ReadXValData($xvalfilename, \%xval_data);
	
#	# remove test complexes that are small
#	foreach my $iter (keys %xval_data) {
#		foreach my $comp (keys %{$xval_data{$iter}}) {
#			if (scalar keys %{$complexes{$comp}} <= 3) {
#				delete $xval_data{$iter}{$comp};
#			}
#		}
#	}
	
	for (my $iter=0; $iter<$xval_iters; $iter++) {
		
		print "---------------- Xval iteration $iter ---------------\n";
		
		# initialize %curr_train_complexes to the train complexes from xval_data
		my %curr_train_complexes = (); 
		foreach my $comp (keys %complexes) {
			if (!defined $xval_data{$iter}) { die; }
			if (!defined $xval_data{$iter}{$comp}) {
				foreach my $prot (keys %{$complexes{$comp}}) {
					$curr_train_complexes{$comp}{$prot} = 1;
				}
			}
		}
		
		my %currdata = ();
		CreateCurrData(\%inputdata, \%currdata, \%curr_train_complexes);
					
		
		# Use Isolatedness with EM
		if ($approach eq "sss_isoem") {
				
			# feature selection 
			my %finalselectedfeats = ();
			my $done = 0;
			while ($done==0) {
				my %selectedfeats = ();
				LearnFeatureSelection(\%currdata, \%selectedfeats, "sss");
				SelectFeatures(\%currdata, \%selectedfeats, "sss", \%curr_train_complexes);
				if (scalar keys %selectedfeats == scalar keys %finalselectedfeats) {
					$done = 1;
				}
				else {
					%finalselectedfeats = %selectedfeats;
				}
			}
			print "After feature selection, num instances = ".(scalar keys %{$currdata{INSTANCES}}).", NUMCOMPLEXEDGES = $currdata{NUMCOMPLEXEDGES}, NUMPROTEINS = $currdata{NUMPROTEINS}, NUM_NONCOMP_BORDER = $currdata{NUM_NONCOMP_BORDER}\n";
			
			# ----- Do discretization using training set -----
			LearnDiscretizationRanges(\%currdata, \%{$disc_cutpoints{$iter}}, "sss");
			DiscretizeData(\%currdata, \%{$disc_cutpoints{$iter}});
			
			# calculate likelihoods
			CalcLikelihoodsSSS(\%currdata, \%{$likelihoods{$iter}});
			print "Calculated likelihoods\n\n";
			ScoreEdgesPosteriorProbSSS(\%currdata, \%{$likelihoods{$iter}});
			
			# print scores
#			PrintScoresSSS(("$outputbasefilename isoemiter1 scored_edges iter".$iter.".txt"),\%currdata);
			
			# calculate isolatedness_2
			for (my $curr_isoemiter=2; $curr_isoemiter <= $NUM_ISOEMITER; $curr_isoemiter++) {
			
				my %neighbours = (); # $neighbours{$prot1}{$prot2} = 1, iff $prot1 and $prot2 are neighbours
				CalcNeighbours(\%currdata, \%neighbours, $ISO_THR);
				my %triangles = (); # $ppitriangles{$prot1|$prot2|$prot3} = 1, iff they form a PPI triangle
				CalcTriangles(\%currdata, \%triangles, $ISO_THR);
				
				CalcIso23combined(\%currdata, \%neighbours, \%triangles, $ISO_THR);
#				CalcIso3(\%currdata, \%neighbours, \%triangles, $ISO_THR);
				LearnDiscretizationRanges(\%currdata, \%{$disc_cutpoints{$iter}}, "sss", "ISO2");
#				LearnDiscretizationRanges(\%currdata, \%{$disc_cutpoints{$iter}}, "sss", "ISO3");
				DiscretizeData(\%currdata, \%{$disc_cutpoints{$iter}}, "ISO2");
#				DiscretizeData(\%currdata, \%{$disc_cutpoints{$iter}}, "ISO3");
				CalcLikelihoodsSSS(\%currdata, \%{$likelihoods{$iter}}, "ISO2");
#				CalcLikelihoodsSSS(\%currdata, \%{$likelihoods{$iter}}, "ISO3");
				print "ISO EM uter $curr_isoemiter, calculated likelihoods\n\n";
				ScoreEdgesPosteriorProbSSS(\%currdata, \%{$likelihoods{$iter}});
				PrintScoresSSS(("$outputbasefilename isoemiter$curr_isoemiter scored_edges iter".$iter.".txt"),\%currdata);
			}
			
		} # end Use Isolatedness with EM
		
		# Approach is not Isolatedness with EM
		else {
		
			# ----- Do feature selection using training set -----
			if ($approach eq "s" || $approach eq "sss") {
				my %finalselectedfeats = ();
				my $done = 0;
				while ($done==0) {
					my %selectedfeats = ();
					LearnFeatureSelection(\%currdata, \%selectedfeats, $approach);
					SelectFeatures(\%currdata, \%selectedfeats, $approach, \%curr_train_complexes);
					if (scalar keys %selectedfeats == scalar keys %finalselectedfeats) {
						$done = 1;
					}
					else {
						%finalselectedfeats = %selectedfeats;
					}
				}
				print "After feature selection, num instances = ".(scalar keys %{$currdata{INSTANCES}}).", NUMCOMPLEXEDGES = $currdata{NUMCOMPLEXEDGES}, NUMPROTEINS = $currdata{NUMPROTEINS}, NUM_NONCOMP_BORDER = $currdata{NUM_NONCOMP_BORDER}\n";
			#	PrintArff(("$outputfilebasename iter".$iter." train".".arff"), \%traindata);
			#	PrintArff(("$outputfilebasename iter".$iter." test".".arff"), \%testdata);
			
				# ----- Do discretization using training set -----
				LearnDiscretizationRanges(\%currdata, \%{$disc_cutpoints{$iter}}, $approach);
				DiscretizeData(\%currdata, \%{$disc_cutpoints{$iter}});
		#		PrintArff(("$outputfilebasename iter".$iter." disc train".".arff"), \%traindata);
		#		PrintArff(("$outputfilebasename iter".$iter." disc test".".arff"), \%testdata);
		
		
				# ----- Learn likelihood parameters using training set -----
				if ($approach eq "sss") {
					CalcLikelihoodsSSS(\%currdata, \%{$likelihoods{$iter}});
					print "Calculated likelihoods\n\n";
					ScoreEdgesPosteriorProbSSS(\%currdata, \%{$likelihoods{$iter}});
					print "Scored edges\n\n";
				}
				else {
					CalcLikelihoods(\%currdata, \%{$likelihoods{$iter}});
					print "Calculated likelihoods\n\n";
					if ($USE_POSTERIOR_RATIO==0) {
						ScoreEdgesPosteriorProb(\%currdata, \%{$likelihoods{$iter}});
					}
					else {
						ScoreEdgesPosteriorRatio(\%currdata, \%{$likelihoods{$iter}});
					}
					print "Scored edges\n\n";
				}
				
			}
		
			elsif ($approach eq "t") {
				ScoreEdgesPPIValue(\%currdata);
				print "Scored edges\n\n";
			}
			elsif ($approach eq "r") {
				ScoreEdgesStringValue(\%currdata);
				print "Scored edges\n\n";
			}
			elsif ($approach eq "n") {
				ScoreEdgesNoWeight(\%currdata);
				print "Scored edges\n\n";
			}
		#	elsif ($approach eq "stringonly") {
		#		ScoreEdgesString(\%testdata);
		#		print "Scored edges\n\n";
		#	}
			
			if ($approach eq "sss") {
				PrintScoresSSS(("$outputbasefilename scored_edges iter".$iter.".txt"),\%currdata);
			}
			else {
				# ----- Print scored edges for CMC -----
				PrintScoresCMC(("$outputbasefilename scored_edges iter".$iter.".txt"),\%currdata);		
	#			PrintEdgesLikelihoodRatios(("$outputbasefilename likeratios iter".$iter.".txt"), \%currdata, \%{$likelihoods{$iter}});
				
					
#				# ----- Save test edges' scores for cross-validation -----
#				# $xval_results[i][0] = score
#				# $xval_results[i][1] = complex edge?
#				my @xval_results;
#				foreach my $key (keys %{$currdata{INSTANCES}}) {
#					# don't evaluate on edges that are known to be co-complex
#					if ($currdata{INSTANCES}{$key}{CLASS}{ISCOMPLEX} == 1) { next; }
#					my @tmparray;
#					$tmparray[0] = $currdata{INSTANCES}{$key}{SCORE};
#					$tmparray[1] = $inputdata{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"};
#					push (@xval_results, \@tmparray);
#				}
#				print "Iteration $iter done\n\n";
#				
#				# set $num_comp_testedges = number of complex test edges, including non-data edges
#				my %comptrainedges = (); # the test edges should exclude training co-complex edges
#				foreach my $comp (keys %complexes) {
#					if (defined $xval_data{$iter}{$comp}) { next; }
#					foreach my $prot1 (keys %{$complexes{$comp}}) {
#						foreach my $prot2 (keys %{$complexes{$comp}}) {
#							if ($prot1 ge $prot2) { next; }
#							$comptrainedges{"$prot1|$prot2"} = 1;
#						}
#					}
#				}
#				my %comptestedges = (); # can't simply use number of prots in test comps, because different test comps may have some overlap in edges
#				foreach my $comp (keys %{$xval_data{$iter}}) {
#					foreach my $prot1 (keys %{$complexes{$comp}}) {
#						foreach my $prot2 (keys %{$complexes{$comp}}) {
#							if ($prot1 ge $prot2) { next; }
#							if (defined $comptrainedges{"$prot1|$prot2"}) { next; }
#							$comptestedges{"$prot1|$prot2"} = 1;
#						}
#					}
#				}
#				my $num_comp_testedges = scalar keys %comptestedges;
#				
#				# set $num_negative_edges = total number of edges - number of complex edges
#				my %allcompedges = ();
#				foreach my $comp (keys %complexes) {
#					foreach my $prot1 (keys %{$complexes{$comp}}) {
#						foreach my $prot2 (keys %{$complexes{$comp}}) {
#							if ($prot1 ge $prot2) { next; }
#							$allcompedges{"$prot1|$prot2"} = 1;
#						}
#					}
#				}
#				my $num_negative_edges = $total_edges_space - scalar keys %allcompedges;
#				
#				# print xval_results 
#				open (XVALRESULTS_FILE, ">$outputbasefilename xvalresults iter$iter.txt");
#				print XVALRESULTS_FILE "$num_comp_testedges\n";
#				print XVALRESULTS_FILE "$num_negative_edges\n";
#				foreach (@xval_results) {
#					print XVALRESULTS_FILE "$$_[0]\t$$_[1]\n";
#				}
#				close XVALRESULTS_FILE;
			}
		} # end Approach is not Isolatedness with EM
				
		
		
		# clear memory (sometimes run out of memory)
		undef %currdata;
	} # end foreach xval_iter
	#
	## print likelihoods learned in all iterations
	#PrintLikelihoods(\%likelihoods, \%disc_cutpoints, "ITERS");
	#
	
	
} # end mode eq "x"





# read from xval results file. 
elsif ($mode eq "e") {
		
	# $xval_data{$iter}{$complexid} = 1, if $complexid is a TEST complex in $iter
	my %xval_data = (); 
	my $xval_iters = ReadXValData($xvalfilename, \%xval_data);
	my %xval_results_iters; # $xval_results_iters{$iter}[i][0] = score, [1] = iscomplex?
	my %num_complex_edges_in_test_iters; # $num_complex_edges_in_test_iters{$iter} = num test complex edges in iter $iter
	my %num_neg_edges_iters; # $num_neg_edgeS_iters{$iter} = num negative edges in iter $iter
	
	ReadXValResultsIters($xvalresultsbasefilename, \%xval_results_iters, \%num_complex_edges_in_test_iters, \%num_neg_edges_iters, $xval_iters);
	
	my $avg_auc = 0;
	my $se_auc = 0;
	my $avg_aucroc = 0;
	my $se_aucroc = 0;
	for (my $iter=0; $iter<$xval_iters; $iter++) {
		print "Iter $iter:\n";
		my $auc = 0;
		my $fscore = 0;
		my $aucroc = 0;
		CalcPrecVsRecall (\@{$xval_results_iters{$iter}}, $num_complex_edges_in_test_iters{$iter}, $num_neg_edges_iters{$iter}, \$auc, \$fscore, \$aucroc);
		$avg_auc += $auc;
		$se_auc += $auc ** 2;
		$avg_aucroc += $aucroc;
		$se_aucroc += $aucroc**2;
	}
	$avg_auc /= $xval_iters;
	$se_auc /= $xval_iters;
	$se_auc = $se_auc - $avg_auc**2;
	$se_auc = ($se_auc * $xval_iters / ($xval_iters-1)) ** .5;
	$se_auc = $se_auc / ($xval_iters ** .5);
	$avg_aucroc /= $xval_iters;
	$se_aucroc /= $xval_iters;
	$se_aucroc = $se_aucroc - $avg_aucroc**2;
	$se_aucroc = ($se_aucroc * $xval_iters / ($xval_iters-1)) ** .5;
	$se_aucroc = $se_aucroc / ($xval_iters ** .5);
	
	print "\n";
	print "AUCmean\t$avg_auc\n";
	print "AUCse\t$se_auc\n";
	print "AUCROCmean\t$avg_aucroc\n";
	print "AUCROCse\t$se_aucroc\n\n";
	
	my $auc = 0;
	my $fscore = 0;
	my $aucroc = 0;
	CalcPrecVsRecallAllIters (\%xval_results_iters, \%num_complex_edges_in_test_iters, \%num_neg_edges_iters, \$auc, \$fscore, \$aucroc);

#	# calculate overal precision-recall graph
#	my @prec_rec_results_overall = (); # @prec_rec_results_overall[i]{REC}, {PREC}, {N}, where recalls should be in steps of threshold
#	for (my $iter=0; $iter<$xval_iters; $iter++) {
#		my $threshold = 0.01;
#		for (my $idx=0; $idx<scalar @{$prec_rec_results{$iter}}; $idx++) {
#			if ($prec_rec_results{$iter}[$idx]{REC} >= $threshold) {
#				my %tmphash = ();
#				$tmphash{$REC} = 
#				$prec_rec_results_overall[$idx]{REC} = $threshold;
#				$prec_rec_results_overall[$idx]{PREC} += $prec_rec_results{$iter}[$idx]{PREC};
#			
#	
	
	
}





# ================ Data utility functions (data reading, filtering data, printing, etc) =====================================
#--------------------- Read raw data -----------------------
sub ReadData($$$$) {

	my $dataref = $_[0];
	my $complexesref = $_[1];
	my $inputfilename = $_[2];
	my $complexfilename = $_[3];
	
	my %all_proteins = ();

	# $$dataref{SUBTYPES}{$subtype}{TYPE} = "real"
	# $$dataref{INSTANCES}{$key}{DATA}{$type} = $score, $$dataref{INSTANCES}{$key}{CLASS}{"COMPLEXES"}{$complexid} = 1/0, $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"} = 1/0, {COCOMP2}, {COCOMP3}, {COCOMP4}

	open(INPUTFILE, "$inputfilename") || die ("Cannot open input data file $inputfilename");
	open(COMPLEXFILE, "$complexfilename") || die ("Cannot open complex file $complexfilename");

	foreach my $line (<INPUTFILE>) { 
		chomp($line);
	  (my $id_a, my $id_b, my $subtype, my $score) = split(/\t/, $line);
	  if ($subtype eq '') {
	  	$subtype = 'Undefined';
	  }
	  if ($score eq '') {
	  	$score = 1;
	  }
	  if ($score eq "NaN") {
	  	$score = 0;
	  }
	  $id_a = uc($id_a);
	  $id_b = uc($id_b);
	  $all_proteins{$id_a} = 1;
	  $all_proteins{$id_b} = 1;
	  if ($id_a eq $id_b) {
	  	next;
	  }
	  if ($id_a eq '' || $id_b eq '') {
	   	die ($line);
	  }
	  my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
	  	  
	  # if this interaction & subtype already exists, keep the higher absolute one
		if (defined $$dataref{INSTANCES}{$key} && defined $$dataref{INSTANCES}{$key}{DATA} && defined $$dataref{INSTANCES}{$key}{DATA}{$subtype} && abs($$dataref{INSTANCES}{$key}{DATA}{$subtype}) > abs($score)) {
			next;
		}
	
		# add the interaction
	  $$dataref{SUBTYPES}{$subtype}{TYPE} = "real";
	  $$dataref{INSTANCES}{$key}{DATA}{$subtype} = $score;
	  $$dataref{INSTANCES}{$key}{CLASS}{"COMPLEXES"} = ();
	  $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"} = 0;
	  $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS} = 0;
	}
	close(INPUTFILE);
	
	# remove non-PPI edges
	foreach my $key(keys %{$$dataref{INSTANCES}}) {
		if (!defined $$dataref{INSTANCES}{$key}{DATA}{PPIREL} || $$dataref{INSTANCES}{$key}{DATA}{PPIREL} == 0) {
			delete $$dataref{INSTANCES}{$key};
		}		
	}
	
	# remove edges whose PPI reliability is not above the given threshold
	if ($ONLY_PPI_THR != 0) {
		foreach my $key(keys %{$$dataref{INSTANCES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{PPIREL} || $$dataref{INSTANCES}{$key}{DATA}{PPIREL} < $ONLY_PPI_THR) {
				delete $$dataref{INSTANCES}{$key};
			}		
		}
	}
	
	# ------------- Set 'missing' attributes of all instances to 0 ------
	foreach my $key(keys %{$$dataref{INSTANCES}}) {
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$$dataref{INSTANCES}{$key}{DATA}{$subtype} = 0;
			}
		}
	}
  
	#---------------------  read complexes -----------------------
	foreach my $line (<COMPLEXFILE>) {
		chomp($line);
		(my $id, my $comp) = split(/\t/, $line);
		$id = uc($id);
		$$complexesref{$comp}{$id} = 1;
	}
	my %all_complex_edges = ();
	my %all_cocomp2_edges = ();
	my %all_cocomp4_edges = ();
	foreach my $complex_id (keys %$complexesref) {
		my $compsize = scalar keys %{$$complexesref{$complex_id}};
		my $cocomp_size_class = ($compsize==2 || $compsize==3) ? 2 : 4;
		foreach my $complex_int1 (keys %{$$complexesref{$complex_id}}) {
			foreach my $complex_int2 (keys %{$$complexesref{$complex_id}}) {
	  
	  		$all_proteins{$complex_int1} = 1;
	  		$all_proteins{$complex_int2} = 1;
	  		
				if ($complex_int1 ne $complex_int2) {
					my $key = ($complex_int1 lt $complex_int2)?"$complex_int1|$complex_int2":"$complex_int2|$complex_int1";
					$all_complex_edges{$key} = 1;
					if ($cocomp_size_class==4) { 
						$all_cocomp4_edges{$key} = 1; 
						delete $all_cocomp2_edges{$key};
					}
					elsif ($cocomp_size_class==2 && !defined $all_cocomp4_edges{$key}) {
						$all_cocomp2_edges{$key} = 1; 
					}

					if (defined $$dataref{INSTANCES}{$key}) {					
						$$dataref{INSTANCES}{$key}{CLASS}{"COMPLEXES"}{$complex_id} = 1;
						$$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"} = 1;
						if ($cocomp_size_class > $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS}) {
							$$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS} = $cocomp_size_class;
						}
					}
				}
			}
		}
	}
	close(COMPLEXFILE);
	$$dataref{NUMPROTEINS} = scalar keys %all_proteins;
	$$dataref{NUMCOMPLEXEDGES} = scalar keys %all_complex_edges;
	$$dataref{NUMCOCOMP2EDGES} = scalar keys %all_cocomp2_edges;
	$$dataref{NUMCOCOMP4EDGES} = scalar keys %all_cocomp4_edges;
	

	#--------------- Collect information about data ------------
	my $total_possible_edges = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
	my $total_possible_complex_edges = $$dataref{NUMCOMPLEXEDGES};
	my $total_data_edges = scalar keys %{$$dataref{INSTANCES}};
	my $total_complex_data_edges = 0;
	my $total_noncomplex_data_edges = 0;
	my %subtypes_counts = (); # $subtypes_counts{$subtype}{NUMEDGES}, $subtypes_counts{$subtype}{PROTS}, $subtypes_counts{$subtype}{NUMCOMPEDGES}
	my %data_prots = ();
	foreach my $key(keys %{$$dataref{INSTANCES}}) {
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) { die; }
			if ($$dataref{INSTANCES}{$key}{DATA}{$subtype}!=0) {
				$subtypes_counts{$subtype}{NUMEDGES}++;
				(my $id_a, my $id_b) = split(/\|/, $key);
				$subtypes_counts{$subtype}{PROTS}{$id_a} = 1;
				$subtypes_counts{$subtype}{PROTS}{$id_b} = 1;
				if ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
					$subtypes_counts{$subtype}{NUMCOMPEDGES}++;
				}
			}
		}
		(my $id_a, my $id_b) = split(/\|/, $key);
		$data_prots{$id_a} = 1;
		$data_prots{$id_b} = 1;
		if ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
			$total_complex_data_edges++;
		}
		elsif ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==0) {
			$total_noncomplex_data_edges++;
		}
		else {
			die("Error: edge that is not data nor complex!\n");
		}
	}
	my $total_complex_nondata_edges = $total_possible_complex_edges - $total_complex_data_edges;
	my $num_complexes_gte4 = 0;
	foreach my $comp (keys %$complexesref) {
		if (scalar keys %{$$complexesref{$comp}} >= 4) {
			$num_complexes_gte4 ++;
		}
	}

	my %smallcomp_edges = ();
	my %largecomp_edges = ();
	foreach my $complex_id (keys %$complexesref) {
		my $compsize = scalar keys %{$$complexesref{$complex_id}};
		foreach my $complex_int1 (keys %{$$complexesref{$complex_id}}) {
			foreach my $complex_int2 (keys %{$$complexesref{$complex_id}}) {
	  		if ($complex_int1 ge $complex_int2) { next; }
				my $key = "$complex_int1|$complex_int2";
				if ($compsize <= 3) {
					$smallcomp_edges{$key} = 1;
				}
				else {
					$largecomp_edges{$key} = 1;
				}
			}
		}
	}
	my $num_small_large_overlap = 0;
	foreach my $edge (keys %smallcomp_edges) {
		if (defined $largecomp_edges{$edge}) {
			$num_small_large_overlap++;
		}
	}
	print "Read data\n";
	print "Total data edges = $total_data_edges\n";
	print "Total complex data edges = $total_complex_data_edges\n";
	print "Total noncomplex data edges = $total_noncomplex_data_edges\n";
	print "Total possible edges (including nondata edges) = $total_possible_edges\n";	
	print "Total possible complex edges (including nondata complex edges) = $total_possible_complex_edges\n";	
	print "Total possible cocomp2 edges (including nondata complex edges) = $$dataref{NUMCOCOMP2EDGES}\n";	
	print "Total possible cocomp4 edges (including nondata complex edges) = $$dataref{NUMCOCOMP4EDGES}\n";	
	print "Num edges both small and large comp = $num_small_large_overlap\n";
	print "Total complex non-data edges = $total_complex_nondata_edges\n";
	print "Num prots in edges = ".(scalar keys %data_prots)."\n";
	print "Num total prots (incl those without edges) = $$dataref{NUMPROTEINS}\n";
	print "Num complexes = ".(scalar keys %$complexesref)."\n";
	print "Num complexes size >= 4 = $num_complexes_gte4\n";
	print "Subtype counts:\n";
	foreach my $subtype (sort keys %subtypes_counts) {
		print "$subtype\tNum edges: $subtypes_counts{$subtype}{NUMEDGES}\tNum comp edges: $subtypes_counts{$subtype}{NUMCOMPEDGES}\tPrecision (% comp): ".($subtypes_counts{$subtype}{NUMCOMPEDGES}/$subtypes_counts{$subtype}{NUMEDGES})."\tRecall (% comp covered): ".($subtypes_counts{$subtype}{NUMCOMPEDGES}/$total_possible_complex_edges)."\tNum prots: ".(scalar keys %{$subtypes_counts{$subtype}{PROTS}})."\n";
	}
	print "OVERALL\tNum edges: $total_data_edges\tNum comp edges: $total_complex_data_edges\tPrecision (% comp): ".($total_complex_data_edges/$total_data_edges)."\tRecall (% comp covered): ".($total_complex_data_edges/$total_possible_complex_edges)."\tNum prots: ".(scalar keys %data_prots)."\n";
	print "\n";
		
	
#	my $total_edges = scalar keys %{$$dataref{INSTANCES}};
#	my $total_complex_edges = 0;
#	my $total_data_edges = 0;
#	my $total_complex_data_edges = 0;
#	my $total_noncomplex_data_edges = 0;
#	my $total_nondata_complex_edges = 0;
#	my %subtypes_counts = (); # $subtypes_counts{$subtype}{NUMEDGES}, $subtypes_counts{$subtype}{PROTS}, $subtypes_counts{$subtype}{NUMCOMPEDGES}
#	my %all_prots = ();
#	foreach my $key(keys %{$$dataref{INSTANCES}}) {
#		my $has_data = 0;
#		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
#			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype} && $$dataref{INSTANCES}{$key}{DATA}{$subtype}!=0) {
#				$has_data = 1;
#				$subtypes_counts{$subtype}{NUMEDGES}++;
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				$subtypes_counts{$subtype}{PROTS}{$id_a} = 1;
#				$subtypes_counts{$subtype}{PROTS}{$id_b} = 1;
#				if ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
#					$subtypes_counts{$subtype}{NUMCOMPEDGES}++;
#				}
#			}
#		}
#		if ($has_data==1) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			$all_prots{$id_a} = 1;
#			$all_prots{$id_b} = 1;
#		}
#		if ($has_data==1 && $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
#			$total_complex_edges++;
#			$total_data_edges++;
#			$total_complex_data_edges++;
#		}
#		elsif ($has_data==1 && $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==0) {
#			$total_data_edges++;
#			$total_noncomplex_data_edges++;
#		}
#		elsif ($has_data==0 && $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
#			$total_complex_edges++;
#			$total_nondata_complex_edges++;
#		}
#		else {
#			die("Error: edge that is not data nor complex!\n");
#		}
#	}
#	my $num_complexes_gte4 = 0;
#	foreach my $comp (keys %$complexesref) {
#		if (scalar keys %{$$complexesref{$comp}} >= 4) {
#			$num_complexes_gte4 ++;
#		}
#	}
#	print "Read data\n";
#	print "Total edges stored = $total_edges\n";
#	print "Total complex edges stored = $total_complex_edges\n";
#	print "Total data edges = $total_data_edges\n";
#	print "Total complex data edges = $total_complex_data_edges\n";
#	print "Total noncomplex data edges = $total_noncomplex_data_edges\n";
#	print "Total nondata complex edges = $total_nondata_complex_edges\n\n";
#	print "Num prots in edges = ".(scalar keys %all_prots)."\n";
#	print "Num total prots (incl those without edges) = $$dataref{NUMPROTEINS}\n";
#	print "Num total possible edges = ".($$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2)."\n";
#	print "Num total complex edges (incl those not stored) = $$dataref{NUMCOMPLEXEDGES}\n";
#	print "Num complexes = ".(scalar keys %$complexesref)."\n";
#	print "Num complexes size >= 4 = $num_complexes_gte4\n";
#	print "Subtype counts:\n";
##	foreach my $subtype (sort keys %subtypes_counts) {
##		print "$subtype\tNum edges: $subtypes_counts{$subtype}{NUMEDGES}\tNum comp edges: $subtypes_counts{$subtype}{NUMCOMPEDGES}\tPercent comp: ".($subtypes_counts{$subtype}{NUMCOMPEDGES}/$subtypes_counts{$subtype}{NUMEDGES})."\tNum prots: ".(scalar keys %{$subtypes_counts{$subtype}{PROTS}})."\n";
##	}
#		
#	foreach my $subtype (sort keys %subtypes_counts) {
#		print "$subtype\tNum edges: $subtypes_counts{$subtype}{NUMEDGES}\tNum comp edges: $subtypes_counts{$subtype}{NUMCOMPEDGES}\tPrecision (% comp): ".($subtypes_counts{$subtype}{NUMCOMPEDGES}/$subtypes_counts{$subtype}{NUMEDGES})."\tRecall (% comp covered): ".($subtypes_counts{$subtype}{NUMCOMPEDGES}/$total_possible_complex_edges)."\tNum prots: ".(scalar keys %{$subtypes_counts{$subtype}{PROTS}})."\n";
#	}
#	print "OVERALL\tNum edges: $total_data_edges\tNum comp edges: $total_complex_data_edges\tPrecision (% comp): ".($total_complex_data_edges/$total_data_edges)."\tRecall (% comp covered): ".($total_complex_data_edges/$total_possible_complex_edges)."\tNum prots: ".(scalar keys %data_prots)."\n";
#	print "\n";
		
	
}


sub ReadXValResultsIters {
	my $xvalresultsbasefilename = $_[0];
	my $xvalres_ref = $_[1];
	my $numcompedges_ref = $_[2];
	my $numnegedges_ref = $_[3];
	my $num_iters = $_[4];
	for (my $iter=0; $iter<$num_iters; $iter++) {
		open (XVALRESULTS_FILE, "$xvalresultsbasefilename xvalresults iter$iter.txt") || die "File $xvalresultsbasefilename xvalresults iter$iter.txt";
		my $line = <XVALRESULTS_FILE>;
		chomp $line;
		$$numcompedges_ref{$iter} = $line+0;
		$line = <XVALRESULTS_FILE>;
		chomp $line;
		$$numnegedges_ref{$iter} = $line+0;
		while ($line = <XVALRESULTS_FILE>) {
			chomp $line;
			(my $score, my $iscomplex) = split(/\t/, $line);
			my @tmparray = ($score+0, $iscomplex+0);
			push (@{$$xvalres_ref{$iter}}, \@tmparray);
		}
		close XVALRESULTS_FILE;
		print "Read xval results for iter $iter, results size = ".(scalar @{$$xvalres_ref{$iter}}).", num complex edges = $$numcompedges_ref{$iter}, num negative edges = $$numnegedges_ref{$iter}\n";
		
	}	
}

sub ReadXValResults($$$) {
	my $filename = $_[0];
	my $xvalresultsref = $_[1];
	my $numcomplexedgesref = $_[2];
	open (XVALRESULTS_FILE, $filename) || die $!;
	my $line = <XVALRESULTS_FILE>;
	chomp $line;
	if ($$numcomplexedgesref!=0 && $$numcomplexedgesref!=$line) { die; }
	$$numcomplexedgesref = $line + 0;
	while ($line = <XVALRESULTS_FILE>) {
		chomp $line;
		(my $score, my $iscomplex) = split(/\t/, $line);
		# use these when not enough memory to calculate prec-recall for all edges! Do both separately then combine results manually
		#if ($score < 0.00001) { next; }
#		if ($score >= 0.00001) { next; }
		my @tmparray = ($score+0, $iscomplex+0);
		push (@{$xvalresultsref}, \@tmparray);
	}
	close XVALRESULTS_FILE;
	print "Read xval results, results size = ".(scalar @{$xvalresultsref}).", num complex edges = $$numcomplexedgesref\n";
}


#--------------------- Filter data by approach -----------------------
# Should only filter continous data, NOT discretized data, because discretization sets an absent data score to 0 and then bins it.
# Thus after discretization the instance is considered to have the data score, so filtering by essential subtypes will not work.
sub FilterDataByApproach($$$) {
	my $dataref = $_[0];
	my $compsref = $_[2];
	my %edges_to_keep = ();
	my %data_to_use = (); 
	
	if (scalar @_ != 3) { die; }
	
	if ($_[1] eq "s" || $_[1] eq "sss" || $_[1] eq "sss_isoem") {
		$edges_to_keep{"PPIREL"} = 1;
		$edges_to_keep{"STRING"} = 1; 
		$edges_to_keep{"PUBMED"} = 1;
		$edges_to_keep{"PPIREL_DEG"} = 1;
		$edges_to_keep{"STRING_DEG"} = 1; 
		$edges_to_keep{"PUBMED_DEG"} = 1;
		$edges_to_keep{"PPIREL_NBC"} = 1;
		$edges_to_keep{"STRING_NBC"} = 1; 
		$edges_to_keep{"PUBMED_NBC"} = 1;
		$edges_to_keep{"PPIREL_CD"} = 1;
		$edges_to_keep{"STRING_CD"} = 1; 
		$edges_to_keep{"PUBMED_CD"} = 1;
		
		$data_to_use{"PPIREL"} = 1;
		$data_to_use{"STRING"} = 1; 
		$data_to_use{"PUBMED"} = 1;
		$data_to_use{"PPIREL_DEG"} = 1;
		$data_to_use{"STRING_DEG"} = 1; 
		$data_to_use{"PUBMED_DEG"} = 1;
		$data_to_use{"PPIREL_NBC"} = 1;
		$data_to_use{"STRING_NBC"} = 1; 
		$data_to_use{"PUBMED_NBC"} = 1;
		$data_to_use{"PPIREL_CD"} = 1;
		$data_to_use{"STRING_CD"} = 1; 
		$data_to_use{"PUBMED_CD"} = 1;
	}
	
	
	
	else {
		die "Error in FilterDataByApproach: approach $_[1]";
	}

	
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		if (defined $edges_to_keep{$subtype} && $$dataref{SUBTYPES}{$subtype}{TYPE} != "real") {
			die ("Error in FilterDataByApproach: keeping edges which are not real-type means all edges will be kept");
		}
	}
	
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		if (!defined $data_to_use{$subtype}) {
			delete $$dataref{SUBTYPES}{$subtype};
		}
	}
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		
		# check if this edge is to be kept
		my $keep = 1;
		if (scalar keys %edges_to_keep > 0) { 
			$keep = 0;
			foreach my $subtype (keys %edges_to_keep) {
				if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype} && $$dataref{INSTANCES}{$key}{DATA}{$subtype}!=0) {
					$keep = 1;
					last;
				}
			}
		}
		if ($keep==0) { 
			delete $$dataref{INSTANCES}{$key};
			next;
		}
		
		# delete data that is not used
		foreach my $subtype (keys %{$$dataref{INSTANCES}{$key}{DATA}}) {
			if (!defined $data_to_use{$subtype}) {
				delete $$dataref{INSTANCES}{$key}{DATA}{$subtype};
			}
		}
	}
	
	# calculate NUMPROTEINS and NUMCOMPLEXEDGES
	my %all_proteins = ();
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		(my $id1, my $id2) = split (/\|/, $key);
		$all_proteins{$id1} = 1;
		$all_proteins{$id2} = 1;
	}
	foreach my $comp (keys %{$compsref}) {
		foreach my $prot (keys %{$$compsref{$comp}}) {
			$all_proteins{$prot} = 1;
		}
	}
	my %all_compedges = ();
	foreach my $comp (keys %{$compsref}) {
		foreach my $prot1 (keys %{$$compsref{$comp}}) {
			foreach my $prot2 (keys %{$$compsref{$comp}}) {
				if ($prot2 le $prot1) { next; }
				$all_compedges{"$prot1|$prot2"} = 1;
			}
		}
	}
	$$dataref{NUMPROTEINS} = scalar keys %all_proteins;
	$$dataref{NUMCOMPLEXEDGES} = scalar keys %all_compedges;
					
				
	print "Filtered data\n";
	print "Subtypes:\n";
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		print "$subtype\t$$dataref{SUBTYPES}{$subtype}{TYPE}\n";
	}
	print "Num instances remaining: ".(scalar keys %{$$dataref{INSTANCES}})."\n\n";
	
}



#--------------------- print in ARFF format -----------------------
sub PrintArff($$) {
	
	open OUTPUTFILE, (">$_[0]");
	my $inputdataref = $_[1];
		
	# print header
	print OUTPUTFILE "\@relation \'complex_edge\'\n";
	foreach my $subtype (sort keys %{$$inputdataref{SUBTYPES}}) {
		print OUTPUTFILE "\@attribute \'$subtype\' ".($$inputdataref{SUBTYPES}{$subtype}{TYPE})."\n";
	}
	print OUTPUTFILE "\@attribute \'Class\' {\'complex\',\'noncomp\'}\n";
	print OUTPUTFILE "\@data\n";
	
	# print data
	foreach my $key (sort keys %{$$inputdataref{INSTANCES}}) {		
		foreach my $subtype (sort keys %{$$inputdataref{SUBTYPES}}) {
			if (defined $$inputdataref{INSTANCES}{$key}{DATA}{$subtype}) {
				print OUTPUTFILE ($$inputdataref{INSTANCES}{$key}{DATA}{$subtype}).", ";
			}
			else {
				print OUTPUTFILE "0, ";
			}
		}
		if ($$inputdataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
			print OUTPUTFILE "\'complex\'";
		}
		else {
			print OUTPUTFILE "\'noncomp\'";
		}
		print OUTPUTFILE "\n";
	}
	close OUTPUTFILE;

}

sub PrintTabDelimited($$) {
	
	open OUTPUTFILE, (">$_[0]");
	my $inputdataref = $_[1];
		
	# print header
	foreach my $subtype (sort keys %{$$inputdataref{SUBTYPES}}) {
		print OUTPUTFILE "$subtype\t";
	}
	print OUTPUTFILE "\n";
	
	# print data
	foreach my $key (sort keys %{$$inputdataref{INSTANCES}}) {		
		foreach my $subtype (sort keys %{$$inputdataref{SUBTYPES}}) {
			if (defined $$inputdataref{INSTANCES}{$key}{DATA}{$subtype}) {
				print OUTPUTFILE ($$inputdataref{INSTANCES}{$key}{DATA}{$subtype})."\t";
			}
			else {
				print OUTPUTFILE "0\t";
			}
		}
		if ($$inputdataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
			print OUTPUTFILE "complex";
		}
		else {
			print OUTPUTFILE "noncomp";
		}
		print OUTPUTFILE "\n";
	}
	close OUTPUTFILE;

}
# ------------------ Print likelihoods --------------
sub PrintLikelihoods($$$) {
	
	print "in PrintLikelihoods\n";
	
	# $$likelihoodsref{$iter}{$subtype}{$bin}{"COMPLEX"} = complex edges likelihood for this subtype and bin
	# $$likelihoodsref{$iter}{$subtype}{$bin}{"NONCOMP"} = noncomp edges likelihood for this subtype and bin
	# $$likelihoodsref{$iter}{$subtype}{$bin}{"ALL"} = all edges likelihood for this subtype and bin
	# $$likelihoodsref{$iter}{$subtype}{$bin}{"LIKERATIO"} = likelihood ratio of complex edge vs all edge
	# $$likelihoodsref{$iter}{"PRIORCOMP"} = prior prob of complex edge vs all edge
	my $likelihoodsref = $_[0];
	my $cutsref = $_[1]; # @{$$cutsref{$iter}{$subtype}} is an array of cutpoints for this subtype
	my $mode = $_[2];
	
	# Print likelihoods learned from ALL training data (no iterations)
	if ($mode eq "ALL") {
		foreach my $subtype (sort keys %{$$cutsref{"ALL"}}) {
			print "$subtype\n";
			print "score\tlikelihood\n";
			my $firstvalue = 0;
			my $lastvalue = 1;
			if ($subtype eq "Expression_" || $subtype eq "Expression_Rosetta") {
				$firstvalue = -1;
			}
			if ($subtype eq "COEXPR") {
				$lastvalue = 4;
			}
			if ($subtype eq "gravy") {
				$lastvalue = 3;
			}
			if (!defined $$likelihoodsref{"ALL"}{$subtype} || (scalar keys %{$$likelihoodsref{"ALL"}{$subtype}})!=(scalar @{$$cutsref{"ALL"}{$subtype}} + 1)) { 
				die "Error in PrintLikelihoods: feature cutpoints and bin number dont match";
			}
			my $cutpointsref = $$cutsref{"ALL"}{$subtype};
			for (my $cutidx=0; $cutidx < scalar @{$$cutsref{"ALL"}{$subtype}}; $cutidx++) {
				my $currscore = $$cutpointsref[$cutidx];				
				my $bin = 0;
				while ($bin < scalar @{$cutpointsref} && $currscore > ($$cutpointsref[$bin]+.0000000000001)) {
					$bin++;
				}
				my $likelihood = $$likelihoodsref{"ALL"}{$subtype}{$bin}{"LIKERATIO"};
				if ($cutidx==0) {
					print "$firstvalue\t$likelihood\n";
					print "$currscore\t$likelihood\n";
				}
				else {
					print (($$cutpointsref[$cutidx-1])."\t$likelihood\n");
					print "$currscore\t$likelihood\n";
				}
			}
			my $lastbin = scalar @{$$cutsref{"ALL"}{$subtype}};
			my $likelihood = $$likelihoodsref{"ALL"}{$subtype}{$lastbin}{"LIKERATIO"};
			print (($$cutsref{"ALL"}{$subtype}[$lastbin-1])."\t$likelihood\n");
			print "$lastvalue\t$likelihood\n";
			print "\n";
		} # end foreach subtype
		
		
	}
	
	# Print averaged likelihoods learned from iterations
	elsif ($mode eq "ITERS") {

		# $subtype_cutpoints{$subtype} is a sorted array of all cutpoints for this subtype
		my %subtype_cutpoints = ();
		foreach my $iter (keys %{$likelihoodsref}) {
			if ($iter eq "ALL") { next; }
			foreach my $subtype (keys %{$$likelihoodsref{$iter}}) {
				foreach my $cutpoint (@{$$cutsref{$iter}{$subtype}}) {
					$subtype_cutpoints{$subtype}{$cutpoint} = 1;
				}
			}
		}
		foreach my $subtype (keys %subtype_cutpoints) {
			my @tmp = sort {$a <=> $b} keys %{$subtype_cutpoints{$subtype}};
			$subtype_cutpoints{$subtype} = \@tmp;
		}
		
		
		foreach my $subtype (sort keys %subtype_cutpoints) {
			print "$subtype\n";
			print "score\tlikelihood\tSE\tlower\tupper\tN\n";
			my $firstvalue = 0;
			my $lastvalue = 1;
			if ($subtype eq "Expression_" || $subtype eq "Expression_Rosetta") {
				$firstvalue = -1;
			}
			if ($subtype eq "gravy") {
				$lastvalue = 3;
			}
#			if ($subtype eq "STRING") {
#				$lastvalue = 1000;
#			}
			for (my $cutidx=0; $cutidx < scalar @{$subtype_cutpoints{$subtype}}; $cutidx++) {
				my $currscore = $subtype_cutpoints{$subtype}[$cutidx];
				my @likvals;
				foreach my $iter (keys %{$likelihoodsref}) {
					if ($iter eq "ALL") { next; }
					if (!defined $$likelihoodsref{$iter}{$subtype}) { next; }
					my $cutpointsref = $$cutsref{$iter}{$subtype};
					my $bin = 0;
					while ($bin < scalar @{$cutpointsref} && $currscore > ($$cutpointsref[$bin]+.0000000000001)) {
						$bin++;
					}
					push (@likvals, $$likelihoodsref{$iter}{$subtype}{$bin}{"LIKERATIO"});
				}
				my $avglikelihood = 0;
				my $standarderr = 0;
				my $numsamples = scalar @likvals;
				foreach (@likvals) {
					$avglikelihood += $_;
				}
				if ($numsamples > 0) { 
					$avglikelihood /= $numsamples; 
				}
				if ($numsamples > 1) {
					foreach (@likvals) {
						$standarderr += ($_ - $avglikelihood)**2;
					}
					$standarderr /= ($numsamples-1);
					$standarderr /= $numsamples;
					$standarderr = $standarderr ** .5;
				}
				if ($cutidx==0) {
					print "$firstvalue\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n";
					print "$currscore\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n";
				}
				else {
					print (($subtype_cutpoints{$subtype}[$cutidx-1])."\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n");
					print "$currscore\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n";
				}
			}
			my @likvals;
			foreach my $iter (keys %{$likelihoodsref}) {
				if ($iter eq "ALL") { next; }
				if (!defined $$likelihoodsref{$iter}{$subtype}) { next; }
				my $lastbin = scalar @{$$cutsref{$iter}{$subtype}};
				push (@likvals, $$likelihoodsref{$iter}{$subtype}{$lastbin}{"LIKERATIO"});
			}
			my $avglikelihood = 0;
			my $standarderr = 0;
			my $numsamples = scalar @likvals;
			foreach (@likvals) {
				$avglikelihood += $_;
			}
			if ($numsamples > 0) {
				$avglikelihood /= $numsamples; 
			}
			if ($numsamples > 1) {
				foreach (@likvals) {
					$standarderr += ($_ - $avglikelihood)**2;
				}
				$standarderr /= ($numsamples-1);
				$standarderr /= $numsamples;
				$standarderr = $standarderr ** .5;
			}
			print (($subtype_cutpoints{$subtype}[scalar @{$subtype_cutpoints{$subtype}} - 1])."\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n");
			print "$lastvalue\t$avglikelihood\t$standarderr\t".($avglikelihood-$standarderr)."\t".($avglikelihood+$standarderr)."\t$numsamples\n";
			print "\n";
		} # end foreach subtype
	} # end mode eq "ITERS"
			
	else {
		die "Invalid mode in PrintLikelihoods";
	}
}

# ================== Isolatedness EM functions ============
sub CalcNeighbours {
	my $dataref = $_[0];
	my $neighboursref = $_[1];
	my $iso_thr = $_[2];
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} + $$dataref{INSTANCES}{$key}{SCORE}{4} < $iso_thr) { next; }
		(my $prot1, my $prot2) = split(/\|/, $key);
		$$neighboursref{$prot1}{$prot2} = 1;
		$$neighboursref{$prot2}{$prot1} = 1;
	}
	print "Calculated neighbours (threshold $iso_thr), size = ".(scalar keys %$neighboursref)."\n";
}

sub CalcTriangles {
	my $dataref = $_[0];
	my $trianglesref = $_[1];
	my $iso_thr = $_[2];
	
	my %nbs = ();
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} < $iso_thr) { next; }
		(my $prot1, my $prot2) = split(/\|/, $key);
		$nbs{$prot1}{$prot2} = 1;
		$nbs{$prot2}{$prot1} = 1;
	}
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} < $iso_thr) { next; }
		(my $prot1, my $prot2) = split(/\|/, $key);
		if (scalar keys %{$nbs{$prot1}} > scalar keys %{$nbs{$prot2}}) {
			my $tmpprot = $prot1;
			$prot1 = $prot2;
			$prot2 = $tmpprot;
		}
		foreach my $nb (keys %{$nbs{$prot1}}) {
			if ($nb eq $prot2) { next; }
			if (defined $nbs{$nb}{$prot2}) {
				my $triplet = join('|', sort ($prot1, $prot2, $nb));
				my $edge2 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				my $edge3 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$dataref{INSTANCES}{$edge2}{SCORE}{2} < $iso_thr) { die; }
				if ($$dataref{INSTANCES}{$edge3}{SCORE}{2} < $iso_thr) { die; }
				$$trianglesref{$triplet} = 1;
			}
		}
	}	
	print "Calculated triangles (threshold $iso_thr), size = ".(scalar keys %$trianglesref)."\n";
}

##--------------------- Calculate isolatedness, add it to currdata -----------------------
#sub CalcIso {
#	my $dataref = $_[0];
#	my $neighboursref = $_[1];
##	my $trianglesref = $_[2];
##	print "In CalcISO...\n";
#
#	$$dataref{SUBTYPES}{"ISO"}{TYPE} = "real";
#	foreach my $key (keys %{$$dataref{INSTANCES}}) {
##		print "$key, score $$dataref{INSTANCES}{$key}{SCORE}{2}\n";
#		$$dataref{INSTANCES}{$key}{DATA}{"ISO"} = $$dataref{INSTANCES}{$key}{SCORE}{2};   # use score for co-complex2
#		(my $prot1, my $prot2) = split(/\|/, $key);
#		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
#			if ($nb eq $prot2) { next; }
#			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
##			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
#			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
#				$$dataref{INSTANCES}{$key}{DATA}{"ISO"} *= (1 - .9999);
#			}
#			else {
#				$$dataref{INSTANCES}{$key}{DATA}{"ISO"} *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
#			}					
#		}
#		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
#			if ($nb eq $prot1) { next; }
#			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
##			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
#			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
#				$$dataref{INSTANCES}{$key}{DATA}{"ISO"} *= (1 - .9999);
#			}
#			else {
#				$$dataref{INSTANCES}{$key}{DATA}{"ISO"} *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
#			}					
#		}
##		print "Final score 		$$dataref{INSTANCES}{$key}{DATA}{ISO}\n";
#	}
#	
#}

# iso2(x,y) = score(x,y) * min(1-score(a,b)), where (a,b) <- all external edges
sub CalcIso2_min {
	my $dataref = $_[0];
	my $neighboursref = $_[1];
	my $iso_thr = $_[2];
#	print "In CalcISO...\n";

	$$dataref{SUBTYPES}{"ISO2"}{TYPE} = "real";
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
#		print "$key, score $$dataref{INSTANCES}{$key}{SCORE}{2}\n";
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} < $iso_thr) { 
			$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = 0;
			next; 
		}
		$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = $$dataref{INSTANCES}{$key}{SCORE}{2};   # use score for co-complex2
		(my $prot1, my $prot2) = split(/\|/, $key);
		my $highest_ext_score = 0;
		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
			if ($nb eq $prot2) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} > $highest_ext_score) {
				$highest_ext_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2};
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
			if ($nb eq $prot1) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} > $highest_ext_score) {
				$highest_ext_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2};
			}					
		}
		if ($highest_ext_score == 1) {
			$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - .9999);
		}
		else {
			$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - $highest_ext_score);
		}					
#		print "Final score 		$$dataref{INSTANCES}{$key}{DATA}{ISO2}\n";
	}
	
}



sub CalcIso2 {
	my $dataref = $_[0];
	my $neighboursref = $_[1];
	my $iso_thr = $_[2];
#	print "In CalcISO...\n";

	$$dataref{SUBTYPES}{"ISO2"}{TYPE} = "real";
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
#		print "$key, score $$dataref{INSTANCES}{$key}{SCORE}{2}\n";
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} < $iso_thr) { 
			$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = 0;
			next; 
		}
		$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = $$dataref{INSTANCES}{$key}{SCORE}{2};   # use score for co-complex2
		(my $prot1, my $prot2) = split(/\|/, $key);
		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
			if ($nb eq $prot2) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - .9999);
			}
			else {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
			if ($nb eq $prot1) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - .9999);
			}
			else {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
			}					
		}
#		print "Final score 		$$dataref{INSTANCES}{$key}{DATA}{ISO2}\n";
	}
	
}

sub CalcIso23combined {
	my $dataref = $_[0];
	my $neighboursref = $_[1];
	my $trianglesref = $_[2];
	my $iso_thr = $_[3];
	print "In CalcISO23combined...\n";

	# calculate ISO2 due to edge	
	$$dataref{SUBTYPES}{"ISO2"}{TYPE} = "real";
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
#		print "$key, score $$dataref{INSTANCES}{$key}{SCORE}{2}\n";
		if ($$dataref{INSTANCES}{$key}{SCORE}{2} < $iso_thr) { 
			$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = 0;
			next; 
		}
		$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} = $$dataref{INSTANCES}{$key}{SCORE}{2};   # use score for co-complex2
		(my $prot1, my $prot2) = split(/\|/, $key);
		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
			if ($nb eq $prot2) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2} + $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
			if ($ext_edge_score < $iso_thr) { die; }
#			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($ext_edge_score > .9999) {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - .9999);
			}
			else {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - $ext_edge_score);
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
			if ($nb eq $prot1) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2} + $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
			if ($ext_edge_score < $iso_thr) { die; }
#			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($ext_edge_score > .9999) {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - .9999);
			}
			else {
				$$dataref{INSTANCES}{$key}{DATA}{"ISO2"} *= (1 - $ext_edge_score);
			}					
		}
#		print "Final score 		$$dataref{INSTANCES}{$key}{DATA}{ISO2}\n";
	}
	
	# for each triangle, calculate the ISO3, add it to its edges' ISO2
	foreach my $triangle (keys %{$trianglesref}) {
		(my $prot1, my $prot2, my $prot3) = split(/\|/, $triangle);
		my $edge1 = "$prot1|$prot2";
		my $edge2 = "$prot1|$prot3";
		my $edge3 = "$prot2|$prot3";
#		if ($$dataref{INSTANCES}{$edge1}{DATA}{"PPI"}==0 || $$dataref{INSTANCES}{$edge2}{DATA}{"PPI"}==0 || $$dataref{INSTANCES}{$edge3}{DATA}{"PPI"}==0) { next; }
		if (!defined $$dataref{INSTANCES}{$edge1}{SCORE}{2} || !defined $$dataref{INSTANCES}{$edge2}{SCORE}{2} || !defined $$dataref{INSTANCES}{$edge3}{SCORE}{2}) { die; }
		if ($$dataref{INSTANCES}{$edge1}{SCORE}{2} < $iso_thr || 
				$$dataref{INSTANCES}{$edge2}{SCORE}{2} < $iso_thr || 
				$$dataref{INSTANCES}{$edge3}{SCORE}{2} < $iso_thr) { die; }
		my $tri_score = $$dataref{INSTANCES}{$edge1}{SCORE}{2} * $$dataref{INSTANCES}{$edge2}{SCORE}{2} * $$dataref{INSTANCES}{$edge3}{SCORE}{2};
		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
			if ($nb eq $prot2 || $nb eq $prot3) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}+$$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
			if ($ext_edge_score < $iso_thr) { die; }
#			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($ext_edge_score == .9999) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $ext_edge_score);
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
			if ($nb eq $prot1 || $nb eq $prot3) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}+$$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
			if ($ext_edge_score < $iso_thr) { die; }
#			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($ext_edge_score == .9999) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $ext_edge_score);
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot3}}) {
			if ($nb eq $prot1 || $nb eq $prot2) { next; }
			my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}+$$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
			if ($ext_edge_score < $iso_thr) { die; }
#			my $ext_edge_score = $$dataref{INSTANCES}{$ext_edge}{SCORE}{4};
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($ext_edge_score == .9999) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $ext_edge_score);
			}					
		}
		$$dataref{INSTANCES}{$edge1}{DATA}{"ISO2"} += $tri_score;
		$$dataref{INSTANCES}{$edge2}{DATA}{"ISO2"} += $tri_score;
		$$dataref{INSTANCES}{$edge3}{DATA}{"ISO2"} += $tri_score;
#		print "Final score $tri_score\n";
		
	}
	
	
	
}


sub CalcIso3_old {
	my $dataref = $_[0];
	my $neighboursref = $_[1];
	my $trianglesref = $_[2];
	my $iso_thr = $_[3];
	print "In CalcISO3...\n";

	$$dataref{SUBTYPES}{"ISO3"}{TYPE} = "real";
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		$$dataref{INSTANCES}{$key}{DATA}{"ISO3"} = 0;   
	}
	
	foreach my $triangle (keys %{$trianglesref}) {
		(my $prot1, my $prot2, my $prot3) = split(/\|/, $triangle);
		my $edge1 = "$prot1|$prot2";
		my $edge2 = "$prot1|$prot3";
		my $edge3 = "$prot2|$prot3";
#		if ($$dataref{INSTANCES}{$edge1}{DATA}{"PPI"}==0 || $$dataref{INSTANCES}{$edge2}{DATA}{"PPI"}==0 || $$dataref{INSTANCES}{$edge3}{DATA}{"PPI"}==0) { next; }
		if (!defined $$dataref{INSTANCES}{$edge1}{SCORE}{2} || !defined $$dataref{INSTANCES}{$edge2}{SCORE}{2} || !defined $$dataref{INSTANCES}{$edge3}{SCORE}{2}) { die; }
		if ($$dataref{INSTANCES}{$edge1}{SCORE}{2} < $iso_thr || $$dataref{INSTANCES}{$edge2}{SCORE}{2} < $iso_thr || $$dataref{INSTANCES}{$edge3}{SCORE}{2} < $iso_thr) { die; }
		my $tri_score = ($$dataref{INSTANCES}{$edge1}{SCORE}{2} + $$dataref{INSTANCES}{$edge2}{SCORE}{2} + $$dataref{INSTANCES}{$edge3}{SCORE}{2}) / 3;
		foreach my $nb (keys %{$$neighboursref{$prot1}}) {
			if ($nb eq $prot2 || $nb eq $prot3) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} < $iso_thr) { die; }
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot2}}) {
			if ($nb eq $prot1 || $nb eq $prot3) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} < $iso_thr) { die; }
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
			}					
		}
		foreach my $nb (keys %{$$neighboursref{$prot3}}) {
			if ($nb eq $prot1 || $nb eq $prot2) { next; }
			my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} < $iso_thr) { die; }
#			print "Ext edge $ext_edge, score $$dataref{INSTANCES}{$ext_edge}{SCORE}{2}\n";
			if ($$dataref{INSTANCES}{$ext_edge}{SCORE}{2} == 1) {
				$tri_score *= (1 - .9999);
			}
			else {
				$tri_score *= (1 - $$dataref{INSTANCES}{$ext_edge}{SCORE}{2});
			}					
		}
#		print "Final score $tri_score\n";
		if ($tri_score > $$dataref{INSTANCES}{$edge1}{DATA}{"ISO3"}) {
			$$dataref{INSTANCES}{$edge1}{DATA}{"ISO3"} = $tri_score;
		}
		if ($tri_score > $$dataref{INSTANCES}{$edge2}{DATA}{"ISO3"}) {
			$$dataref{INSTANCES}{$edge2}{DATA}{"ISO3"} = $tri_score;
		}
		if ($tri_score > $$dataref{INSTANCES}{$edge3}{DATA}{"ISO3"}) {
			$$dataref{INSTANCES}{$edge3}{DATA}{"ISO3"} = $tri_score;
		}
	}
	
}
# ================ Modelling (feature selection, discretization, parameter learning) functions ====================================
# ------------------------ Feature selection, only works on non-real data ----------------------
sub LearnFeatureSelection ($$$) {
	my $dataref = $_[0];
	my $selectedfeats = $_[1];
	my $approach = $_[2];
	my $total_possible_edges = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
	
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		
		if ($$dataref{SUBTYPES}{$subtype}{TYPE} != "real") { next; }
		
		# if any instances has no value for this subtype, set its value to 0
		foreach my $key (keys %{$$dataref{INSTANCES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$$dataref{INSTANCES}{$key}{DATA}{$subtype} = 0;
			}
		}
		my @instances;
		my @classes;
		foreach my $key (sort {$$dataref{INSTANCES}{$a}{DATA}{$subtype} <=> $$dataref{INSTANCES}{$b}{DATA}{$subtype}} keys %{$$dataref{INSTANCES}}) {
			if ($COUNT_ZERO_EDGES==2 && $$dataref{INSTANCES}{$key}{DATA}{$subtype}==0) { next; }
			push (@instances, $$dataref{INSTANCES}{$key}{DATA}{$subtype});
			if ($approach eq "sss") {
				push (@classes, $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS});
			}
			else {
				push (@classes, $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"});
			}
		}
		
		# if there are implicit edges, then ensure @instances has at least one 0 value!
		my $implicit_edges = 0;
		if ($COUNT_ZERO_EDGES==1 ||$COUNT_ZERO_EDGES==3) {
			$implicit_edges = $total_possible_edges - (scalar @instances);
			if ($implicit_edges > 0) {
				my $has_zero = 0;
				my $zero_position = -1; # if no 0, where to insert 0
				for (my $ii=0; $ii < scalar @instances; $ii++) {
					if ($instances[$ii]==0) {
						$has_zero = 1;
						last;
					}
					if ($instances[$ii] > 0) {
						$zero_position = $ii;
						last;
					}
				}
				if ($has_zero==0) {
					if ($zero_position==-1) {
						$zero_position = scalar @instances;
					}
					splice (@instances, $zero_position, 0, (0));
					splice (@classes, $zero_position, 0, (0));
					$implicit_edges--;
				}
			}
		}
		
		if (FeatureSelectionMDL(\@instances, \@classes, 0, scalar @instances - 1, $implicit_edges) == 1) {
			$$selectedfeats{$subtype} = 1;
		}
	}

	print "Selected features:\n";
	foreach my $subtype (sort keys %{$selectedfeats}) {
		print "$subtype\n";
	}
	print "\n";

}

sub SelectFeatures ($$$$) {
	my $dataref = $_[0];
	my $selectedfeats = $_[1];	
	my $approach = $_[2];
	my $compsref = $_[3];
	
	
	my %edges_to_keep = ();
	if ($approach eq "sss" || $approach eq "s") {
		$edges_to_keep{"PPIREL"} = 1;
		$edges_to_keep{"STRING"} = 1; 
		$edges_to_keep{"PUBMED"} = 1;
		$edges_to_keep{"PPIREL_DEG"} = 1;
		$edges_to_keep{"STRING_DEG"} = 1; 
		$edges_to_keep{"PUBMED_DEG"} = 1;
		$edges_to_keep{"PPIREL_NBC"} = 1;
		$edges_to_keep{"STRING_NBC"} = 1; 
		$edges_to_keep{"PUBMED_NBC"} = 1;
		$edges_to_keep{"PPIREL_CD"} = 1;
		$edges_to_keep{"STRING_CD"} = 1; 
		$edges_to_keep{"PUBMED_CD"} = 1;
	}	
	
	else {
		die "Error in SelectFeatures: approach $approach";
	}
	
	# delete the subtype definitions
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		if (!defined $$selectedfeats{$subtype}) {
			delete $$dataref{SUBTYPES}{$subtype};
		}
	}
	foreach my $subtype (keys %edges_to_keep) {
		if (!defined $$selectedfeats{$subtype}) {
			delete $edges_to_keep{$subtype};
		}
	}
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		# delete the subtype from each edge
		foreach my $subtype (keys %{$$dataref{INSTANCES}{$key}{DATA}}) {
			if (!defined $$selectedfeats{$subtype}) {
				delete $$dataref{INSTANCES}{$key}{DATA}{$subtype};
			}
		}
		
		# check if this edge is to be kept
		my $keep = 0;
		foreach my $subtype (keys %edges_to_keep) {
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype} && $$dataref{INSTANCES}{$key}{DATA}{$subtype}!=0) {
				$keep = 1;
				last;
			}
		}
		if ($keep==0) { 
			delete $$dataref{INSTANCES}{$key};
		}
	}
	
	# calculate NUMPROTEINS and NUMCOMPLEXEDGES
	my %all_proteins = ();
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		(my $id1, my $id2) = split (/\|/, $key);
		$all_proteins{$id1} = 1;
		$all_proteins{$id2} = 1;
	}
	foreach my $comp (keys %{$compsref}) {
		foreach my $prot (keys %{$$compsref{$comp}}) {
			$all_proteins{$prot} = 1;
		}
	}
	my %all_compedges = ();
	foreach my $comp (keys %{$compsref}) {
		foreach my $prot1 (keys %{$$compsref{$comp}}) {
			foreach my $prot2 (keys %{$$compsref{$comp}}) {
				if ($prot2 le $prot1) { next; }
				$all_compedges{"$prot1|$prot2"} = 1;
			}
		}
	}
	$$dataref{NUMPROTEINS} = scalar keys %all_proteins;
	$$dataref{NUMCOMPLEXEDGES} = scalar keys %all_compedges;
				
}

	
# ------------------------ Discretization -------------------------------------------------
sub LearnDiscretizationRanges ($$$) {
	my $dataref = $_[0];
	my $allcutpoints = $_[1];
	my $approach = $_[2];
	my $disc_type = $_[3];
	my $total_possible_edges = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
	
	foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
		if (defined $disc_type && $subtype ne $disc_type) { next; }
		if ($$dataref{SUBTYPES}{$subtype}{TYPE} != "real") { next; }
		
		# if any instances has no value for this subtype, set its value to 0
		foreach my $key (keys %{$$dataref{INSTANCES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$$dataref{INSTANCES}{$key}{DATA}{$subtype} = 0;
			}
		}
		my @instances;
		my @classes;
		foreach my $key (sort {$$dataref{INSTANCES}{$a}{DATA}{$subtype} <=> $$dataref{INSTANCES}{$b}{DATA}{$subtype}} keys %{$$dataref{INSTANCES}}) {
			if ($COUNT_ZERO_EDGES==2 && $$dataref{INSTANCES}{$key}{DATA}{$subtype}==0) { next; }
			push (@instances, $$dataref{INSTANCES}{$key}{DATA}{$subtype});
			if ($approach eq "sss") {
				push (@classes, $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS});
			}
			else {
				push (@classes, $$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"});
			}
		}
#		print "subtype $subtype, num instances ".(scalar @instances).", num classes ".(scalar @classes)."\n";
				
		# if there are implicit edges, then ensure @instances has at least one 0 value!
		my $implicit_edges = 0;
		if ($COUNT_ZERO_EDGES==1 || $COUNT_ZERO_EDGES==3) {
			$implicit_edges = $total_possible_edges - (scalar @instances);
			if ($implicit_edges > 0) {
				my $has_zero = 0;
				my $zero_position = -1; # if no 0, where to insert 0
				for (my $ii=0; $ii < scalar @instances; $ii++) {
					if ($instances[$ii]==0) {
						$has_zero = 1;
						last;
					}
					if ($instances[$ii] > 0) {
						$zero_position = $ii;
						last;
					}
				}
				if ($has_zero==0) {
					if ($zero_position==-1) {
						$zero_position = scalar @instances;
					}
					splice (@instances, $zero_position, 0, (0));
					splice (@classes, $zero_position, 0, (0));
					$implicit_edges--;
				}
			}
		}
		
		my $cutpoints = DiscretizeMDL(\@instances, \@classes, 0, scalar @instances - 1, $implicit_edges);
		$$allcutpoints{$subtype} = $cutpoints;
		
		# if want to ignore all zeros, create a new cutpoint between 0 and the lowest value
		if ($COUNT_ZERO_EDGES==2) {
			
		}
	}
	
	print "Discretization cutpoints:\n";
	foreach my $subtype (sort keys %{$allcutpoints}) {
		print "$subtype cutpoints:\n";
		foreach (@{$$allcutpoints{$subtype}}) {
			print "$_\n";
		}
	}
	print "\n";
}


sub DiscretizeData ($$) {
	my $dataref = $_[0];
	my $allcutpoints = $_[1];
	my $disc_type = $_[2];
		
	foreach my $subtype (keys %{$allcutpoints}) {
		if (defined $disc_type && $subtype ne $disc_type) { next; }
		if ($$dataref{SUBTYPES}{$subtype}{TYPE} != "real") { next; }
		
		# if any instances has no value for this subtype, set its value to 0
		foreach my $key (keys %{$$dataref{INSTANCES}}) {
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype} && $subtype eq "gravy") {
				die ("In LearnFeatureSelection: gravy subtype with no score!");
			}				
			if (!defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$$dataref{INSTANCES}{$key}{DATA}{$subtype} = 0;
			}
		}
		
		my $cutpoints = $$allcutpoints{$subtype};
		my @subtypetype;
		for (my $count=0; $count<=scalar @$cutpoints; $count++) {
			push @subtypetype, $count;
		}
		$$dataref{SUBTYPES}{$subtype}{TYPE} = "{".(join(", ", @subtypetype))."}";	
			
		foreach my $key (keys %{$$dataref{INSTANCES}}) {
			my $score = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$score = $$dataref{INSTANCES}{$key}{DATA}{$subtype};
			}
			my $bin = 0;
			while ($bin < scalar @$cutpoints && $score > $$cutpoints[$bin]) {
					$bin++;
			}
			$$dataref{INSTANCES}{$key}{DATA}{$subtype} = $bin;
		}
	}	
}

# instances: reference to array of instances, sorted in ascending order
# classes: reference to array of classes, corresponding to instances
# start_idx, end_idx: inclusive start and end indices to look at
sub DiscretizeMDL () {
	if (scalar @_ != 5 && scalar @_ != 7) { die; }
	
	my $instances = $_[0];
	my $classes = $_[1];
	my $start_idx = $_[2];
	my $end_idx = $_[3];
	my $implicit_edges = $_[4]; 

	if ($implicit_edges > 0 && ($$instances[$start_idx]>0 || $$instances[$end_idx]<0)) {
		die "Error in DiscretizeMDL: implicit_edges > 0, but range does not contain 0!";
	}
	
#	print "----- Discretize MDL -----\n";
#	print "start idx = $start_idx, end idx = $end_idx\n";

	my $classes_count;
	my $instance_value_classes;
	if (scalar @_ ==7) {
		# classes_count: reference to hash with number of classes in the instances in the range of interest
		# instance_value_classes: the (possibly multiple) class values at each instance value. Used for finding boundary points
		$classes_count = $_[5];
		$instance_value_classes = $_[6];
	}	
	
	# count number of classes
	if (scalar keys %{$classes_count} == 0) {
#		print "calculating classes_count\n";
		for (my $idx=$start_idx; $idx<=$end_idx; $idx++) {
			$$classes_count{$$classes[$idx]}++;
		}
		$$classes_count{0} += $implicit_edges;
	}
	my $num_classes = 0;
	foreach my $class (keys %{$classes_count}) {
		if ($$classes_count{$class} > 0) {
			$num_classes++;
		}
	}
#	print "num classes = $num_classes\n";
	if ($num_classes==1 || $$instances[$start_idx]==$$instances[$end_idx]) {
		my @allcuts = (-1);
		return \@allcuts;
	}	
	my $num_total = $end_idx - $start_idx + 1 + $implicit_edges;
#	print "classes count:\n";
#	foreach my $class (sort keys %{$classes_count}) {
#		print "class $class\t$$classes_count{$class}\n";
#	}
	if (scalar keys %{$instance_value_classes} == 0) {
#		print "calculating instance_value_classes\n";
		for (my $idx=$start_idx; $idx <= $end_idx; $idx++) {
			$$instance_value_classes{$$instances[$idx]}{$$classes[$idx]} = 1;
		}
		if ($implicit_edges>0) {
			$$instance_value_classes{0}{0} = 1;
		}
	}
#	foreach my $inst (sort keys %{$instance_value_classes}) {
#		foreach my $class (sort keys %{$$instance_value_classes{$inst}}) {
#			print "instance $inst class $class\n";
#		}
#	}

	# calculate entropy of entire set
	my $s_entropy = 0;
	foreach my $class (keys %{$classes_count}) {
		if ($$classes_count{$class}==0) { next; }
		$s_entropy -= $$classes_count{$class}/$num_total*log($$classes_count{$class}/$num_total)/LOG2;
	}
#	print "set entropy = $s_entropy\n";

	# iterate through instances in order of its subtype value
	my %s1_classes_count = ();
	my %s2_classes_count = ();
	foreach my $class (keys %{$classes_count}) {
		$s1_classes_count{$class} = 0;
		$s2_classes_count{$class} = $$classes_count{$class};
	}
	my $best_partition_entropy = 99999999;
	my $best_partition_k1; # number of classes in subset 1
	my $best_partition_k2; # number of classes in subset 2
	my $best_partition_s1_entropy; # entropy in subset 1
	my $best_partition_s2_entropy; # entropy in subset 2
	my $best_partition_index;
	my %best_partition_s1_classes_count;
	my %best_partition_s2_classes_count;
	for (my $idx=$start_idx; $idx<$end_idx; $idx++) {
			
		$s1_classes_count{$$classes[$idx]}++;
		$s2_classes_count{$$classes[$idx]}--;
		
		if ($implicit_edges>0 && $$instances[$idx]<=0 && $$instances[$idx+1]>0) {
			$s1_classes_count{0} += $implicit_edges;
			$s2_classes_count{0} -= $implicit_edges;
		}
#		print "\nidx = $idx\n";
#		print "s1 classes count:\n";
#		foreach my $class (sort keys %s1_classes_count) {
#			print "class $class\t$s1_classes_count{$class}\n";
#		}
#		print "s2 classes count:\n";
#		foreach my $class (sort keys %s2_classes_count) {
#			print "class $class\t$s2_classes_count{$class}\n";
#		}
			
		
		# don't look at non-boundary-points
		if ($$instances[$idx]==$$instances[$idx+1]) { next; }
		my $is_boundary_point = 0;
		foreach my $c1 (keys %{$$instance_value_classes{$$instances[$idx]}}) {
			foreach my $c2 (keys %{$$instance_value_classes{$$instances[$idx+1]}}) {
#				print "instance $$instances[$idx] has class $c1, instance ".$$instances[$idx+1]." has class $c2\n";
				if ($c1 != $c2) { 
					$is_boundary_point = 1;
					last;
				}
			}
			if ($is_boundary_point==1) { last; }
		}
		if ($is_boundary_point==0) { next; }
#		print "\nIndex $idx, instance $$instances[$idx], is boundary point\n";
			
		# calculate the partition class information entropy E(A, T; S)
		my $s1_size = $idx - $start_idx + 1;
		if ($$instances[$idx]>0 || ($$instances[$idx]<=0 && $$instances[$idx+1]>0)) {
			$s1_size += $implicit_edges;
		}
		my $s2_size = $num_total - $s1_size;
		my $s1_entropy = 0;
		foreach my $class (keys %s1_classes_count) {
			if ($s1_classes_count{$class}==0) { next; }
			$s1_entropy -= $s1_classes_count{$class}/$s1_size*log($s1_classes_count{$class}/$s1_size)/LOG2;
		}
		my $s2_entropy = 0;
		foreach my $class (keys %s2_classes_count) {
			if ($s2_classes_count{$class}==0) { next; }
			$s2_entropy -= $s2_classes_count{$class}/$s2_size*log($s2_classes_count{$class}/$s2_size)/LOG2;
		}
		my $partition_entropy = $s1_size/$num_total * $s1_entropy + $s2_size/$num_total * $s2_entropy;
#		print "s1 size = $s1_size, s2 size = $s2_size, s1_entropy = $s1_entropy, s2_entropy = $s2_entropy, partition entropy = $partition_entropy\n";
		
		# save the best partition so far
		if ($partition_entropy < $best_partition_entropy) {
			$best_partition_entropy = $partition_entropy;
			$best_partition_k1 = 0;
			foreach my $class (keys %s1_classes_count) {
				if ($s1_classes_count{$class}>0) { $best_partition_k1++; }
			}
			$best_partition_k2 = 0;
			foreach my $class (keys %s2_classes_count) {
				if ($s2_classes_count{$class}>0) { $best_partition_k2++; }
			}
			$best_partition_s1_entropy = $s1_entropy;
			$best_partition_s2_entropy = $s2_entropy;
			$best_partition_index = $idx;
			%best_partition_s1_classes_count = %s1_classes_count;
			%best_partition_s2_classes_count = %s2_classes_count;
#			print "Best partition so far!\n";
#			print "Partition entropy = $best_partition_entropy, k1 = $best_partition_k1, k2 = $best_partition_k2, s1_entropy = $best_partition_s1_entropy, s2_entropy = $best_partition_s2_entropy, index = $best_partition_index\n";
#			print "best partition s1 classes count:\n";
#			foreach my $class (sort keys %best_partition_s1_classes_count) {
#				print "class $class $best_partition_s1_classes_count{$class}\n";
#			}
#			print "best partition s2 classes count:\n";
#			foreach my $class (sort keys %best_partition_s2_classes_count) {
#				print "class $class $best_partition_s2_classes_count{$class}\n";
#			}
		}
	}
		
	# decide to cut or not
	my $gain = $s_entropy - $best_partition_entropy;
	my $delta = log(3**$num_classes-2)/LOG2 - ($num_classes*$s_entropy - $best_partition_k1*$best_partition_s1_entropy - $best_partition_k2*$best_partition_s2_entropy);
	
#	print "cut at ".(($$instances[$best_partition_index]+$$instances[$best_partition_index+1])/2).", index = $best_partition_index? gain = $gain, delta = $delta, RHS = ".(log($num_total-1)/(LOG2*$num_total) + $delta/$num_total)."\n";
	
	if ($gain > $MDL_ALPHA * (log($num_total-1)/(LOG2*$num_total) + $delta/$num_total)) {		
#		print "cut\n";
		my $s1_cuts = DiscretizeMDL($instances, $classes, $start_idx, $best_partition_index, ($$instances[$start_idx]<=0 && $$instances[$best_partition_index]>=0)?$implicit_edges:0, \%best_partition_s1_classes_count, $instance_value_classes);
		my $s2_cuts = DiscretizeMDL($instances, $classes, $best_partition_index+1, $end_idx, ($$instances[$best_partition_index+1]<=0 && $$instances[$end_idx]>=0)?$implicit_edges:0, \%best_partition_s2_classes_count, $instance_value_classes);
		my @allcuts;
		if ($$s1_cuts[0] != -1) {
			@allcuts = (@allcuts, @$s1_cuts);
		}
		@allcuts = (@allcuts, ($$instances[$best_partition_index]+$$instances[$best_partition_index+1])/2);
		if ($$s2_cuts[0] != -1) {
			@allcuts = (@allcuts, @$s2_cuts);
		}
		return \@allcuts;
	}
	else {		
#		print "no cut\n";		
		my @allcuts = (-1);
		return \@allcuts;
	}
}




# instances: reference to array of instances, sorted in ascending order
# classes: reference to array of classes, corresponding to instances
# start_idx, end_idx: inclusive start and end indices to look at
# implicit_edges: the number of implicit edges (not represented in instances because there are too many!). These have instance and class values = 0
sub FeatureSelectionMDL () {
	my $instances = $_[0];
	my $classes = $_[1];
	my $start_idx = $_[2];
	my $end_idx = $_[3];
	my $implicit_edges = $_[4];
	
	if (scalar @_ != 5) { die; }
	
	if ($implicit_edges > 0 && ($$instances[$start_idx]>0 || $$instances[$end_idx]<0)) {
		die "Error in FeatureSelectionMDL: implicit_edges > 0, but range does not contain 0!";
	}

#	print "----- Discretize MDL -----\n";
#	print "start idx = $start_idx, end idx = $end_idx\n";

	my $classes_count;
	my $instance_value_classes;
	for (my $idx=$start_idx; $idx<=$end_idx; $idx++) {
		$$classes_count{$$classes[$idx]}++;
	}
	if ($implicit_edges > 0) {
		$$classes_count{0} += $implicit_edges;
	}
	my $num_classes = 0;
	foreach my $class (keys %{$classes_count}) {
		if ($$classes_count{$class} > 0) {
			$num_classes++;
		}
	}
#	print "num classes = $num_classes\n";
	if ($num_classes==1 || $$instances[$start_idx]==$$instances[$end_idx]) {
		return 0;
	}	
	my $num_total = $end_idx - $start_idx + 1 + $implicit_edges;	
#	print "classes count:\n";
#	foreach my $class (sort keys %{$classes_count}) {
#		print "class $class\t$$classes_count{$class}\n";
#	}
	for (my $idx=$start_idx; $idx <= $end_idx; $idx++) {
		$$instance_value_classes{$$instances[$idx]}{$$classes[$idx]} = 1;
	}
	if ($implicit_edges > 0) {
		$$instance_value_classes{0}{0} = 1;
	}
#	foreach my $inst (sort keys %{$instance_value_classes}) {
#		foreach my $class (sort keys %{$$instance_value_classes{$inst}}) {
#			print "instance $inst class $class\n";
#		}
#	}

	# calculate entropy of entire set
	my $s_entropy = 0;
	foreach my $class (keys %{$classes_count}) {
		if ($$classes_count{$class}==0) { next; }
		$s_entropy -= $$classes_count{$class}/$num_total*log($$classes_count{$class}/$num_total)/LOG2;
	}
#	print "set entropy = $s_entropy\n";



	# iterate through instances in order of its subtype value
	my %s1_classes_count = ();
	my %s2_classes_count = ();
	foreach my $class (keys %{$classes_count}) {
		$s1_classes_count{$class} = 0;
		$s2_classes_count{$class} = $$classes_count{$class};
	}
	my $best_partition_entropy = 99999999;
	my $best_partition_k1; # number of classes in subset 1
	my $best_partition_k2; # number of classes in subset 2
	my $best_partition_s1_entropy; # entropy in subset 1
	my $best_partition_s2_entropy; # entropy in subset 2
	my $best_partition_index;
	my %best_partition_s1_classes_count;
	my %best_partition_s2_classes_count;
	for (my $idx=$start_idx; $idx<$end_idx; $idx++) {
			
		$s1_classes_count{$$classes[$idx]}++;
		$s2_classes_count{$$classes[$idx]}--;
		
		if ($implicit_edges>0 && $$instances[$idx]<=0 && $$instances[$idx+1]>0) {
			$s1_classes_count{0} += $implicit_edges;
			$s2_classes_count{0} -= $implicit_edges;
		}
#		print "\nidx = $idx\n";
#		print "s1 classes count:\n";
#		foreach my $class (sort keys %s1_classes_count) {
#			print "class $class\t$s1_classes_count{$class}\n";
#		}
#		print "s2 classes count:\n";
#		foreach my $class (sort keys %s2_classes_count) {
#			print "class $class\t$s2_classes_count{$class}\n";
#		}
			
		
		# don't look at non-boundary-points
		if ($$instances[$idx]==$$instances[$idx+1]) { next; }
		my $is_boundary_point = 0;
		foreach my $c1 (keys %{$$instance_value_classes{$$instances[$idx]}}) {
			foreach my $c2 (keys %{$$instance_value_classes{$$instances[$idx+1]}}) {
#				print "instance $$instances[$idx] has class $c1, instance ".$$instances[$idx+1]." has class $c2\n";
				if ($c1 != $c2) { 
					$is_boundary_point = 1;
					last;
				}
			}
			if ($is_boundary_point==1) { last; }
		}
		if ($is_boundary_point==0) { next; }
#		print "\nIndex $idx, instance $$instances[$idx], is boundary point\n";
			
		# calculate the partition class information entropy E(A, T; S)
		my $s1_size = $idx - $start_idx + 1;
		if ($$instances[$idx]>0 || ($$instances[$idx]<=0 && $$instances[$idx+1]>0)) {
			$s1_size += $implicit_edges;
		}
		my $s2_size = $num_total - $s1_size;
		my $s1_entropy = 0;
		foreach my $class (keys %s1_classes_count) {
			if ($s1_classes_count{$class}==0) { next; }
			$s1_entropy -= $s1_classes_count{$class}/$s1_size*log($s1_classes_count{$class}/$s1_size)/LOG2;
		}
		my $s2_entropy = 0;
		foreach my $class (keys %s2_classes_count) {
			if ($s2_classes_count{$class}==0) { next; }
			$s2_entropy -= $s2_classes_count{$class}/$s2_size*log($s2_classes_count{$class}/$s2_size)/LOG2;
		}
		my $partition_entropy = $s1_size/$num_total * $s1_entropy + $s2_size/$num_total * $s2_entropy;
#		print "s1 size = $s1_size, s2 size = $s2_size, s1_entropy = $s1_entropy, s2_entropy = $s2_entropy, partition entropy = $partition_entropy\n";
		
		# save the best partition so far
		if ($partition_entropy < $best_partition_entropy) {
			$best_partition_entropy = $partition_entropy;
			$best_partition_k1 = 0;
			foreach my $class (keys %s1_classes_count) {
				if ($s1_classes_count{$class}>0) { $best_partition_k1++; }
			}
			$best_partition_k2 = 0;
			foreach my $class (keys %s2_classes_count) {
				if ($s2_classes_count{$class}>0) { $best_partition_k2++; }
			}
			$best_partition_s1_entropy = $s1_entropy;
			$best_partition_s2_entropy = $s2_entropy;
			$best_partition_index = $idx;
			%best_partition_s1_classes_count = %s1_classes_count;
			%best_partition_s2_classes_count = %s2_classes_count;
#			print "Best partition so far!\n";
#			print "Partition entropy = $best_partition_entropy, k1 = $best_partition_k1, k2 = $best_partition_k2, s1_entropy = $best_partition_s1_entropy, s2_entropy = $best_partition_s2_entropy, index = $best_partition_index\n";
#			print "best partition s1 classes count:\n";
#			foreach my $class (sort keys %best_partition_s1_classes_count) {
#				print "class $class $best_partition_s1_classes_count{$class}\n";
#			}
#			print "best partition s2 classes count:\n";
#			foreach my $class (sort keys %best_partition_s2_classes_count) {
#				print "class $class $best_partition_s2_classes_count{$class}\n";
#			}
		}
	}
		
	# decide to cut or not
	my $gain = $s_entropy - $best_partition_entropy;
	my $delta = log(3**$num_classes-2)/LOG2 - ($num_classes*$s_entropy - $best_partition_k1*$best_partition_s1_entropy - $best_partition_k2*$best_partition_s2_entropy);
	
#	print "cut at ".(($$instances[$best_partition_index]+$$instances[$best_partition_index+1])/2).", index = $best_partition_index? gain = $gain, delta = $delta, RHS = ".(log($num_total-1)/(LOG2*$num_total) + $delta/$num_total)."\n";
	
	if ($gain > $MDL_ALPHA * (log($num_total-1)/(LOG2*$num_total) + $delta/$num_total)) {
		return 1;
	}
	else {		
		return 0;
	}
}



#--------------------- Calculate likelihoods -----------------------
sub CalcLikelihoods($$$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	
	# $counts{$subtype}{$bin}{"COMPLEX"} = # samples with that subtype & bin value and are co-complex edges
	# $counts{$subtype}{$bin}{"NONCOMP"} = # samples with that subtype & bin value and are non-co-complex edges
	# $counts{$subtype}{$bin}{"NONCOMP_BORDER"} = # samples with that subtype & bin value and are non-co-complex edges, but at least one interactor is in a complex
	# $counts{$subtype}{$bin}{"NONCOMP_NONBORDER"} = # samples with that subtype & bin value and are non-co-complex edges, and no interactor is in a complex
	# $counts{$subtype}{$bin}{"ALL"} = # samples with that subtype & bin value
	# Bin 0 also represents all implicit edges
	my %counts = ();
	my $complex_edges = 0;
	my $noncomp_edges = 0;
	my $noncomp_border_edges = 0;	
	my $noncomp_nonborder_edges = 0;
	my $all_edges = 0;
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		# counts for "COMPLEX" model
		if ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
			$complex_edges++;
			foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
				my $bin = 0;
				if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
					$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
				}
				$counts{$subtype}{$bin}{"COMPLEX"}++;
			}
		}
		# counts for "NONCOMP" model
		else {
			$noncomp_edges++;
			foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
				my $bin = 0;
				if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
					$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
				}
				$counts{$subtype}{$bin}{"NONCOMP"}++;
				
			}
			
		}
		# counts for "ALL" model (complex and noncomp)
		$all_edges++;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$counts{$subtype}{$bin}{"ALL"}++;
		}
	}
	
	if ($COUNT_ZERO_EDGES==1 || $COUNT_ZERO_EDGES==3) {
		# set counts for bin 0 to include implicit edges
		foreach my $subtype (keys %counts) {
			$counts{$subtype}{0}{"COMPLEX"} = $$dataref{NUMCOMPLEXEDGES};
#			$counts{$subtype}{0}{"COMPLEX"} = 11295;
			$counts{$subtype}{0}{"ALL"} = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
			$counts{$subtype}{0}{"NONCOMP"} = $counts{$subtype}{0}{"ALL"} - $counts{$subtype}{0}{"COMPLEX"};
			foreach my $bin (keys %{$counts{$subtype}}) {
				if ($bin==0) { next; }
				$counts{$subtype}{0}{"COMPLEX"} -= $counts{$subtype}{$bin}{"COMPLEX"};
				$counts{$subtype}{0}{"ALL"} -= $counts{$subtype}{$bin}{"ALL"};
				$counts{$subtype}{0}{"NONCOMP"} -= $counts{$subtype}{$bin}{"NONCOMP"};
			}
		}
		$all_edges = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
		$complex_edges = $$dataref{NUMCOMPLEXEDGES};
		#$complex_edges = 11295;
		$noncomp_edges = $all_edges - $complex_edges;
	}
	
	# check that the counts sum up correctly
	foreach my $subtype (keys %counts) {
		my $comp_edges_check = 0;
		my $noncomp_edges_check = 0;
		my $all_edges_check = 0;
		my $noncomp_border_edges_check = 0;
		my $noncomp_nonborder_edges_check = 0;
		foreach my $bin (keys %{$counts{$subtype}}) {
			if ($counts{$subtype}{$bin}{"ALL"} != $counts{$subtype}{$bin}{"COMPLEX"} + $counts{$subtype}{$bin}{"NONCOMP"}) {
				die; 
			}	
			$comp_edges_check += $counts{$subtype}{$bin}{"COMPLEX"};
			$noncomp_edges_check += $counts{$subtype}{$bin}{"NONCOMP"};
			$all_edges_check += $counts{$subtype}{$bin}{"ALL"};
		}
		if ($comp_edges_check != $complex_edges) { die; }
		if ($noncomp_edges_check != $noncomp_edges) { die; }
		if ($all_edges_check != $all_edges) { die; }
	}
			
	
	# laplace smoothing
	foreach my $subtype (keys %counts) {
		foreach my $bin (keys %{$counts{$subtype}}) {
			$counts{$subtype}{$bin}{"COMPLEX"}++;
			$counts{$subtype}{$bin}{"NONCOMP"}++;
			$counts{$subtype}{$bin}{"ALL"}+=2;
			$complex_edges++;
			$noncomp_edges++;
			$all_edges+=2;
		}
	}
	
#	# laplace smoothing (wrong)
#	$complex_edges++;
#	$noncomp_edges++;
#	$all_edges++;
#	foreach my $subtype (keys %counts) {
#		foreach my $bin (keys %{$counts{$subtype}}) {
#			$counts{$subtype}{$bin}{"COMPLEX"}++;
#			$counts{$subtype}{$bin}{"NONCOMP"}++;
#			$counts{$subtype}{$bin}{"ALL"}++;
#		}
#	}
		
	# calculate likelihoods
	foreach my $subtype (keys %counts) {
		foreach my $bin (keys %{$counts{$subtype}}) {
			if ($CUMULATIVE_PROB_LIKELIHOODS==0) {
				$$likelihoodsref{$subtype}{$bin}{"COMPLEX"} = $counts{$subtype}{$bin}{"COMPLEX"} / $complex_edges;
				$$likelihoodsref{$subtype}{$bin}{"NONCOMP"} = $counts{$subtype}{$bin}{"NONCOMP"} / $noncomp_edges;
			}
			else {
				my $cumul_comp_count = 0;
				my $cumul_noncomp_count = 0;
				for (my $binidx=$bin; $binidx < scalar keys %{$counts{$subtype}}; $binidx++) {
					$cumul_comp_count += $counts{$subtype}{$binidx}{"COMPLEX"};
					$cumul_noncomp_count += $counts{$subtype}{$binidx}{"NONCOMP"};
				}
				$$likelihoodsref{$subtype}{$bin}{"COMPLEX"} = $cumul_comp_count / $complex_edges;
				$$likelihoodsref{$subtype}{$bin}{"NONCOMP"} = $cumul_noncomp_count / $noncomp_edges;
			}
			$$likelihoodsref{$subtype}{$bin}{"ALL"} = $counts{$subtype}{$bin}{"ALL"} / $all_edges;
			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"} = $$likelihoodsref{$subtype}{$bin}{"COMPLEX"} / $$likelihoodsref{$subtype}{$bin}{"NONCOMP"};
		}
		if ($COUNT_ZERO_EDGES==3) {
			$$likelihoodsref{$subtype}{0}{"COMPLEX"} = 1;
			$$likelihoodsref{$subtype}{0}{"NONCOMP"} = 1;
			$$likelihoodsref{$subtype}{0}{"ALL"} = 1;
			$$likelihoodsref{$subtype}{0}{"LIKERATIO"} = $$likelihoodsref{$subtype}{0}{"COMPLEX"} / $$likelihoodsref{$subtype}{0}{"NONCOMP"};
		}
	}
	
	# calculate prior probabilities
	$$likelihoodsref{"PRIORCOMP"} = $complex_edges / $all_edges;  # multiply by a factor to estimate prior probability as higher than it actually is
	if ($$likelihoodsref{"PRIORCOMP"} > 1) { die; }
	$$likelihoodsref{"PRIORNONCOMP"} = 1 - $$likelihoodsref{"PRIORCOMP"};
		
	# print
	print "Complex edges = $complex_edges, noncomp edges = $noncomp_edges, all edges = $all_edges\n";
	print "Likelihoods:\n";
	print "Feat\tBin\tNf,c/Nc\tNf,~c/N~c\tNf/N\tLikRat\tLikRat*Prior\n";
	foreach my $subtype (sort keys %counts) {
		foreach my $bin (sort {$a<=>$b} keys %{$counts{$subtype}}) {
			print "$subtype\t$bin\t".
				($counts{$subtype}{$bin}{"COMPLEX"})."/".$complex_edges."=".($$likelihoodsref{$subtype}{$bin}{"COMPLEX"})."\t".
				($counts{$subtype}{$bin}{"NONCOMP"})."/".$noncomp_edges."=".($$likelihoodsref{$subtype}{$bin}{"NONCOMP"})."\t".
				($counts{$subtype}{$bin}{"ALL"})."/".$all_edges."=".($counts{$subtype}{$bin}{"ALL"} / $all_edges)."\t".
				($$likelihoodsref{$subtype}{$bin}{"LIKERATIO"})."\t".($$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}*$$likelihoodsref{"PRIORCOMP"})."\n";
		}
	}
	print "\n";
	
}



# $likelihoods{$subtype}{$bin}{0} = P(C=0|f)
# $likelihoods{$subtype}{$bin}{2} = P(C=2|f)
# $likelihoods{$subtype}{$bin}{4} = P(C=4|f)
# $likelihoods{$subtype}{$bin}{ALL} = P(f)
# $likelihoods{$subtype}{$bin}{LIKERATIO}{2} = likelihood ratio for class 2 = P(C=2|f) / P(C!=2|f)
# $likelihoods{$subtype}{$bin}{LIKERATIO}{4} = likelihood ratio for class 4 = P(C=4|f) / P(C!=4|f)
# $likelihoods{PRIORCOMP}{$subtype}{0} = P(C=0)
sub CalcLikelihoodsSSS($$$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	my $type_calc = $_[2];
	
	# $counts{$subtype}{$bin}{0} = # samples with that subtype & bin value and COCOMP_SIZE_CLASS = 0
	# $counts{$subtype}{$bin}{2} = # samples with that subtype & bin value and COCOMP_SIZE_CLASS = 2
	# $counts{$subtype}{$bin}{4} = # samples with that subtype & bin value and COCOMP_SIZE_CLASS = 4
	# $counts{$subtype}{$bin}{"ALL"} = # samples with that subtype & bin value
	# Bin 0 also represents all implicit edges
	my %counts = ();
	# $counts_priors{0} = # samples with COCOMP_SIZE_CLASS = 0
	# $counts_priors{2} = # samples with COCOMP_SIZE_CLASS = 2
	# $counts_priors{4} = # samples with COCOMP_SIZE_CLASS = 4
	# $counts_priors{ALL} = # all samples
	my %counts_priors = ();
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		
		if (!defined $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS}) {
			die "$key";
		}
		# counts for the cocomp classes
		my $cocomp_class = $$dataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS};
		$counts_priors{$cocomp_class}++;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$counts{$subtype}{$bin}{$cocomp_class}++;
		}
		
		# counts for "ALL" 
		$counts_priors{ALL}++;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$counts{$subtype}{$bin}{ALL}++;
		}
	}
	
	if ($COUNT_ZERO_EDGES==1 || $COUNT_ZERO_EDGES==3) {
		# set counts for bin 0 to include implicit edges
		foreach my $subtype (keys %counts) {
			$counts{$subtype}{0}{2} = $$dataref{NUMCOCOMP2EDGES};
			$counts{$subtype}{0}{4} = $$dataref{NUMCOCOMP4EDGES};
			$counts{$subtype}{0}{"ALL"} = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
			$counts{$subtype}{0}{0} = $counts{$subtype}{0}{"ALL"} - $counts{$subtype}{0}{2} - $counts{$subtype}{0}{4};
			
			foreach my $bin (keys %{$counts{$subtype}}) {
				if ($bin==0) { next; }
				$counts{$subtype}{0}{2} -= $counts{$subtype}{$bin}{2};
				$counts{$subtype}{0}{4} -= $counts{$subtype}{$bin}{4};
				$counts{$subtype}{0}{"ALL"} -= $counts{$subtype}{$bin}{"ALL"};
				$counts{$subtype}{0}{0} -= $counts{$subtype}{$bin}{0};
			}
		}
		$counts_priors{ALL} = $$dataref{NUMPROTEINS} * ($$dataref{NUMPROTEINS}-1) / 2;
		$counts_priors{2} = $$dataref{NUMCOCOMP2EDGES};
		$counts_priors{4} = $$dataref{NUMCOCOMP4EDGES};
		$counts_priors{0} = $counts_priors{ALL} - $counts_priors{2} - $counts_priors{4};
	}
	
	# check that the counts sum up correctly
	foreach my $subtype (keys %counts) {
		my $cocomp0_check = 0;
		my $cocomp2_check = 0;
		my $cocomp4_check = 0;
		my $all_edges_check = 0;
		foreach my $bin (keys %{$counts{$subtype}}) {
			if ($counts{$subtype}{$bin}{"ALL"} != $counts{$subtype}{$bin}{0} + $counts{$subtype}{$bin}{2} + $counts{$subtype}{$bin}{4}) {
				die; 
			}	
			$cocomp0_check += $counts{$subtype}{$bin}{0};
			$cocomp2_check += $counts{$subtype}{$bin}{2};
			$cocomp4_check += $counts{$subtype}{$bin}{4};
			$all_edges_check += $counts{$subtype}{$bin}{"ALL"};
		}
		if ($cocomp0_check != $counts_priors{0}) { die; }
		if ($cocomp2_check != $counts_priors{2}) { die; }
		if ($cocomp4_check != $counts_priors{4}) { die; }
		if ($all_edges_check != $counts_priors{ALL}) { die; }
	}
			
	
	# laplace smoothing
	# since different amounts are added to the priors of the various subtypes (dependent on number of bins in each subtype), have different prior counts for each subtype
	my %counts_priors_smoothed = ();
	foreach my $subtype (keys %counts) {
		$counts_priors_smoothed{$subtype}{0} = $counts_priors{0};
		$counts_priors_smoothed{$subtype}{2} = $counts_priors{2};
		$counts_priors_smoothed{$subtype}{4} = $counts_priors{4};
		$counts_priors_smoothed{$subtype}{ALL} = $counts_priors{ALL};
		foreach my $bin (keys %{$counts{$subtype}}) {
			$counts{$subtype}{$bin}{0}++;
			$counts{$subtype}{$bin}{2}++;
			$counts{$subtype}{$bin}{4}++;
			$counts{$subtype}{$bin}{"ALL"}+=3;
			$counts_priors_smoothed{$subtype}{0}++;
			$counts_priors_smoothed{$subtype}{2}++;
			$counts_priors_smoothed{$subtype}{4}++;
			$counts_priors_smoothed{$subtype}{ALL}+=3;
		}
	}
	
		
			
	# calculate likelihoods
	foreach my $subtype (keys %counts) {
		foreach my $bin (keys %{$counts{$subtype}}) {
			$$likelihoodsref{$subtype}{$bin}{0} = $counts{$subtype}{$bin}{0} / $counts_priors_smoothed{$subtype}{0};
			$$likelihoodsref{$subtype}{$bin}{2} = $counts{$subtype}{$bin}{2} / $counts_priors_smoothed{$subtype}{2};
			$$likelihoodsref{$subtype}{$bin}{4} = $counts{$subtype}{$bin}{4} / $counts_priors_smoothed{$subtype}{4};
			$$likelihoodsref{$subtype}{$bin}{"ALL"} = $counts{$subtype}{$bin}{"ALL"} / $counts_priors_smoothed{$subtype}{ALL};
#			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin}{2} / $$likelihoodsref{$subtype}{$bin}{0};
#			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin}{4} / $$likelihoodsref{$subtype}{$bin}{0};
			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin}{2} / (($counts{$subtype}{$bin}{0}+$counts{$subtype}{$bin}{4})/($counts_priors_smoothed{$subtype}{0}+$counts_priors_smoothed{$subtype}{4}));
			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin}{4} / (($counts{$subtype}{$bin}{0}+$counts{$subtype}{$bin}{2})/($counts_priors_smoothed{$subtype}{0}+$counts_priors_smoothed{$subtype}{2}));
		}
#		if ($COUNT_ZERO_EDGES==3) {
#			$$likelihoodsref{$subtype}{0}{0} = 1;
#			$$likelihoodsref{$subtype}{0}{2} = 1;
#			$$likelihoodsref{$subtype}{0}{4} = 1;
#			$$likelihoodsref{$subtype}{0}{"ALL"} = 1;
#			$$likelihoodsref{$subtype}{0}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{0}{2} / $$likelihoodsref{$subtype}{0}{0};
##			$$likelihoodsref{$subtype}{0}{LIKERATIO"}{4} = $$likelihoodsref{$subtype}{0}{4} / $$likelihoodsref{$subtype}{0}{0};
#			$$likelihoodsref{$subtype}{0}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{0}{2} / $$likelihoodsref{$subtype}{0}{0};
#			$$likelihoodsref{$subtype}{0}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{0}{4} / $$likelihoodsref{$subtype}{0}{0};
#		}
		# calculate prior probabilities
#		$$likelihoodsref{"PRIORCOMP"}{$subtype}{0} = $counts_priors_smoothed{$subtype}{0} / $counts_priors_smoothed{$subtype}{ALL};
#		$$likelihoodsref{"PRIORCOMP"}{$subtype}{2} = $counts_priors_smoothed{$subtype}{2} / $counts_priors_smoothed{$subtype}{ALL};
#		$$likelihoodsref{"PRIORCOMP"}{$subtype}{4} = $counts_priors_smoothed{$subtype}{4} / $counts_priors_smoothed{$subtype}{ALL};
	}
	# calculate prior probabilities
	$$likelihoodsref{"PRIORCOMP"}{0} = $counts_priors{0} / $counts_priors{ALL};
	$$likelihoodsref{"PRIORCOMP"}{2} = $counts_priors{2} / $counts_priors{ALL};
	$$likelihoodsref{"PRIORCOMP"}{4} = $counts_priors{4} / $counts_priors{ALL};
	
		
	# print
	print "Cocomp0 edges = $counts_priors{0}, Cocomp2 edges = $counts_priors{2}, Cocomp4 edges = $counts_priors{4}, all edges = $counts_priors{ALL}, \n";
	print "Likelihoods:\n";
	print "Feat\tBin\tNf,c0/Nc0\tNf,c2/Nc2\tNf,c4/Nc4\tNf/N\tLikRatC2\n";#\tLikRat\tLikRat*Prior\n";
	foreach my $subtype (sort keys %counts) {
		foreach my $bin (sort {$a<=>$b} keys %{$counts{$subtype}}) {
			print "$subtype\t$bin\t".
				($counts{$subtype}{$bin}{0})."/".$counts_priors_smoothed{$subtype}{0}."=".($$likelihoodsref{$subtype}{$bin}{0})."\t".
				($counts{$subtype}{$bin}{2})."/".$counts_priors_smoothed{$subtype}{2}."=".($$likelihoodsref{$subtype}{$bin}{2})."\t".
				($counts{$subtype}{$bin}{4})."/".$counts_priors_smoothed{$subtype}{4}."=".($$likelihoodsref{$subtype}{$bin}{4})."\t".
				($counts{$subtype}{$bin}{"ALL"})."/".$counts_priors_smoothed{$subtype}{ALL}."=".($counts{$subtype}{$bin}{"ALL"} / $counts_priors_smoothed{$subtype}{ALL})."\t".
				($$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2}).
				"\n";
		}
	}
	print "\n";
	
	
	# smoothen likelihoods: likelihoods must be non-decreasing, for PPIREL, STRING, PUBMED
	if ($SMOOTHEN_LIKELIHOODS==1) {
		SmoothenLikelihoods (\%counts, $likelihoodsref, "PPIREL", \%counts_priors_smoothed, 2);
		SmoothenLikelihoods (\%counts, $likelihoodsref, "STRING", \%counts_priors_smoothed, 2);
		SmoothenLikelihoods (\%counts, $likelihoodsref, "PUBMED", \%counts_priors_smoothed, 2);
		print "Final smoothened likelihoods (likelihoods must be non-decreasing)\n";
		print "Cocomp0 edges = $counts_priors{0}, Cocomp2 edges = $counts_priors{2}, Cocomp4 edges = $counts_priors{4}, all edges = $counts_priors{ALL}, \n";
		print "Likelihoods:\n";
		print "Feat\tBin\tNf,c0/Nc0\tNf,c2/Nc2\tNf,c4/Nc4\tNf/N\tLikRatC2\n";#\tLikRat\tLikRat*Prior\n";
		foreach my $subtype (sort keys %counts) {
			foreach my $bin (sort {$a<=>$b} keys %{$counts{$subtype}}) {
				print "$subtype\t$bin\t".
					($counts{$subtype}{$bin}{0})."/".$counts_priors_smoothed{$subtype}{0}."=".($$likelihoodsref{$subtype}{$bin}{0})."\t".
					($counts{$subtype}{$bin}{2})."/".$counts_priors_smoothed{$subtype}{2}."=".($$likelihoodsref{$subtype}{$bin}{2})."\t".
					($counts{$subtype}{$bin}{4})."/".$counts_priors_smoothed{$subtype}{4}."=".($$likelihoodsref{$subtype}{$bin}{4})."\t".
					($counts{$subtype}{$bin}{"ALL"})."/".$counts_priors_smoothed{$subtype}{ALL}."=".($counts{$subtype}{$bin}{"ALL"} / $counts_priors_smoothed{$subtype}{ALL})."\t".
					($$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2}).
					"\n";
			}
		}
		print "\n";
	} # end if SMOOTHEN_LIKELIHOODS==1
}



# smoothen likelihoods (ensure non-decreasing) for a particular class 
sub SmoothenLikelihoods {
	my $countsref = $_[0];
	my $likelihoodsref = $_[1];
	my $subtype = $_[2];
	my $countspriorsref = $_[3];	
	my $class = $_[4];

	my $finished = 0;
	my %merged_bins = (); # $merged_bins{1}{1} = 1, {1}{2} = 1, iff bins 1 and 2 were merged
	foreach my $bin (sort {$a<=>$b} keys %{$$countsref{$subtype}}) {
		$merged_bins{$bin}{$bin} = 1;
	}
	while ($finished==0) {
		$finished = 1;
		foreach my $bin (sort {$a<=>$b} keys %{$$countsref{$subtype}}) {
			if ($bin==0) { next; }
			if ($$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{$class} >= $$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{$class}) { next; }
			$$countsref{$subtype}{$bin}{0} += $$countsref{$subtype}{$bin-1}{0};
			$$countsref{$subtype}{$bin-1}{0} = $$countsref{$subtype}{$bin}{0};
			$$countsref{$subtype}{$bin}{2} += $$countsref{$subtype}{$bin-1}{2};
			$$countsref{$subtype}{$bin-1}{2} = $$countsref{$subtype}{$bin}{2};
			$$countsref{$subtype}{$bin}{4} += $$countsref{$subtype}{$bin-1}{4};
			$$countsref{$subtype}{$bin-1}{4} = $$countsref{$subtype}{$bin}{4};
			$$countsref{$subtype}{$bin}{"ALL"} += $$countsref{$subtype}{$bin-1}{"ALL"};
			$$countsref{$subtype}{$bin-1}{"ALL"} = $$countsref{$subtype}{$bin}{"ALL"};
			$$likelihoodsref{$subtype}{$bin}{0} = $$countsref{$subtype}{$bin}{0} / $$countspriorsref{$subtype}{0};
			$$likelihoodsref{$subtype}{$bin-1}{0} = $$countsref{$subtype}{$bin-1}{0} / $$countspriorsref{$subtype}{0};
			$$likelihoodsref{$subtype}{$bin}{2} = $$countsref{$subtype}{$bin}{2} / $$countspriorsref{$subtype}{2};
			$$likelihoodsref{$subtype}{$bin-1}{2} = $$countsref{$subtype}{$bin-1}{2} / $$countspriorsref{$subtype}{2};
			$$likelihoodsref{$subtype}{$bin}{4} = $$countsref{$subtype}{$bin}{4} / $$countspriorsref{$subtype}{4};
			$$likelihoodsref{$subtype}{$bin-1}{4} = $$countsref{$subtype}{$bin-1}{4} / $$countspriorsref{$subtype}{4};
			$$likelihoodsref{$subtype}{$bin}{ALL} = $$countsref{$subtype}{$bin}{ALL} / $$countspriorsref{$subtype}{ALL};
			$$likelihoodsref{$subtype}{$bin-1}{ALL} = $$countsref{$subtype}{$bin-1}{ALL} / $$countspriorsref{$subtype}{ALL};
#			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin}{2} / $$likelihoodsref{$subtype}{$bin}{0};
#			$$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin-1}{2} / $$likelihoodsref{$subtype}{$bin-1}{0};
			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin}{2} / (($$countsref{$subtype}{$bin}{0}+$$countsref{$subtype}{$bin}{4})/($$countspriorsref{$subtype}{0}+$$countspriorsref{$subtype}{4}));
			$$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin-1}{2} / (($$countsref{$subtype}{$bin-1}{0}+$$countsref{$subtype}{$bin-1}{4})/($$countspriorsref{$subtype}{0}+$$countspriorsref{$subtype}{4}));
			$$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin}{4} / (($$countsref{$subtype}{$bin}{0}+$$countsref{$subtype}{$bin}{2})/($$countspriorsref{$subtype}{0}+$$countspriorsref{$subtype}{2}));
			$$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin-1}{4} / (($$countsref{$subtype}{$bin-1}{0}+$$countsref{$subtype}{$bin-1}{2})/($$countspriorsref{$subtype}{0}+$$countspriorsref{$subtype}{2}));
			foreach my $bb (keys %{$merged_bins{$bin}}) {
				if ($bb == $bin) { next; }
				$$countsref{$subtype}{$bb}{0} = $$countsref{$subtype}{$bin}{0};
				$$countsref{$subtype}{$bb}{2} = $$countsref{$subtype}{$bin}{2};
				$$countsref{$subtype}{$bb}{4} = $$countsref{$subtype}{$bin}{4};
				$$countsref{$subtype}{$bb}{"ALL"} = $$countsref{$subtype}{$bin}{"ALL"};
				$$likelihoodsref{$subtype}{$bb}{0} = $$likelihoodsref{$subtype}{$bin}{0};
				$$likelihoodsref{$subtype}{$bb}{2} = $$likelihoodsref{$subtype}{$bin}{2};
				$$likelihoodsref{$subtype}{$bb}{4} = $$likelihoodsref{$subtype}{$bin}{4};
				$$likelihoodsref{$subtype}{$bb}{"ALL"} = $$likelihoodsref{$subtype}{$bin}{"ALL"};
				$$likelihoodsref{$subtype}{$bb}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{2};
				$$likelihoodsref{$subtype}{$bb}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin}{"LIKERATIO"}{4};
			}
			foreach my $bb (keys %{$merged_bins{$bin-1}}) {
				if ($bb == $bin-1) { next; }
				$$countsref{$subtype}{$bb}{0} = $$countsref{$subtype}{$bin-1}{0};
				$$countsref{$subtype}{$bb}{2} = $$countsref{$subtype}{$bin-1}{2};
				$$countsref{$subtype}{$bb}{4} = $$countsref{$subtype}{$bin-1}{4};
				$$countsref{$subtype}{$bb}{"ALL"} = $$countsref{$subtype}{$bin-1}{"ALL"};
				$$likelihoodsref{$subtype}{$bb}{0} = $$likelihoodsref{$subtype}{$bin-1}{0};
				$$likelihoodsref{$subtype}{$bb}{2} = $$likelihoodsref{$subtype}{$bin-1}{2};
				$$likelihoodsref{$subtype}{$bb}{4} = $$likelihoodsref{$subtype}{$bin-1}{4};
				$$likelihoodsref{$subtype}{$bb}{"ALL"} = $$likelihoodsref{$subtype}{$bin-1}{"ALL"};
				$$likelihoodsref{$subtype}{$bb}{"LIKERATIO"}{2} = $$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{2};
				$$likelihoodsref{$subtype}{$bb}{"LIKERATIO"}{4} = $$likelihoodsref{$subtype}{$bin-1}{"LIKERATIO"}{4};
			}
			foreach my $bb (keys %{$merged_bins{$bin}}) {
				foreach my $bb1 (keys %{$merged_bins{$bin-1}}) {
					$merged_bins{$bb}{$bb1} = 1;
				}
			}
			foreach my $bb1 (keys %{$merged_bins{$bin-1}}) {
				foreach my $bb (keys %{$merged_bins{$bin}}) {
					$merged_bins{$bb1}{$bb} = 1;
				}
			}
			$finished = 0;
		} # end foreach bin
	} # end while $finished==0
	
}


# ======================== Functions to score and print scores =======================

#--------------------- Score edges by posterior prob -----------------------
sub ScoreEdgesPosteriorProb($$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		my $numer = 1;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$numer *= $$likelihoodsref{$subtype}{$bin}{"COMPLEX"};
		}
		$numer = $numer ** $DEPENDENCE_FACTOR;
		$numer *= $$likelihoodsref{"PRIORCOMP"};
		my $denom = 1;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$denom *= $$likelihoodsref{$subtype}{$bin}{"NONCOMP"};
		}
		$denom = $denom ** $DEPENDENCE_FACTOR;
		$denom *= $$likelihoodsref{"PRIORNONCOMP"};
		$denom += $numer;
		if ($denom==0) { die; }
		my $score = $numer / $denom;
		if ($score==0 || $score>1) { 
			die "Edge $key, numer = $numer, denom = $denom, score = $score"; 
		}
		$$dataref{INSTANCES}{$key}{SCORE} = $score;
	}
}


# $likelihoods{$subtype}{$bin}{0} = P(f|C=0)
# $likelihoods{$subtype}{$bin}{ALL} = P(f)
# $likelihoods{PRIORCOMP}{0} = P(C=0)
sub ScoreEdgesPosteriorProbSSS($$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	
	my %numers = (); # $numers{$class} = P(f1|$class) * P(f2|$class) * ... * P($class)
	my @subtypes_sorted = sort keys %{$$dataref{SUBTYPES}};
	my @classes = (0, 2, 4);
	
	print "Calculating scores...\n";
#	print "Pair";
#	foreach my $class (@classes) {
#		foreach my $subtype (@subtypes_sorted) {
#			print "\t$subtype $class";
#		}
#		print "\t$class prior";
#	}
#	print "\n";
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
#		print "$key";
		my $denom = 0;
		foreach my $class (@classes) {
			my $numer = 1;
			foreach my $subtype (@subtypes_sorted) {
				my $bin = 0;
				if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
					$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
				}
				$numer *= $$likelihoodsref{$subtype}{$bin}{$class};
#				print "\t$$likelihoodsref{$subtype}{$bin}{$class}";
			}
			if ($SSS_PRIORS{$class}==-1) {
				$numer *= $$likelihoodsref{"PRIORCOMP"}{$class};
			}
			else {
				$numer *= $SSS_PRIORS{$class};
			}
#			print "\t$$likelihoodsref{PRIORCOMP}{$class}";
			$numers{$class} = $numer;
			$denom += $numer;
		}
#		print "\n";
		if ($denom==0) { die; }
		foreach my $class (@classes) {
			my $score = $numers{$class} / $denom;
			if ($score==0 || $score>1) { 
				die "Edge $key, class = $class, numer = $numers{$class}, denom = $denom, score = $score";  
			}
			$$dataref{INSTANCES}{$key}{SCORE}{$class} = $score;
		}
	}
}

#--------------------- Score edges by posterior prob -----------------------
sub ScoreEdgesPosteriorRatio($$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		my $numer = 1;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$numer *= $$likelihoodsref{$subtype}{$bin}{"COMPLEX"};
		}
		$numer = $numer ** $DEPENDENCE_FACTOR;
		$numer *= $$likelihoodsref{"PRIORCOMP"};
		my $denom = 1;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$denom *= $$likelihoodsref{$subtype}{$bin}{"NONCOMP"};
		}
		$denom = $denom ** $DEPENDENCE_FACTOR;
		$denom *= $$likelihoodsref{"PRIORNONCOMP"};
		if ($denom==0) { die; }
		my $score = $numer / $denom;
		$$dataref{INSTANCES}{$key}{SCORE} = $score;
	}
}

#--------------------- Score edges by likelihood ratio -----------------------
sub ScoreEdgesLikelihoodRatio($$) {
	my $dataref = $_[0];
	my $likelihoodsref = $_[1];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		my $score = 1;
		foreach my $subtype (keys %{$$dataref{SUBTYPES}}) {
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			$score *= $$likelihoodsref{$subtype}{$bin}{"LIKERATIO"};
		}
		$$dataref{INSTANCES}{$key}{SCORE} = $score;
		if ($score==0) {
			print "$key score = 0!!\n";
		}
	}
}

#--------------------- Score edges by their PPI values -----------------------
sub ScoreEdgesPPIValue($) {
	my $dataref = $_[0];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		my $score = 0;
		if (defined $$dataref{INSTANCES}{$key}{DATA}{"L1"} && $$dataref{INSTANCES}{$key}{DATA}{"L1"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"L1"};
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"L2"} && $$dataref{INSTANCES}{$key}{DATA}{"L2"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"L2"};
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"L12"} && $$dataref{INSTANCES}{$key}{DATA}{"L12"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"L12"};
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"PPI"} && $$dataref{INSTANCES}{$key}{DATA}{"PPI"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"PPI"};
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"PPIL2"} && $$dataref{INSTANCES}{$key}{DATA}{"PPIL2"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"PPIL2"};
		}
		else {
			die "Error: scoring edges by PPI values, by $key has no PPI value\n";
		}
		$$dataref{INSTANCES}{$key}{SCORE} = $score;
		if ($score==0) {
			die "$key score = 0!!\n";
		}
	}
}

#--------------------- Score edges with 1 -----------------------
sub ScoreEdgesNoWeight($) {
	my $dataref = $_[0];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		if (defined $$dataref{INSTANCES}{$key}{DATA}{"L12"} && $$dataref{INSTANCES}{$key}{DATA}{"L12"}!=0) {
			$$dataref{INSTANCES}{$key}{SCORE} = 1;
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"L1"} && $$dataref{INSTANCES}{$key}{DATA}{"L1"}!=0) {
			$$dataref{INSTANCES}{$key}{SCORE} = 1;
		}
		elsif (defined $$dataref{INSTANCES}{$key}{DATA}{"PPI"} && $$dataref{INSTANCES}{$key}{DATA}{"PPI"}!=0) {
			$$dataref{INSTANCES}{$key}{SCORE} = 1;
		}
		else {
			die "Error: scoring edges with 1, but $key has no PPI value\n";
		}
	}
}

#--------------------- Score edges by STRING values -----------------------
sub ScoreEdgesStringValue($) {
	my $dataref = $_[0];
	
	foreach my $key (keys %{$$dataref{INSTANCES}}) {
		my $score = 0;
		if (defined $$dataref{INSTANCES}{$key}{DATA}{"STRING"} && $$dataref{INSTANCES}{$key}{DATA}{"STRING"}!=0) {
			$score = $$dataref{INSTANCES}{$key}{DATA}{"STRING"};
		}
		else {
			die "Error: Scoring by String value, but no STRING score! $key";
		}
		$$dataref{INSTANCES}{$key}{SCORE} = $score;
		if ($score==0) {
			die "$key score = 0!!\n";
		}
		if ($score > 1) { die; }
	}
}


#--------------------- Print scores ----------------------
sub PrintScores($) {
	my $outputfile = $_[0];
	my $dataref = $_[1];
	open (OUTPUTFILE, ">$outputfile");
	
	foreach my $key (sort keys %{$$dataref{INSTANCES}}) {		
		print OUTPUTFILE "$key, ";
		foreach my $subtype (sort keys %{$$dataref{SUBTYPES}}) {
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				print OUTPUTFILE ($$dataref{INSTANCES}{$key}{DATA}{$subtype}).", ";
			}
			else {
				print OUTPUTFILE "0, ";
			}
		}
		if ($$dataref{INSTANCES}{$key}{CLASS}{"ISCOMPLEX"}==1) {
			print OUTPUTFILE "\'complex\'";
		}
		else {
			print OUTPUTFILE "\'noncomp\'";
		}
		if (defined $$dataref{INSTANCES}{$key}{SCORE}) {
			print OUTPUTFILE ", $$dataref{INSTANCES}{$key}{SCORE}";
		}
		else {
			print OUTPUTFILE ", 0";
		}
		print OUTPUTFILE "\n";
	}
	close OUTPUTFILE;	
}


#--------------------- Print scores in CMC format ----------------------
sub PrintScoresCMC($) {
	my $outputfile = $_[0];
	my $dataref = $_[1];
	open (OUTPUTFILE, ">$outputfile");
	
	foreach my $key (sort keys %{$$dataref{INSTANCES}}) {
		(my $id_a, my $id_b) = split(/\|/, $key);
		if (defined $$dataref{INSTANCES}{$key}{SCORE}) {
			print OUTPUTFILE "$id_a $id_b $$dataref{INSTANCES}{$key}{SCORE}\n";
		}
		else {
			print OUTPUTFILE "$id_a $id_b 0\n";
		}
	}
	close OUTPUTFILE;	
}


sub PrintScoresSSS($) {
	my $outputfile = $_[0];
	my $dataref = $_[1];
	open (OUTPUTFILE, ">$outputfile");

	my @pairs_sorted = sort keys %{$$dataref{INSTANCES}};
	if (!defined $$dataref{INSTANCES}{$pairs_sorted[0]}{SCORE}) { die "$pairs_sorted[0] no score"; }
	my %classes = %{$$dataref{INSTANCES}{$pairs_sorted[0]}{SCORE}};
	my $num_classes = scalar keys %classes;

	foreach my $key (@pairs_sorted) {
		(my $id_a, my $id_b) = split(/\|/, $key);
		if (!defined $$dataref{INSTANCES}{$key}{SCORE}) { die; }
		if (scalar keys %{$$dataref{INSTANCES}{$key}{SCORE}} != $num_classes) { die "$key has different number of classes than $num_classes"; }
		print OUTPUTFILE "$id_a\t$id_b";
		foreach my $class (sort keys %classes) {
			print OUTPUTFILE "\t$$dataref{INSTANCES}{$key}{SCORE}{$class}";
		}
		print OUTPUTFILE "\n";
	}
	close OUTPUTFILE;	
}


#--------------------- Print likelihood ratios and scores ----------------------
sub PrintEdgesLikelihoodRatios($) {
	my $outputfile = $_[0];
	my $dataref = $_[1];
	my $likelihoodsref = $_[2];
	open (OUTPUTFILE, ">$outputfile");
	
				
	foreach my $subtype (sort keys %{$$dataref{SUBTYPES}}) {
		foreach my $key (sort keys %{$$dataref{INSTANCES}}) {
			(my $id_a, my $id_b) = split(/\|/, $key);
			my $bin = 0;
			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
			}
			if ($bin==0) { next; }
			my $likerat = $$likelihoodsref{$subtype}{$bin}{"LIKERATIO"};
			print OUTPUTFILE "$id_a\t$id_b\t$subtype\t$likerat\n";
		}		
	}
	
	foreach my $key (sort keys %{$$dataref{INSTANCES}}) {
		(my $id_a, my $id_b) = split(/\|/, $key);
		my $score = 0;
		if (defined $$dataref{INSTANCES}{$key}{SCORE}) {
			$score = $$dataref{INSTANCES}{$key}{SCORE};
		}
		print OUTPUTFILE "$id_a\t$id_b\tSCORE\t$score\n";
	}
	close OUTPUTFILE;	
}
#sub PrintEdgesLikelihoodRatios($) {
#	my $outputfile = $_[0];
#	my $dataref = $_[1];
#	my $likelihoodsref = $_[2];
#	open (OUTPUTFILE, ">$outputfile");
#	
#	print OUTPUTFILE "int1\tint2\t";
#	foreach my $subtype (sort keys %{$$dataref{SUBTYPES}}) {
#		print OUTPUTFILE "$subtype\t";
#	}
#	print OUTPUTFILE "score\n";
#			
#	foreach my $key (sort keys %{$$dataref{INSTANCES}}) {
#		(my $id_a, my $id_b) = split(/\|/, $key);
#		print OUTPUTFILE "$id_a\t$id_b\t";
#		
#		foreach my $subtype (sort keys %{$$dataref{SUBTYPES}}) {
#			my $bin = 0;
#			if (defined $$dataref{INSTANCES}{$key}{DATA}{$subtype}) {
#				$bin = $$dataref{INSTANCES}{$key}{DATA}{$subtype} + 0;
#			}
#			my $likerat = $$likelihoodsref{$subtype}{$bin}{"LIKERATIO"};
#			print OUTPUTFILE "$likerat\t";
#		}		
#			
#		if (defined $$dataref{INSTANCES}{$key}{SCORE}) {
#			print OUTPUTFILE "$$dataref{INSTANCES}{$key}{SCORE}\n";
#		}
#		else {
#			print OUTPUTFILE "0\n";
#		}
#	}
#	close OUTPUTFILE;	
#}


sub CreateCurrData ($$$$$) {
	my $inputdataref = $_[0];
	my $currdataref = $_[1];
	my $traincompsref = $_[2];
	
	# copy inputdata to currdata
	foreach my $subtype (keys %{$$inputdataref{SUBTYPES}}) {
		$$currdataref{SUBTYPES}{$subtype}{TYPE} = $$inputdataref{SUBTYPES}{$subtype}{TYPE};
	}		
	foreach my $key (keys %{$$inputdataref{INSTANCES}}) {			
		foreach my $subtype (keys %{$$inputdataref{SUBTYPES}}) {
			if (!defined $$inputdataref{INSTANCES}{$key}{DATA}{$subtype}) { die; }
			$$currdataref{INSTANCES}{$key}{DATA}{$subtype} = $$inputdataref{INSTANCES}{$key}{DATA}{$subtype};
		}
		$$currdataref{INSTANCES}{$key}{CLASS}{ISCOMPLEX} = 0;
		$$currdataref{INSTANCES}{$key}{CLASS}{COCOMP_SIZE_CLASS} = 0;
	}
	
	# set ISCOMPLEX and COCOMP_SIZE_CLASS
	foreach my $comp (keys %{$traincompsref}) {
		my $compsize = scalar keys %{$$traincompsref{$comp}};
		my $cocomp_size_class = ($compsize==2 || $compsize==3) ? 2 : 4;
		foreach my $prot1 (keys %{$$traincompsref{$comp}}) {
			foreach my $prot2 (keys %{$$traincompsref{$comp}}) {
				if ($prot1 ge $prot2) { next; }
				if (defined $$currdataref{INSTANCES}{"$prot1|$prot2"}) {
					$$currdataref{INSTANCES}{"$prot1|$prot2"}{CLASS}{ISCOMPLEX} = 1;
					if ($cocomp_size_class > $$currdataref{INSTANCES}{"$prot1|$prot2"}{CLASS}{COCOMP_SIZE_CLASS}) {
						$$currdataref{INSTANCES}{"$prot1|$prot2"}{CLASS}{COCOMP_SIZE_CLASS} = $cocomp_size_class;
					}
				}
			}
		}
	}	
	
	
	# calculate NUMPROTEINS and NUMCOMPLEXEDGES
	my %all_proteins = ();
	foreach my $key (keys %{$$currdataref{INSTANCES}}) {
		(my $id1, my $id2) = split (/\|/, $key);
		$all_proteins{$id1} = 1;
		$all_proteins{$id2} = 1;
	}
	foreach my $comp (keys %{$traincompsref}) {
		foreach my $prot (keys %{$$traincompsref{$comp}}) {
			$all_proteins{$prot} = 1;
		}
	}
	my %all_compedges = ();
	my %all_cocomp2_edges = ();
	my %all_cocomp4_edges = ();
	foreach my $comp (keys %{$traincompsref}) {
		my $compsize = scalar keys %{$$traincompsref{$comp}};
		my $cocomp_size_class = ($compsize==2 || $compsize==3) ? 2 : 4;
		foreach my $prot1 (keys %{$$traincompsref{$comp}}) {
			foreach my $prot2 (keys %{$$traincompsref{$comp}}) {
				if ($prot2 le $prot1) { next; }
				$all_compedges{"$prot1|$prot2"} = 1;
				if ($cocomp_size_class==4) { 
					$all_cocomp4_edges{"$prot1|$prot2"} = 1; 
					delete $all_cocomp2_edges{"$prot1|$prot2"};
				}
				elsif ($cocomp_size_class==2 && !defined $all_cocomp4_edges{"$prot1|$prot2"}) {
					$all_cocomp2_edges{"$prot1|$prot2"} = 1; 
				}
			}
		}
	}
	$$currdataref{NUMPROTEINS} = scalar keys %all_proteins;
	$$currdataref{NUMCOMPLEXEDGES} = scalar keys %all_compedges;
	$$currdataref{NUMCOCOMP2EDGES} = scalar keys %all_cocomp2_edges;
	$$currdataref{NUMCOCOMP4EDGES} = scalar keys %all_cocomp4_edges;
				
	
	
	# verify and print some stats
	my $num_curr_comp_edges = 0;
	foreach my $key (keys %{$$currdataref{INSTANCES}}) {			
		if ($$currdataref{INSTANCES}{$key}{CLASS}{ISCOMPLEX}==1) { $num_curr_comp_edges++; }
	}
	print "CreateCurrData:\tinputdata ".(scalar keys %{$$inputdataref{INSTANCES}})." instances, $$inputdataref{NUMPROTEINS} numproteins, $$inputdataref{NUMCOMPLEXEDGES} numcomplexedges (incl nondata), $$inputdataref{NUM_NONCOMP_BORDER} num_noncomp_border\n";
	print "\t\t\t\t\t\t\t\tcurrdata: ".(scalar keys %{$$currdataref{INSTANCES}})." instances, $$currdataref{NUMPROTEINS} numproteins, $$currdataref{NUMCOMPLEXEDGES} numcomplexedges (incl nondata), $$currdataref{NUM_NONCOMP_BORDER} num_noncomp_border\n";
	print "\t\t\t\t\t\t\t\tNum current train complexes: ".(scalar keys %{$traincompsref}).", num current train complex data edges = $num_curr_comp_edges\n";
	print "\n";
}
			
			

sub ReadXValData ($$) {
	my $filename = $_[0];
	my $xvaldataref = $_[1]; # $xval_date{$iter}{$complexid} = 1, if $complexid is a TEST complex in $iter
	
	open (XVALFILE, "$filename") || die ("Cannot open xval sampling file $filename");
	my $curriter;
	foreach my $line (<XVALFILE>) {
		chomp $line;
		my @toks = split(/\t/, $line);
		if ($toks[0] eq "iter") {
			$curriter = $toks[1] + 0;
		}
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



	
## Should not add any no

# -------- Evaluate crossvalidation results
sub CalcPrecVsRecall ($$) {
	my $results_ref = $_[0];
	my $num_complex_edges = $_[1];
	my $num_neg_edges = $_[2];
	my $auc_ref = $_[3];
	my $fscore_ref = $_[4];
	my $aucroc_ref = $_[5];
		
	print "Score threshold\tPredictions\tRecall\tPrecision\tF1\tFPR\n";
	my $predicts = 0;		
	my $corrects= 0;
	my $lastvalue = -1;
	my $threshold = 0.01;
	my $recall = 0;
	my $auc = 0;
	my $auc_prevrecall = 0;
	my $aucroc = 0;
	my $aucroc_prevFPR = 0;
	@$results_ref = sort {$$b[0] <=> $$a[0]} @$results_ref;
	foreach (@$results_ref) {
		if ($lastvalue != $$_[0]) {
	   	$recall = $corrects/$num_complex_edges;
	   	if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)    {
	   		my $FPR = ($predicts - $corrects) / $num_neg_edges;
	   		my $precision = $corrects/$predicts;
	   		my $fscore = 2 * $precision * $recall / ($precision + $recall);
	    	print "$lastvalue\t$predicts\t$recall\t$precision\t$fscore\t$FPR\n";
	    	while($threshold < $recall) {
	     		$threshold += 0.01;
	    	}
			  # accumulate AUC
			  $auc += ($recall - $auc_prevrecall) * ($precision);
			  $auc_prevrecall = $recall;
			  $aucroc += ($FPR - $aucroc_prevFPR) * ($recall + $auc_prevrecall)/2;
			  $aucroc_prevFPR = $FPR;
	   	}
	   	$lastvalue = $$_[0];
	  }
	  $corrects += $$_[1];
	  $predicts++;
	}
	if ($predicts > 0) {
		$recall = $corrects/$num_complex_edges;
	  my $FPR = ($predicts - $corrects) / $num_neg_edges;
	  my $precision = $corrects/$predicts;
	  my $fscore = 2 * $precision * $recall / ($precision + $recall);
	  $$fscore_ref = $fscore;
		print "$lastvalue\t$predicts\t$recall\t$precision\t$fscore\t$FPR\n";
	}
	else {
		$recall = $corrects/$num_complex_edges;
	  my $FPR = ($predicts - $corrects) / $num_neg_edges;
	  my $fscore = 0;
	  $$fscore_ref = $fscore;
		print "$lastvalue\t$predicts\t$recall\t0\t0\t$FPR\n";
	}
	# accumulate AUC
	my $precision = $corrects/$predicts;
	my $FPR = ($predicts - $corrects) / $num_neg_edges;
	$auc += ($recall - $auc_prevrecall) * $precision;
	$aucroc += (1 - $FPR) * (1 + $recall)/2;
	print "\nAUC\t$auc\n";
	print "ROCAUC\t$aucroc\n";
	print "\n\n";
	$$auc_ref = $auc;
	$$aucroc_ref = $aucroc;
	
}
	

sub CalcPrecVsRecallAllIters ($$) {
	my $results_ref = $_[0]; #$xval_results_iters{$iter}[i][0] = score, [1] = iscomplex?
	my $num_complex_edges_ref = $_[1];
	my $num_neg_edges_ref = $_[2];
	my $auc_ref = $_[3];
	my $fscore_ref = $_[4];
	my $aucroc_ref = $_[5];
	
	my @total_results = ();
	foreach my $iter (keys %{$results_ref}) {
		foreach (@{$$results_ref{$iter}}) {
			my @tmparray = ();
			$tmparray[0] = $$_[0];
			$tmparray[1] = $$_[1];
			push (@total_results, \@tmparray);
		}
	}
	my $total_num_complex_edges = 0;
	foreach my $iter (keys %{$num_complex_edges_ref}) {
		$total_num_complex_edges += $$num_complex_edges_ref{$iter};
	}
	my $total_num_neg_edges = 0;
	foreach my $iter (keys %{$num_neg_edges_ref}) {
		$total_num_neg_edges += $$num_neg_edges_ref{$iter};
	}
		
	print "Score threshold\tPredictions\tRecall\tPrecision\tF1\tFPR\n";
	my $predicts = 0;		
	my $corrects= 0;
	my $lastvalue = -1;
	my $threshold = 0.01;
	my $recall = 0;
	my $auc = 0;
	my $auc_prevrecall = 0;
	my $aucroc = 0;
	my $aucroc_prevFPR = 0;
	@total_results = sort {$$b[0] <=> $$a[0]} @total_results;
	foreach (@total_results) {
		if ($lastvalue != $$_[0]) {
	   	$recall = $corrects/$total_num_complex_edges;
	   	if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)    {
	   		my $FPR = ($predicts - $corrects) / $total_num_neg_edges;
	   		my $precision = $corrects/$predicts;
	   		my $fscore = 2 * $precision * $recall / ($precision + $recall);
	    	print "$lastvalue\t$predicts\t$recall\t$precision\t$fscore\t$FPR\n";
	    	while($threshold < $recall) {
	     		$threshold += 0.01;
	    	}
			  # accumulate AUC
			  $auc += ($recall - $auc_prevrecall) * ($precision);
			  $auc_prevrecall = $recall;
			  $aucroc += ($FPR - $aucroc_prevFPR) * ($recall + $auc_prevrecall)/2;
			  $aucroc_prevFPR = $FPR;
	   	}
	   	$lastvalue = $$_[0];
	  }
	  $corrects += $$_[1];
	  $predicts++;
	}
	if ($predicts > 0) {
		$recall = $corrects/$total_num_complex_edges;
	  my $FPR = ($predicts - $corrects) / $total_num_neg_edges;
	  my $precision = $corrects/$predicts;
	  my $fscore = 2 * $precision * $recall / ($precision + $recall);
	  $$fscore_ref = $fscore;
		print "$lastvalue\t$predicts\t$recall\t$precision\t$fscore\t$FPR\n";
	}
	else {
		$recall = $corrects/$total_num_complex_edges;
	  my $FPR = ($predicts - $corrects) / $total_num_neg_edges;
	  my $fscore = 0;
	  $$fscore_ref = $fscore;
		print "$lastvalue\t$predicts\t$recall\t0\t0\t$FPR\n";
	}
	# accumulate AUC
	my $precision = $corrects/$predicts;
	my $FPR = ($predicts - $corrects) / $total_num_neg_edges;
	$auc += ($recall - $auc_prevrecall) * $precision;
	$aucroc += (1 - $FPR) * (1 + $recall)/2;
	print "\nAUC\t$auc\n";
	print "ROCAUC\t$aucroc\n";
	print "\n\n";
	$$auc_ref = $auc;
	$$aucroc_ref = $aucroc;
	
}
	


