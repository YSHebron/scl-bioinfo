use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use List::Util 'min';
use POSIX;
use Getopt::Std;



# inputs:
# -i <input filename of scored edges. If 4 scores per line, taken as size-specific scoring.>
# -k <number of top scoring small-comp edges to keep>
# -p <features of edges. If defined, will use ppi edges only for finding triangles/edges and for calculating isolatedness>
# -q <score threshold for ppi edges, if using ppi (-p used). Use either -q or -r.>
# -r <number of top ppi edges to use, if using ppi (-p used). Use either -q or -r.>
# -c <use cohesive factor to score? 1 = yes (default), 0 = no, 2 = only use smallscore (skip largescore)>
# -o <output clusters filename>

my %argopts;
if (! getopts('i:k:o:p:q:r:c:', \%argopts)) { die "invalid arguments"; }

# 1, 0, 0 = default for addtri
my %ADDED_TRIANGLES = ();
my $ADD_TRIANGLES = 1;
my $addtri_penalty = 1; # 1 = no penalty, "PPI" = penalize with PPI score
#my $addtri_penalty = "PPI";
my $addtri_missingedge_PPIthr = 0; # since PPI data minimum is .8, only values above .8 would make a difference 
my $addtri_missingedge_SMthr = 0;
my $addtri_presentedges_SMthr = 0;  

my $SCORE_COFAC = 1;
if (defined $argopts{'c'}) {
	$SCORE_COFAC = $argopts{'c'}+0;
}
	

my $CLUS_SCORE_THR = 1e-8;	
my $SM_SCORE_THR = 0;


my $NUM_LARGECOMP_EDGES = 10000;
my $num_edges = 9999999999;
if (defined $argopts{'k'}) {
	#$NUM_LARGECOMP_EDGES = $argopts{'k'}+0;
	$num_edges = $argopts{'k'}+0;
}


# Default = mode 1
# 0 = clusters themselves must be PPI edges, during scoring neighbouring edges may not be PPI edges 
# 1 = clusters themselves must be PPI edges, during scoring neighbouring edges must be PPI edges
# 2 = clusters themselves must be PPI edges, during scoring neighbouring edges must be PPI for small scores, neighbouring edges may not be PPI for large scores
# 3 = clusters themselves need not be PPI edges, during scoring neighbouring edges must be PPI edges
my $USE_PPI_MODE = 1;  

my $inputfilename = $argopts{'i'};
my $largecompfilename = $argopts{'j'};





# Keep only PPI edges
my $USE_PPI = 0;   
my $ppi_thresh = 0;
my $num_top_ppis = 999999999;
my %EDGES_PPI = (); # PPI edges used for finding small comps. $EDGES_PPI{$prot1|$prot2} = ppi score
my %NEIGHBOURS_PPI = ();
my %EDGES_ALLPPI = (); # all PPI edges. Used to add missing edges to complete triangles.
if (defined $argopts{'p'}) {
	$USE_PPI = 1;
	open (PPIFILE, $argopts{'p'}) || die $!;
	foreach my $line (<PPIFILE>) {
		chomp $line;
		my @toks = split(/\s/, $line);
		my $prot1 = $toks[0];
		my $prot2 = $toks[1];
		if ($prot1 eq $prot2) { next; }			
		my $edge = $prot1 lt $prot2 ? "$prot1|$prot2" : "$prot2|$prot1";
		my $score;
		if (scalar @toks == 3) {
			$score = $toks[2]+0;
		}
		elsif (scalar @toks == 4) {
			my $datatype = $toks[2];
			if ($datatype ne "PPIREL") { next; }
			$score = $toks[3]+0;
		}
		$EDGES_ALLPPI{$edge} = $score+0;
	}
	close PPIFILE;
	
	# keep only top PPIs
	if (defined $argopts{'q'} && defined $argopts{'r'}) { die; }
	if (defined $argopts{'q'}) {
		$ppi_thresh = $argopts{'q'}+0;
		foreach my $edge (keys %EDGES_ALLPPI) {
			if ($EDGES_ALLPPI{$edge} >= $ppi_thresh) {
				my $score = $EDGES_ALLPPI{$edge};
				(my $prot1, my $prot2) = split(/\|/, $edge);
				$EDGES_PPI{$edge} = $score+0;
				$NEIGHBOURS_PPI{$prot1}{$prot2} = $score+0;
				$NEIGHBOURS_PPI{$prot2}{$prot1} = $score+0;				
			}
		}
	}
	elsif (defined $argopts{'r'}) {
		$num_top_ppis = $argopts{'r'}+0;
		foreach my $edge (sort {$EDGES_ALLPPI{$b} <=> $EDGES_ALLPPI{$a}} keys %EDGES_ALLPPI) {
			if (scalar keys %EDGES_PPI >= $num_top_ppis) { last; }
			my $score = $EDGES_ALLPPI{$edge};
			(my $prot1, my $prot2) = split(/\|/, $edge);
			$EDGES_PPI{$edge} = $score+0;
			$NEIGHBOURS_PPI{$prot1}{$prot2} = $score+0;
			$NEIGHBOURS_PPI{$prot2}{$prot1} = $score+0;				
		}
	}
	else {
		foreach my $edge (keys %EDGES_ALLPPI) {
			my $score = $EDGES_ALLPPI{$edge};
			(my $prot1, my $prot2) = split(/\|/, $edge);
			$EDGES_PPI{$edge} = $score+0;
			$NEIGHBOURS_PPI{$prot1}{$prot2} = $score+0;
			$NEIGHBOURS_PPI{$prot2}{$prot1} = $score+0;				
		}
	}
	
	
	print "Num PPIs read for finding small comps = ".(scalar keys %EDGES_PPI).", all PPIs = ".(scalar keys %EDGES_ALLPPI)."\n";
}



# read the scored edges
# %scored_edges{$prot1|$prot2} = class2 score (cocomp2/3)
my %scored_edges = ();
my %scored_edges_largecomp = (); # %scored_edges_largecomp{$prot1|$prot2} = largecomp score
open (INPUTFILE, "$inputfilename") || die $!;
my $num_edges_read = 0;
foreach my $line (<INPUTFILE>) {
	chomp $line;
	my @toks = split(/\s/, $line);
	$num_edges_read++;
	my ($id_a, $id_b, $score0, $score1, $score2, $score3, $score4, $score23);
	if (scalar @toks == 6) {
		($id_a, $id_b, $score0, $score2, $score3, $score4) = @toks;
		$score1 = $score2+$score3+$score4;
		$score23 = $score2 + $score3;
	}
	elsif (scalar @toks == 5) {
		($id_a, $id_b, $score0, $score2, $score4) = @toks;
		$score1 = $score2+$score4;
		$score3 = 0;
		$score23 = $score2;			
	}
	elsif (scalar @toks==3) {
		($id_a, $id_b, $score1) = @toks;
		$score2 = $score1;
		$score3 = $score1;
		$score4 = $score1;
		$score23 = $score1;
	}
	else { 	die; 	}
	if ($score4 >= 1e-6) {
		my $edge = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
		$scored_edges_largecomp{$edge} = $score4+0;		
	}
	if ($score23 >= $SM_SCORE_THR) { 
		my $edge = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
		$scored_edges{$edge} = $score23+0;
	}
}
print "num edges read = $num_edges_read, num smallcomp edges kept = ".(scalar keys %scored_edges).", num largecomp edges kept = ".(scalar keys %scored_edges_largecomp)."\n";
close INPUTFILE;
	
	
	
# find the top k edges with highest small comp scores
my %top_edges = (); # $top_edges{$pro1|$prot2} = smallcomp score
foreach my $edge (sort {$scored_edges{$b} <=> $scored_edges{$a}} keys %scored_edges) {
	if (scalar keys %top_edges >= $num_edges) { last; }
	if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2) && !defined $EDGES_PPI{$edge}) { next; }
	$top_edges{$edge} = $scored_edges{$edge};
}
print "num top edges small comps = ".(scalar keys %top_edges)."\n";


# Create one just for size-3 clusters, because when adding triangles we will add missing edges to it
my %top_edges3 = ();
foreach my $edge (keys %top_edges) {
	$top_edges3{$edge} = $top_edges{$edge};
}
print "num top edges for size-3 comps = ".(scalar keys %top_edges3)."\n";
	
# keep only the top $NUM_LARGECOMP_EDGES edges with highest LARGE comp scores
my %top_edges_largecomp = (); # $top_edges_largecomp{$pro1|$prot2} = largecomp score
foreach my $edge (sort {$scored_edges_largecomp{$b} <=> $scored_edges_largecomp{$a}} keys %scored_edges_largecomp) {
	if (scalar keys %top_edges_largecomp >= $NUM_LARGECOMP_EDGES) { last; }
	if ($USE_PPI==1 && $USE_PPI_MODE==1 && !defined $EDGES_PPI{$edge}) { next; }
	$top_edges_largecomp{$edge} = $scored_edges_largecomp{$edge};
}
%scored_edges_largecomp = (); # clear memory
print "num top edges from largecomp = ".(scalar keys %top_edges_largecomp)."\n";

	

# create neighbours
my %neighbours = (); # $neighbour2{$prot1}{$prot2} = 1
foreach my $edge (keys %top_edges) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	$neighbours{$prot1}{$prot2} = 1;
	$neighbours{$prot2}{$prot1} = 1;
}
		
# create neighbours for largecomp
my %neighbours_largecomp = (); # $neighbours_largecomp{$prot1}{$prot2} = 1
foreach my $edge (keys %top_edges_largecomp) {
	(my $prot1, my $prot2) = split(/\|/, $edge);
	$neighbours_largecomp{$prot1}{$prot2} = 1;
	$neighbours_largecomp{$prot2}{$prot1} = 1;
}
	
	
# find all triangles
my %triangles = (); # $triangles{$prot1|$prot2|$prot3} = triangle score
foreach my $edge (keys %top_edges) {
	my $edgescore = $top_edges{$edge};
	(my $prot1, my $prot2) = split(/\|/, $edge);
	if (scalar keys %{$neighbours{$prot1}} > scalar keys %{$neighbours{$prot2}}) {
		my $tmpprot = $prot1;
		$prot1 = $prot2;
		$prot2 = $tmpprot;
	}
	foreach my $nb (keys %{$neighbours{$prot1}}) {
		if ($nb eq $prot2) { next; }
		if (defined $neighbours{$nb}{$prot2}) {
			my $triplet = join('|', sort ($prot1, $prot2, $nb));
			my $edge2 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			my $edge3 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($top_edges{$edge2}==0 || $top_edges{$edge3}==0) { die; }
			my $score= ($edgescore + $top_edges{$edge2} + $top_edges{$edge3})/3; 
			$triangles{$triplet} = $score;
		}
	}
}
print "num triangles = ".(scalar keys %triangles)."\n";
	 
# try to add more triangles: those with 1 missing edge because PPI score not high enough, but its still relatively high, and SM score is high enough
%ADDED_TRIANGLES = (); # keep the PPI reliability of the missing edge in the triangle, to penalize it later
if ($ADD_TRIANGLES==1) {
	foreach my $edge (keys %top_edges) {
		(my $prot1, my $prot2) = split(/\|/, $edge);
		if ($top_edges{$edge} < $addtri_presentedges_SMthr) { next; }  # present edge SM score must exceed threshold
		foreach my $nb (keys %{$neighbours{$prot1}}) {
			if ($nb eq $prot2) { next; }
			if (defined $neighbours{$nb}{$prot2}) { next; }
			my $edge2 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			if ($top_edges{$edge2} < $addtri_presentedges_SMthr) { next; } # present edge SM score must exceed threshold
			my $missing_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if (!defined $scored_edges{$missing_edge} || $scored_edges{$missing_edge} < $addtri_missingedge_SMthr) { next; }  # missing edge SM score must exceed threshold
			if (!defined $EDGES_ALLPPI{$missing_edge} || $EDGES_ALLPPI{$missing_edge} < $addtri_missingedge_PPIthr) { next; } # missing edge PPI score must exceed threshold
			my $triplet = join('|', sort ($prot1, $prot2, $nb));
			my $score = ($top_edges{$edge} + $top_edges{$edge2} + $scored_edges{$missing_edge}) / 3;
			$triangles{$triplet} = $score;
			$top_edges3{$missing_edge} = $scored_edges{$missing_edge};
#			$top_edges3{$missing_edge} = 0;
			$ADDED_TRIANGLES{$triplet}{"PPI"} = $EDGES_ALLPPI{$missing_edge};
			$ADDED_TRIANGLES{$triplet}{"SM"} = $scored_edges{$missing_edge};
		}
		foreach my $nb (keys %{$neighbours{$prot2}}) {
			if ($nb eq $prot1) { next; }
			if (defined $neighbours{$nb}{$prot1}) { next; }
			my $edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($top_edges{$edge2} < $addtri_presentedges_SMthr) { next; } # present edge SM score must exceed threshold
			my $missing_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			if (!defined $scored_edges{$missing_edge} || $scored_edges{$missing_edge} < $addtri_missingedge_SMthr) { next; }  # missing edge SM score must exceed threshold
			if (!defined $EDGES_ALLPPI{$missing_edge} || $EDGES_ALLPPI{$missing_edge} < $addtri_missingedge_PPIthr) { next; } # missing edge PPI score must exceed threshold
			my $triplet = join('|', sort ($prot1, $prot2, $nb));
			my $score = ($top_edges{$edge} + $top_edges{$edge2} + $scored_edges{$missing_edge}) / 3;
			$triangles{$triplet} = $score;
			$top_edges3{$missing_edge} = $scored_edges{$missing_edge};
#			$top_edges3{$missing_edge} = 0;
			$ADDED_TRIANGLES{$triplet}{"PPI"} = $EDGES_ALLPPI{$missing_edge};
			$ADDED_TRIANGLES{$triplet}{"SM"} = $scored_edges{$missing_edge};
		}
	}
	print "Added more triangles by filling missing edges. Added triangles = ".(scalar keys %ADDED_TRIANGLES).", num triangles = ".(scalar keys %triangles).", num top edges for small comps = ".(scalar keys %top_edges)."\n";
}
%scored_edges = (); # clear memory
	
 
my %clusters2 = ();		
my %clusters3 = ();
	 

print "\n\n============ all edges & triangles, recalcprob, including non-triangle ext edges, then score by cohesive factor, USE_PPI = $USE_PPI, ppi_thresh = $ppi_thresh, num_top_ppis = $num_top_ppis, use ppi mode = $USE_PPI_MODE ===================================================\n";

GetClusters2(\%top_edges, \%neighbours, \%top_edges_largecomp, \%neighbours_largecomp, \%clusters2);
GetClusters3(\%top_edges3, \%triangles, \%neighbours, \%top_edges_largecomp, \%neighbours_largecomp, \%clusters3);
	
# print clusters
if (defined $argopts{'o'}) {
	my $outputfile = $argopts{'o'};
	open (OUTPUTCLUSFILE, ">$outputfile") || die $!;
	foreach my $clus (sort {$clusters2{$b}{SCORE} <=> $clusters2{$a}{SCORE}} keys %clusters2) {
		my @prots = sort keys %{$clusters2{$clus}{E}};
		my $protsstring = join (" ", @prots);
		my $size = scalar @prots;
		print OUTPUTCLUSFILE "F$clus(".$size."_$clusters2{$clus}{SCORE}): $protsstring\n";
	}
	foreach my $clus (sort {$clusters3{$b}{SCORE} <=> $clusters3{$a}{SCORE}} keys %clusters3) {
		my @prots = sort keys %{$clusters3{$clus}{E}};
		my $protsstring = join (" ", @prots);
		my $size = scalar @prots;
		print OUTPUTCLUSFILE "F$clus(".$size."_$clusters3{$clus}{SCORE}): $protsstring\n";
	}
	close OUTPUTCLUSFILE;
}
	
 
 
sub GetClusters3 {
	my $edges_ref = $_[0]; 			# $edges{$prot1|$prot2} = score
	my $triangles_ref = $_[1]; 	# $triangles{$prot1|$prot2|$prot3} = score
	my $neighbours_ref = $_[2]; # $neighbours{$prot1}{$prot2} = 1
	my $edges_largecomp_ref = $_[3]; 			# $edges_largecomp{$prot1|$prot2} = score
	my $neighbours_largecomp_ref = $_[4]; # $neighbours_largecomp{$prot1}{$prot2} = 1	
	my $clus_ref = $_[5];  			# clus{$clus}{SCORE} = score, {E}{$prot} = 1
 
	# recalculate the probabilities (smallscore), non-triangle ext edge counted too (missing edge score = 0.001), then score by cohesive factor (small+largescore)
		
	my $MISSING_EDGE_PROB = 0;
	# find triangles
	foreach my $tri (keys %{$triangles_ref}) {
		(my $prot1, my $prot2, my $prot3) = split(/\|/, $tri);
		my $edge1 = "$prot1|$prot2";
		my $edge2 = "$prot1|$prot3";
		my $edge3 = "$prot2|$prot3";
		if (!defined $$edges_ref{$edge1} || !defined $$edges_ref{$edge2} || !defined $$edges_ref{$edge3}) { die; }
		# see if any triangles are incident on edge1
		my $edge1score = $$edges_ref{$edge1};
		foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
			if ($nb eq $prot2 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb})) {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				my $ext_edge2_score = $$edges_ref{$ext_edge2};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
				$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $ext_edge2_score;
				if ($edge1score < 0) {
					$edge1score = 0;
				}
			}
			else {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
				if ($edge1score < 0) {
					$edge1score = 0;
				}
			}
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
			if ($nb eq $prot1 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb})) { next; }
			my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
			my $ext_edge1_score = $$edges_ref{$ext_edge1};
			if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
			$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
			if ($edge1score < 0) {
				$edge1score = 0;
			}
		}
		# see if any triangles are incident on edge2
		my $edge2score = $$edges_ref{$edge2};
		foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
			if ($nb eq $prot2 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot3}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb})) {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				my $ext_edge2_score = $$edges_ref{$ext_edge2};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
				$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $ext_edge2_score;
				if ($edge2score < 0) {
					$edge2score = 0;
				}
			}
			else {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
				if ($edge2score < 0) {
					$edge2score = 0;
				}
			}
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
			if ($nb eq $prot1 || $nb eq $prot2) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb})) {  next; }
			my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
			if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
			my $ext_edge1_score = $$edges_ref{$ext_edge1};
			if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
			$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
			if ($edge2score < 0) {
				$edge2score = 0;
			}
		}
		# see if any triangles are incident on edge3			
		my $edge3score = $$edges_ref{$edge3};
		foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
			if ($nb eq $prot1 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }					
			if (defined $$neighbours_ref{$prot3}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb})) {
				my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				my $ext_edge2_score = $$edges_ref{$ext_edge2};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
				$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $ext_edge2_score;
				if ($edge3score < 0) {
					$edge3score = 0;
				}
			}
			else {
				my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
				if ($edge3score < 0) {
					$edge3score = 0;
				}
			}
		}		
		foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
			if ($nb eq $prot1 || $nb eq $prot2) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb})) { next; }
			my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
			if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
			my $ext_edge1_score = $$edges_ref{$ext_edge1};
			if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
			$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
			if ($edge3score < 0) {
				$edge3score = 0;
			}
		}
		if ($SCORE_COFAC==0) {
			my $triscore = ($edge1score + $edge2score + $edge3score)/3;				
			if (defined $ADDED_TRIANGLES{$tri}) { 	
				if ($addtri_penalty eq "PPI") {
					$triscore *= $ADDED_TRIANGLES{$tri}{PPI};
				}
				else {	$triscore *= $addtri_penalty;  }
			}		
			if ($triscore >= $CLUS_SCORE_THR) { 
				$$clus_ref{$tri}{SCORE} = $triscore;
				$$clus_ref{$tri}{E}{$prot1} = 1;
				$$clus_ref{$tri}{E}{$prot2} = 1;
				$$clus_ref{$tri}{E}{$prot3} = 1;
				next;			
			}
		}
		# sum up external edges			
		my $extscores = 0;
		if ($SCORE_COFAC==2) {
			foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot1}}) {
				if ($nb eq $prot2 || $nb eq $prot3) { next; }
				if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
				my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
				$extscores += $$edges_largecomp_ref{$ext_edge};
			}
			foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot2}}) {
				if ($nb eq $prot1 || $nb eq $prot3) { next; }
				if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
				my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
				$extscores += $$edges_largecomp_ref{$ext_edge};
			}
			foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot3}}) {
				if ($nb eq $prot1 || $nb eq $prot2) { next; }
				if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb}) { next; }
				my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
				if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
				$extscores += $$edges_largecomp_ref{$ext_edge};
			}
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
			if ($nb eq $prot2 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
			$extscores += $$edges_ref{$ext_edge};
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
			if ($nb eq $prot1 || $nb eq $prot3) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
			$extscores += $$edges_ref{$ext_edge};
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
			if ($nb eq $prot1 || $nb eq $prot2) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot3}{$nb}) { next; }
			my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
			if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
			$extscores += $$edges_ref{$ext_edge};
		}
		my $triscore = (($edge1score + $edge2score + $edge3score) / ($edge1score + $edge2score + $edge3score + $extscores)) * (($edge1score + $edge2score + $edge3score)/3);
		if (defined $ADDED_TRIANGLES{$tri}) { 	
			if ($addtri_penalty eq "PPI") {
				$triscore *= $ADDED_TRIANGLES{$tri}{PPI};
			}
			else {	$triscore *= $addtri_penalty;  }
		}
		if ($triscore >= $CLUS_SCORE_THR) { 
			$$clus_ref{$tri}{SCORE} = $triscore;
			$$clus_ref{$tri}{E}{$prot1} = 1;
			$$clus_ref{$tri}{E}{$prot2} = 1;
			$$clus_ref{$tri}{E}{$prot3} = 1;
		}
	}
}


sub GetClusters2 {
	my $edges_ref = $_[0]; 			# $edges{$prot1|$prot2} = score
	my $neighbours_ref = $_[1]; # $neighbours{$prot1}{$prot2} = 1
	my $edges_largecomp_ref = $_[2]; 			# $edges_largecomp{$prot1|$prot2} = score
	my $neighbours_largecomp_ref = $_[3]; # $neighbours_largecomp{$prot1}{$prot2} = 1	
	my $clus_ref = $_[4];  			# clus{$clus}{SCORE} = score, {E}{$prot} = 1
	
	my $MISSING_EDGE_PROB = 0;
	# find edges
	foreach my $edge (keys %{$edges_ref}) {
		(my $prot1, my $prot2) = split(/\|/, $edge);
		my $newscore = $$edges_ref{$edge};
		foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
			if ($nb eq $prot2) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb})) {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				my $ext_edge2_score = $$edges_ref{$ext_edge2};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
				$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $ext_edge2_score;
				if ($newscore < 0) {
					$newscore = 0;
				}
			}
			else {
				my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
				my $ext_edge1_score = $$edges_ref{$ext_edge1};
				if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
				$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
				if ($newscore < 0) {
					$newscore = 0;
				}
			}
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
			if ($nb eq $prot1) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
			if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb})) { next; }
			my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
			my $ext_edge1_score = $$edges_ref{$ext_edge1};
			if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
			$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
			if ($newscore < 0) {
				$newscore = 0;
			}
		}			
		if ($SCORE_COFAC==0) {
			if ($newscore >= $CLUS_SCORE_THR) {
				$$clus_ref{$edge}{SCORE} = $newscore;
				$$clus_ref{$edge}{E}{$prot1} = 1;
				$$clus_ref{$edge}{E}{$prot2} = 1;
				next;			
			}
		}
		# sum up external edges			
		my $extscores = 0;
		if ($SCORE_COFAC==2) {
			foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot1}}) {
				if ($nb eq $prot2) { next; }
				if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
				my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
				if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
				$extscores += $$edges_largecomp_ref{$ext_edge};
			}
			foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot2}}) {
				if ($nb eq $prot1) { next; }
				if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
				my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
				if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
				$extscores += $$edges_largecomp_ref{$ext_edge};
			}
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
			if ($nb eq $prot2) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot1}{$nb}) { next; }
			my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
			if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
			$extscores += $$edges_ref{$ext_edge};
		}
		foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
			if ($nb eq $prot1) { next; }
			if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI{$prot2}{$nb}) { next; }
			my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
			if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
			$extscores += $$edges_ref{$ext_edge};
		}
		my $edgescore = ($newscore / ($newscore + $extscores)) * $newscore;
		if ($edgescore >= $CLUS_SCORE_THR) {
			$$clus_ref{$edge}{SCORE} = $edgescore;
			$$clus_ref{$edge}{E}{$prot1} = 1;
			$$clus_ref{$edge}{E}{$prot2} = 1;
		}
	}
}
		

	

#sub GetClusters {
#	my $edges_ref = $_[0]; 			# $edges{$prot1|$prot2} = score
#	my $triangles_ref = $_[1]; 	# $triangles{$prot1|$prot2|$prot3} = score
#	my $neighbours_ref = $_[2]; # $neighbours{$prot1}{$prot2} = 1
#	my $edges_largecomp_ref = $_[3]; 			# $edges_largecomp{$prot1|$prot2} = score
#	my $neighbours_largecomp_ref = $_[4]; # $neighbours_largecomp{$prot1}{$prot2} = 1	
#	my $clus_ref = $_[5];  			# clus{$clus}{SCORE} = score, {E}{$prot} = 1
#	my $clus2_ref = $_[6];
#	my $clus3_ref = $_[7];
#	my $mode = $_[8];
#	my $only_do_compsize = $_[9]; # "comp2" or "comp3"
#	
#	
#	
#	# recalculate the probabilities (smallscore), non-triangle ext edge counted too (missing edge score = 0.001), then score by cohesive factor (small+largescore)
#	if ($mode eq "recalcprob_nontriext_cohfact") {
#		
#		my $MISSING_EDGE_PROB = 0.001;
#		if (defined $only_do_compsize && $only_do_compsize eq "comp3") {
#			# find triangles
#			foreach my $tri (keys %{$triangles_ref}) {
#				(my $prot1, my $prot2, my $prot3) = split(/\|/, $tri);
#				my $edge1 = "$prot1|$prot2";
#				my $edge2 = "$prot1|$prot3";
#				my $edge3 = "$prot2|$prot3";
#				# see if any triangles are incident on edge1
#				my $edge1score = $$edges_ref{$edge1};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb})) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge1score < 0) {
#							$edge1score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge1score < 0) {
#							$edge1score = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb})) { next; }
#					my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#					my $ext_edge1_score = $$edges_ref{$ext_edge1};
#					if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#					$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
#					if ($edge1score < 0) {
#						$edge1score = 0;
#					}
#				}
#				# see if any triangles are incident on edge2
#				my $edge2score = $$edges_ref{$edge2};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot3}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb})) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge2score < 0) {
#							$edge2score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge2score < 0) {
#							$edge2score = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
#					if ($nb eq $prot1 || $nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb})) {  next; }
#					my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#					if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#					my $ext_edge1_score = $$edges_ref{$ext_edge1};
#					if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#					$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
#					if ($edge2score < 0) {
#						$edge2score = 0;
#					}
#				}
#				# see if any triangles are incident on edge3			
#				my $edge3score = $$edges_ref{$edge3};
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb}) { next; }					
#					if (defined $$neighbours_ref{$prot3}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb})) {
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge3score < 0) {
#							$edge3score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge3score < 0) {
#							$edge3score = 0;
#						}
#					}
#				}		
#				foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
#					if ($nb eq $prot1 || $nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb})) { next; }
#					my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#					if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#					my $ext_edge1_score = $$edges_ref{$ext_edge1};
#					if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#					$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
#					if ($edge3score < 0) {
#						$edge3score = 0;
#					}
#				}
#				# sum up external edges			
#				my $extscores = 0;
#				foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb}) { next; }
#					my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#					if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_largecomp_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot1}{$nb}) { next; }
#					my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#					if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot2}}) {
#					if ($nb eq $prot1 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb}) { next; }
#					my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_largecomp_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1 || $nb eq $prot3) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot2}{$nb}) { next; }
#					my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot3}}) {
#					if ($nb eq $prot1 || $nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb}) { next; }
#					my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#					if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_largecomp_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
#					if ($nb eq $prot1 || $nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI3{$prot3}{$nb}) { next; }
#					my $ext_edge = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#					if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_ref{$ext_edge};
#				}
#				my $triscore = (($edge1score + $edge2score + $edge3score) / ($edge1score + $edge2score + $edge3score + $extscores)) * (($edge1score + $edge2score + $edge3score)/3);;
#				if (defined $ADDED_TRIANGLES{$tri}) { 	$triscore *= $addtri_penalty; }
#				if ($triscore >= $CLUS_SCORE_THR) { 
#					$$clus3_ref{$tri}{SCORE} = $triscore;
#					$$clus3_ref{$tri}{E}{$prot1} = 1;
#					$$clus3_ref{$tri}{E}{$prot2} = 1;
#					$$clus3_ref{$tri}{E}{$prot3} = 1;
#					$$clus_ref{$tri}{SCORE} = $triscore;
#					$$clus_ref{$tri}{E}{$prot1} = 1;
#					$$clus_ref{$tri}{E}{$prot2} = 1;
#					$$clus_ref{$tri}{E}{$prot3} = 1;
#				}
#			}
#		}
#				
#				
#		if (defined $only_do_compsize && $only_do_compsize eq "comp2") {
#			# find edges
#			foreach my $edge (keys %{$edges_ref}) {
#				(my $prot1, my $prot2) = split(/\|/, $edge);
#				my $newscore = $$edges_ref{$edge};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot1}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot2}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot2}{$nb})) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $ext_edge2_score;
#						if ($newscore < 0) {
#							$newscore = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($newscore < 0) {
#							$newscore = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot2}{$nb}) { next; }
#					if (defined $$neighbours_ref{$prot1}{$nb} && !($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot1}{$nb})) { next; }
#					my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#					my $ext_edge1_score = $$edges_ref{$ext_edge1};
#					if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#					$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
#					if ($newscore < 0) {
#						$newscore = 0;
#					}
#				}			
#				# sum up external edges			
#				my $extscores = 0;
#				foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot1}}) {
#					if ($nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot1}{$nb}) { next; }
#					my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#					if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_largecomp_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot1}{$nb}) { next; }
#					my $ext_edge = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#					if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_largecomp_ref{$prot2}}) {
#					if ($nb eq $prot1) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot2}{$nb}) { next; }
#					my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_largecomp_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_largecomp_ref{$ext_edge};
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1) { next; }
#					if ($USE_PPI==1 && ($USE_PPI_MODE==1 || $USE_PPI_MODE==2 || $USE_PPI_MODE==3) && !defined $NEIGHBOURS_PPI2{$prot2}{$nb}) { next; }
#					my $ext_edge = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#					if ($$edges_ref{$ext_edge} == 0) { die "ext_edge $ext_edge"; }
#					$extscores += $$edges_ref{$ext_edge};
#				}
#				my $edgescore = ($newscore / ($newscore + $extscores)) * $newscore;
#				if ($edgescore >= $CLUS_SCORE_THR) {
#					$$clus2_ref{$edge}{SCORE} = $edgescore;
#					$$clus2_ref{$edge}{E}{$prot1} = 1;
#					$$clus2_ref{$edge}{E}{$prot2} = 1;
#					$$clus_ref{$edge}{SCORE} = $edgescore;
#					$$clus_ref{$edge}{E}{$prot1} = 1;
#					$$clus_ref{$edge}{E}{$prot2} = 1;
#				}
#			}
#		}
#	} # end recalcprob_nontriext_cohfact
#	
#	
#	
#	
#	# just use the raw scores
#	elsif ($mode eq "rawscore") {
#		
#		if (defined $only_do_compsize && $only_do_compsize eq "comp3") {
#			# find triangles
#			foreach my $tri (keys %{$triangles_ref}) {
#				(my $prot1, my $prot2, my $prot3) = split(/\|/, $tri);
#				my $edge1 = "$prot1|$prot2";
#				my $edge2 = "$prot1|$prot3";
#				my $edge3 = "$prot2|$prot3";
#				my $edge1score = $$edges_ref{$edge1};
#				my $edge2score = $$edges_ref{$edge2};
#				my $edge3score = $$edges_ref{$edge3};
#				my $triscore = ($edge1score + $edge2score + $edge3score)/3;
#				if (defined $ADDED_TRIANGLES{$tri}) { 	$triscore *= $addtri_penalty; }
#				if ($triscore >= $CLUS_SCORE_THR) { 
#					$$clus3_ref{$tri}{SCORE} = $triscore;
#					$$clus3_ref{$tri}{E}{$prot1} = 1;
#					$$clus3_ref{$tri}{E}{$prot2} = 1;
#					$$clus3_ref{$tri}{E}{$prot3} = 1;
#					$$clus_ref{$tri}{SCORE} = $triscore;
#					$$clus_ref{$tri}{E}{$prot1} = 1;
#					$$clus_ref{$tri}{E}{$prot2} = 1;
#					$$clus_ref{$tri}{E}{$prot3} = 1;
#				}
#			}
#		}
#				
#				
#		if (defined $only_do_compsize && $only_do_compsize eq "comp2") {
#			# find edges
#			foreach my $edge (keys %{$edges_ref}) {
#				if ($USE_PPI>0 && !defined $EDGES_PPI2{$edge}) { next; }
#				(my $prot1, my $prot2) = split(/\|/, $edge);
#				my $newscore = $$edges_ref{$edge};
#				if ($newscore >= $CLUS_SCORE_THR) {
#					$$clus2_ref{$edge}{SCORE} = $newscore;
#					$$clus2_ref{$edge}{E}{$prot1} = 1;
#					$$clus2_ref{$edge}{E}{$prot2} = 1;
#					$$clus_ref{$edge}{SCORE} = $newscore;
#					$$clus_ref{$edge}{E}{$prot1} = 1;
#					$$clus_ref{$edge}{E}{$prot2} = 1;
#				}
#			}
#		}
#	} # end rawscore
#	
#	
#	
#	# recalculate the probabilities (smallscore), non-triangle ext edge counted too (missing edge score = 0.001)
#	elsif ($mode eq "recalcprob_nontriext") {
#		
#		my $MISSING_EDGE_PROB = 0.001;
#		if (defined $only_do_compsize && $only_do_compsize eq "comp3") {
#			# find triangles
#			foreach my $tri (keys %{$triangles_ref}) {
#				(my $prot1, my $prot2, my $prot3) = split(/\|/, $tri);
#				my $edge1 = "$prot1|$prot2";
#				my $edge2 = "$prot1|$prot3";
#				my $edge3 = "$prot2|$prot3";
##				if ($USE_PPI==1 && (!defined $EDGES_PPI3{$edge1} || !defined $EDGES_PPI3{$edge2} || !defined $EDGES_PPI3{$edge3})) { next; }
#				my $edge1score = $$edges_ref{$edge1};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if (defined $$neighbours_ref{$prot2}{$nb}) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge1score < 0) {
#							$edge1score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge1score < 0) {
#							$edge1score = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1 || $nb eq $prot3) { next; }
#					if (!defined $$neighbours_ref{$prot1}{$nb}) { 
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge1score -= $$edges_ref{$edge1} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge1score < 0) {
#							$edge1score = 0;
#						}
#					}
#				}
#				my $edge2score = $$edges_ref{$edge2};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2 || $nb eq $prot3) { next; }
#					if (defined $$neighbours_ref{$prot3}{$nb}) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge2score < 0) {
#							$edge2score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge2score < 0) {
#							$edge2score = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
#					if ($nb eq $prot1) { next; }
#					if ($nb eq $prot2) { next; }
#					if (!defined $$neighbours_ref{$prot1}{$nb}) { 
#						my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge2score -= $$edges_ref{$edge2} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge2score < 0) {
#							$edge2score = 0;
#						}
#					}
#				}
#				
#				my $edge3score = $$edges_ref{$edge3};
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1) { next; }
#					if ($nb eq $prot3) { next; }
#					if (defined $$neighbours_ref{$prot3}{$nb}) {
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						my $ext_edge2 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $ext_edge2_score;
#						if ($edge3score < 0) {
#							$edge3score = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge3score < 0) {
#							$edge3score = 0;
#						}
#					}
#				}		
#				foreach my $nb (keys %{$$neighbours_ref{$prot3}}) {
#					if ($nb eq $prot1) { next; }
#					if ($nb eq $prot2) { next; }
#					if (!defined $$neighbours_ref{$prot2}{$nb}) { 
#						my $ext_edge1 = $prot3 lt $nb ? "$prot3|$nb" : "$nb|$prot3";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$edge3score -= $$edges_ref{$edge3} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($edge3score < 0) {
#							$edge3score = 0;
#						}
#					}
#				}
#				my $triscore = ($edge1score + $edge2score + $edge3score)/3;
#				if (defined $ADDED_TRIANGLES{$tri}) { 	$triscore *= $addtri_penalty; }
#				if ($triscore >= $CLUS_SCORE_THR) { 
#					$$clus3_ref{$tri}{SCORE} = $triscore;
#					$$clus3_ref{$tri}{E}{$prot1} = 1;
#					$$clus3_ref{$tri}{E}{$prot2} = 1;
#					$$clus3_ref{$tri}{E}{$prot3} = 1;
#					$$clus_ref{$tri}{SCORE} = $triscore;
#					$$clus_ref{$tri}{E}{$prot1} = 1;
#					$$clus_ref{$tri}{E}{$prot2} = 1;
#					$$clus_ref{$tri}{E}{$prot3} = 1;
#				}
#			}
#		}
#				
#				
#		if (defined $only_do_compsize && $only_do_compsize eq "comp2") {
#			# find edges
#			foreach my $edge (keys %{$edges_ref}) {
#				if ($USE_PPI>0 && !defined $EDGES_PPI2{$edge}) { next; }
#				(my $prot1, my $prot2) = split(/\|/, $edge);
#				my $newscore = $$edges_ref{$edge};
#				foreach my $nb (keys %{$$neighbours_ref{$prot1}}) {
#					if ($nb eq $prot2) { next; }
#					if (defined $$neighbours_ref{$prot2}{$nb}) {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						my $ext_edge2 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						if ($$edges_ref{$ext_edge2} == 0) { die "ext_edge $ext_edge2"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						my $ext_edge2_score = $$edges_ref{$ext_edge2};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						if ($ext_edge2_score < $MISSING_EDGE_PROB) { $ext_edge2_score = $MISSING_EDGE_PROB; }
#						$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $ext_edge2_score;
#						if ($newscore < 0) {
#							$newscore = 0;
#						}
#					}
#					else {
#						my $ext_edge1 = $prot1 lt $nb ? "$prot1|$nb" : "$nb|$prot1";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($newscore < 0) {
#							$newscore = 0;
#						}
#					}
#				}
#				foreach my $nb (keys %{$$neighbours_ref{$prot2}}) {
#					if ($nb eq $prot1) { next; }
#					if (!defined $$neighbours_ref{$prot1}{$nb}) { 
#						my $ext_edge1 = $prot2 lt $nb ? "$prot2|$nb" : "$nb|$prot2";
#						if ($$edges_ref{$ext_edge1} == 0) { die "ext_edge $ext_edge1"; }
#						my $ext_edge1_score = $$edges_ref{$ext_edge1};
#						if ($ext_edge1_score < $MISSING_EDGE_PROB) { $ext_edge1_score = $MISSING_EDGE_PROB; }
#						$newscore -= $$edges_ref{$edge} * $ext_edge1_score * $MISSING_EDGE_PROB;
#						if ($newscore < 0) {
#							$newscore = 0;
#						}
#					}
#				}			
#				if ($newscore >= $CLUS_SCORE_THR) {
#					$$clus2_ref{$edge}{SCORE} = $newscore;
#					$$clus2_ref{$edge}{E}{$prot1} = 1;
#					$$clus2_ref{$edge}{E}{$prot2} = 1;
#					$$clus_ref{$edge}{SCORE} = $newscore;
#					$$clus_ref{$edge}{E}{$prot1} = 1;
#					$$clus_ref{$edge}{E}{$prot2} = 1;
#				}
#			}
#		}
#	} # end recalcprob_nontriext
#	
#	
#	
#	else { die "mode $mode"; }
#		
#	print "num cluster2 = ".(scalar keys %$clus2_ref).", num cluster3 = ".(scalar keys %$clus3_ref).", num clusters all = ".(scalar keys %$clus_ref)."\n";
#
#}


