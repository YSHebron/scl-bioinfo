use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use POSIX;

use Getopt::Std;




# inputs:
# -i <input clusters file>
# -s <topological scores (from iterative-adjustCD) file>
# -h <hubs file>
# -t <hub add threshold>
# -o <output clusters filename>

my %argopts;
if (! getopts('i:s:t:h:o:', \%argopts)) { die "invalid arguments"; }

my $inputfilename = $argopts{'i'};
my $hubsfilename = $argopts{'h'};
my $hubaddthresh = $argopts{'t'}+0;


# $clusters{$c}{SCORE} = score, $clusters{$c}{ELEMENTS}{$pid} = 1
my %clusters = ();
ReadClusters($inputfilename, \%clusters);



my %hub_prots = (); # $hub_prots{$prot} = 1 if $prot is a hub
my %hub_edges = (); # $hub_edges{$edge} = score
my %hub_nbs = (); # $hub_nbs{$prot}{$hub} = score
open (HUBSFILE, $hubsfilename) || die $!;
my $in_edges = 0;
foreach my $line (<HUBSFILE>) {
	chomp $line;
	if ($line eq "hubs") { 
		$in_edges = 0;
	}
	elsif ($line eq "edges") {
		$in_edges = 1;
	}
	elsif ($in_edges==0) {
		(my $prot, my $deg) = split(/\s/, $line);
		$hub_prots{$prot} = 1;
	}
	else {
		(my $edge, my $score) = split(/\s/, $line);
		(my $prot1, my $prot2) = split(/\|/, $edge);
		if ($prot1 eq $prot2) { next; }
		$hub_edges{$edge} = $score;
		if (defined $hub_prots{$prot1}) {
			$hub_nbs{$prot2}{$prot1} = $score;
		}
		if (defined $hub_prots{$prot2}) {
			$hub_nbs{$prot1}{$prot2} = $score;
		}
		if (!defined $hub_prots{$prot1} && !defined $hub_prots{$prot2}) {
			die "hub edge line $line, prots $prot1 $prot2, but neither prots are hubs";
		}
	}
}
print "Num hub prots = ".(scalar keys %hub_prots).", num hub edges = ".(scalar keys %hub_edges).", num hub nbs = ".(scalar keys %hub_nbs)."\n";


# read topological scores
my %topo_scores = (); # $topo_scores{$prot1|$prot2} = topo score
open (TOPOFILE, $argopts{'s'}) || die $!;
foreach my $line (<TOPOFILE>) {
	chomp $line;
	(my $prot1, my $prot2, my $score) = split(/\s/, $line);
	if ($prot1 eq $prot2) { next; }
	my $edge = $prot1 lt $prot2 ? "$prot1|$prot2" : "$prot2|$prot1";
	$topo_scores{$edge} = $score;
}

my %clusters_added_hubs = ();
my %hubs_added = ();
foreach my $clus (keys %clusters) {
	my %connected_hubs = (); # $connected_hubs{$hub} = average PPI score connecting hub to the cluster
	my %connected_hubs_topo_score = (); # $connected_hubs_topo_score{$hub} = average topological score connecting hub to the cluster
#	print "Clus $clus, size ".(scalar keys %{$clusters{$clus}{ELEMENTS}}).", connected to hubs:\n";
	foreach my $prot (keys %{$clusters{$clus}{ELEMENTS}}) {
		if (defined $hub_nbs{$prot}) {
			foreach my $hub (keys %{$hub_nbs{$prot}}) {
				if (defined $clusters{$clus}{ELEMENTS}{$hub}) { next; }
				if ($prot eq $hub) { die; }
				if (!defined $hub_nbs{$prot}{$hub}) { die; }
				$connected_hubs{$hub} += $hub_nbs{$prot}{$hub};
				my $edg = $prot lt $hub ? "$prot|$hub" : "$hub|$prot";
				$connected_hubs_topo_score{$hub} += $topo_scores{$edg};
#				print "hub $hub, cluster elem $prot, topo score $topo_scores{$edg}, ppi score $hub_nbs{$prot}{$hub}\n";
			}
		}
	}
	foreach my $hub (keys %connected_hubs) {
		$connected_hubs{$hub} /= scalar keys %{$clusters{$clus}{ELEMENTS}};
		$connected_hubs_topo_score{$hub} /= scalar keys %{$clusters{$clus}{ELEMENTS}};
	}
	foreach my $hub (sort {$connected_hubs{$b} <=> $connected_hubs{$a}} keys %connected_hubs) {
		if ($connected_hubs_topo_score{$hub} >= $hubaddthresh) {
			$clusters_added_hubs{$clus} = 1;
			$hubs_added{$clus} = 1;
			my $old_num_prots = scalar keys %{$clusters{$clus}{ELEMENTS}};
			my $old_total_edges = $old_num_prots * ($old_num_prots-1) * .5;
			my $new_total_edges = ($old_num_prots+1) * ($old_num_prots) * .5;
			my $new_score = ($clusters{$clus}{SCORE} * $old_total_edges + $connected_hubs{$hub} * $old_num_prots) / $new_total_edges;
			$clusters{$clus}{ELEMENTS}{$hub} = 1;
			$clusters{$clus}{SCORE} = $new_score;
#			print "$hub\tTopo connection $connected_hubs_topo_score{$hub}\tADDED, PPI connection $connected_hubs{$hub}, new score $clusters{$clus}{SCORE}\n";
		}
		else {
#			print "$hub\tTopo connection $connected_hubs_topo_score{$hub}\n";
		}
	}
	
}
print "Num clusters with hubs added = ".(scalar keys %clusters_added_hubs)."\n";
print "Num hubs added = ".(scalar keys %hubs_added)."\n";
	

# print
open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $clus (sort {$clusters{$b}{SCORE} <=> $clusters{$a}{SCORE}} keys %clusters) {
	my @prots = sort keys %{$clusters{$clus}{ELEMENTS}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
#	if ($size < 4) { next; }
	print OUTPUTFILE "C$clus(".$size."_$clusters{$clus}{SCORE}): $protsstring\n";
}	







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









