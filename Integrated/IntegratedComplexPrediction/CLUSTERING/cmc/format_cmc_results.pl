use strict;
use IO::Handle;


# ARGV[0]: Input data file
# ARGV[1]: Scored edges file
# ARGV[2]: Cluster scoring method (1 = weighted density, 2 = density, 3 = coherence factor)
# ARGV[3]: Output file

open (EDGESFILE, $ARGV[1]) or die $!;
open (OUTPUTFILE, ">$ARGV[3]") or die $!;
OUTPUTFILE->autoflush(1);

my $scoring_method = $ARGV[2] + 0;


# read the scored edges file
my %scorededges;
my %neighbours; # only used if method==3
foreach my $line (<EDGESFILE>) {
	chomp $line;
	(my $id_a, my $id_b, my $score) = split(' ', $line);
	my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
	$scorededges{$key} = $score + 0;
	if ($scoring_method==3) {
		$neighbours{$id_a}{$id_b} = $score+0;
		$neighbours{$id_b}{$id_a} = $score+0;
	}
}
	
	
	
# read clusters
my %clusters = ();
ReadClusters($ARGV[0], \%clusters);
	
	


ScoreClusters(\%clusters, \%scorededges, $scoring_method);
	
# print
foreach my $clus (sort {$clusters{$b}{SCORE} <=> $clusters{$a}{SCORE}} keys %clusters) {
	my @prots = sort keys %{$clusters{$clus}{ELEMS}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
	if ($size < 2) { next; }
	print OUTPUTFILE "C$clus(".$size."_$clusters{$clus}{SCORE}): $protsstring\n";
}	






# -------------------------------------------------------------



# $clusters{$c}{ELEMS}{$pid} = 1
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
		shift @toks;
		foreach (@toks) {
			$$clusters_ref{$clus}{ELEMS}{$_} = 1;
		}
	}		
	close CLUSTERS_FILE;	
}





sub ScoreClusters($$$) {
	my $clustersref = $_[0];
	my $edgesref = $_[1];
	my $method = $_[2]; # 1 = weighted density, 2 = density, 3 = coherence factor
	
	if ($method == 1 || $method == 2) {
		foreach my $clus (keys %$clustersref) {
			my $totalweight = 0;
			my @prots = keys %{$$clustersref{$clus}{ELEMS}};
			for (my $idx1 = 0; $idx1 < scalar @prots; $idx1++) {
				for (my $idx2 = $idx1+1; $idx2 < scalar @prots; $idx2++) {
					my $id_a = $prots[$idx1];
					my $id_b = $prots[$idx2];
					my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
					if ($method==1) {
						if (defined $$edgesref{$key}) {
							$totalweight += $$edgesref{$key};
						}
					}
					elsif ($method==2) {
						if (defined $$edgesref{$key} && $$edgesref{$key} > 0) {
							$totalweight++;
						}
					}
				}
			}
			my $score = $totalweight*2 / ((scalar @prots) * ((scalar @prots) -1));
			$$clustersref{$clus}{SCORE} = $score;
		}
	} # end if method == 1 or 2
	
	elsif ($method == 3) {
		foreach my $clus (keys %$clustersref) {
			my $intweight = 0;
			my @prots = keys %{$$clustersref{$clus}{ELEMS}};
			for (my $idx1 = 0; $idx1 < scalar @prots; $idx1++) {
				for (my $idx2 = $idx1+1; $idx2 < scalar @prots; $idx2++) {
					my $id_a = $prots[$idx1];
					my $id_b = $prots[$idx2];
					my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
					if (defined $$edgesref{$key}) {
						$intweight += $$edgesref{$key};
					}
				}
			}
			# sum external weights
			my $extweight = 0;
			foreach my $prot (@prots) {
				foreach my $nb (keys %{$neighbours{$prot}}) {
					if (defined $$clustersref{$clus}{ELEMS}{$nb}) { next; }
					$extweight += $neighbours{$prot}{$nb};
				}
			}
			$extweight += $intweight;
			my $density = $intweight * 2 / ((scalar @prots) * ((scalar @prots) -1));
			my $score = $intweight / $extweight * $density;
			$$clustersref{$clus}{SCORE} = $score;			
		} 
	} # end if method == 3
}








