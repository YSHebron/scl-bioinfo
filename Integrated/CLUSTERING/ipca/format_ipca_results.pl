use strict;
use IO::Handle;


# ARGV[0]: Input data file
# ARGV[1]: Scored edges file
# ARGV[2]: Cluster scoring method (1 = weighted density, 2 = density)
# ARGV[3]: Output file
# Eg. perl format_ipca_results "clusters_ipca_keep10k_iter0.txt"
#															"scored_edges keep10k iter0.txt"
#															1
#															"clusters ipcaformatted keep10k iter0.txt"

open (INPUTFILE, $ARGV[0]) or die $!;
open (EDGESFILE, $ARGV[1]) or die $!;
open (OUTPUTFILE, ">$ARGV[3]") or die $!;
OUTPUTFILE->autoflush(1);

my $scoring_method = $ARGV[2] + 0;


# read the scored edges file
my %scorededges;
foreach my $line (<EDGESFILE>) {
	chomp $line;
	(my $id_a, my $id_b, my $score) = split(' ', $line);
	my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
	$scorededges{$key} = $score + 0;
}
	
# read complexes
my %clusters;
my $curr_clus;
foreach my $line (<INPUTFILE>) {
	chomp $line;
	my @tokens = split (" ", $line);
	
	if ($tokens[0] eq "Complex") {
		$curr_clus = $tokens[1];
	}
	else {
		$clusters{$curr_clus}{ELEMS}{$tokens[0]} = 1;
	}
}


ScoreClusters(\%clusters, \%scorededges, $scoring_method);
	
# print
foreach my $clus (sort {$clusters{$b}{SCORE} <=> $clusters{$a}{SCORE}} keys %clusters) {
	my @prots = sort keys %{$clusters{$clus}{ELEMS}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
#	if ($size < 4) { next; }
	print OUTPUTFILE "C$clus(".$size."_$clusters{$clus}{SCORE}): $protsstring\n";
}	






# -------------------------------------------------------------

sub ScoreClusters($$$) {
	my $clustersref = $_[0];
	my $edgesref = $_[1];
	my $method = $_[2]; # 1 = weighted density, 2 = density
	
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
}







