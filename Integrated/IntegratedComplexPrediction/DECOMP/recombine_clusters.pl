use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use POSIX;

use Getopt::Std;




# inputs:
# -i <input clusters basefilename>
# -n <number of files>
# -p <ppi edges file, for scoring clusters>
# -o <output clusters filename>

my %argopts;
if (! getopts('i:n:p:o:', \%argopts)) { die "invalid arguments"; }

my $inputbasefilename = $argopts{'i'};
my $outputfilename = $argopts{'o'};
my $numfiles = $argopts{'n'}+0;


# read ppi data
my %ppi_edges = (); # $ppi_edges{$key} = $score
open (PPIFILE, "$argopts{'p'}") || die $!;
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
	$ppi_edges{$key} = $score;
}
	
	

# $combined_clusters{$c}{SCORE} = score, $clusters{$c}{ELEMENTS}{$pid} = 1
my %combined_clusters = ();
for (my $idx=0; $idx<$numfiles; $idx++) {
	my $inputfile = $inputbasefilename . "_" . $idx . ".txt";
	open (INPUTFILE, "$inputfile") || die "$inputfile $!";
	my %clusters = ();
	ReadClusters($inputfile, \%clusters);
	
	# remove duplicates
	foreach my $newclus (keys %clusters) {
		foreach my $clus (keys %combined_clusters) {
			if (EqualSets($clusters{$newclus}{ELEMENTS}, $combined_clusters{$clus}{ELEMENTS})==1) {
				delete $clusters{$newclus};
#				print "Cluster $newclus from $inputfile is same as $clus\n";
				last;
			}
		}
	}
	foreach my $newclus (keys %clusters) {
		my $newclusname = "$newclus"."_$idx";
		$combined_clusters{$newclusname}{SCORE} = $clusters{$newclus}{SCORE};
		foreach my $prot (keys %{$clusters{$newclus}{ELEMENTS}}) {
			$combined_clusters{$newclusname}{ELEMENTS}{$prot} = 1;
		}
	}
	
	close INPUTFILE;
}



ScoreClusters(\%combined_clusters, \%ppi_edges, 1);


# print
open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
foreach my $clus (sort {$combined_clusters{$b}{SCORE} <=> $combined_clusters{$a}{SCORE}} keys %combined_clusters) {
	my @prots = sort keys %{$combined_clusters{$clus}{ELEMENTS}};
	my $protsstring = join (" ", @prots);
	my $size = scalar @prots;
#	if ($size < 4) { next; }
	print OUTPUTFILE "C$clus(".$size."_$combined_clusters{$clus}{SCORE}): $protsstring\n";
}	



sub EqualSets {
	my $clus1_ref = $_[0];
	my $clus2_ref = $_[1];
	if  (scalar keys %{$clus1_ref} != scalar keys %{$clus2_ref}) {
		return 0;
	}
	foreach my $e1 (keys %{$clus1_ref}) {
		if (!defined $$clus2_ref{$e1}) {
			return 0;
		}
	}
	return 1;
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

sub ScoreClusters($$$) {
	my $clustersref = $_[0];
	my $edgesref = $_[1];
	my $method = $_[2]; # 1 = weighted density, 2 = density
	
	foreach my $clus (keys %$clustersref) {
		my $totalweight = 0;
		my @prots = keys %{$$clustersref{$clus}{ELEMENTS}};
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













