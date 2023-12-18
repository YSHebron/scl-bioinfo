use strict;
use IO::Handle;
use List::Util 'shuffle';
use List::Util 'max';
use POSIX;

use Getopt::Std;




# inputs:
# -i <input data file, OR basefilename (see -x)>
# -k <k = number of top edges to keep, or "all">
# -s <scaling mode: "norm" or "ln" or "noweight" or "none" or "topall1" or "scale" or "stretch">
# -t <if scaling mode is "scale", scale by this factor. if "stretch", scale minimum score by this and "stretch" the rest>
# -h <hub removal: "none" or "hub" or "edge">
# -g <hub threshold, default 30>
# -n <Number of iterations. If defined, will iterate over the input file name, eg if input is "swc xval.8.2 scored_edges", will process "swc xval.8.2 scored_edges iter0.txt", "swc xval.8.2 scored_edges iter1.txt", etc.>
# -o <output filename OR basefilename (if -x option used, and basefilename is "swc xval.8.2 scored_edges 10k", will output "swc xval.8.2 scored_edges 10k iter0.txt", "swc xval.8.2 scored_edges 10k iter1.txt", etc>

my %argopts;
if (! getopts('i:k:s:n:h:g:t:o:', \%argopts)) { die "invalid arguments"; }

my $num_iters = 0;
if (defined $argopts{'n'}) {
	$num_iters = $argopts{'n'} + 0;
	#if ($num_iters == 0) { die; }
}

if (!defined $argopts{'s'}) { die; }
my $mode = $argopts{'s'}; 

my $scalefactor = 1;
if (defined $argopts{'t'}) {
	$scalefactor = $argopts{'t'}+0;
}

if (!defined $argopts{'k'}) { die; }
my $num_to_keep = $argopts{'k'};

my $hubremovalmode = "none";
my $hubthresh = 30;
if (defined $argopts{'h'}) {
	$hubremovalmode = $argopts{'h'};
	if ($hubremovalmode ne "none" && $hubremovalmode ne "hub" && $hubremovalmode ne "edge") { die; }
	if (defined $argopts{'g'}) {
		$hubthresh = $argopts{'g'} + 0;
	}
}




if ($num_iters==0) {
	open (INPUTFILE, $argopts{'i'}) || die $!;
	open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
	
	my %data = ();
	foreach my $line (<INPUTFILE>) {
		chomp $line;
		(my $id_a, my $id_b, my $score) = split(/\s/, $line);
		my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
		$data{$key} = $score;
	}
	
	if ($mode eq "norm") {
		# cap at maxscore
		my $maxscore = 100;
		foreach my $key (keys %data) {
			my $score = $data{$key};
			if ($score > $maxscore) {
				$data{$key} = $maxscore;
			}
		}
		# normalize
		my $highest = -9999;
		foreach my $key (keys %data) {
			if ($data{$key} > $highest) {
				$highest = $data{$key};
			}
		}
		foreach my $key (keys %data) {
			$data{$key} = $data{$key} / $highest;
		}
		# create data to be output
		my %outputdata = (); # $outputdata{$key} = score
		if ($hubremovalmode eq "hub") {
#			my $count = 0;
#			my %protdegree = (); # $protdegree{$prot} = degree of $prot
#			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#				$outputdata{$key} = $data{$key};
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				$protdegree{$id_a}++;
#				$protdegree{$id_b}++;
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}
#			my %hubs = ();
#			foreach my $prot (keys %protdegree) {
#				if ($protdegree{$prot} > $hubthresh) {
#					$hubs{$prot} =1 ;
#				}
#			}
#			print "Num hubs found = ".(scalar keys %hubs)."\n";			
		} # end $hubremovalmode eq "hub"
		elsif ($hubremovalmode eq "edge") {
			my $count = 0;
			my %protdegree = (); # $protdegree{$prot} = degree of $prot
			my $num_hub_edges = 0;
			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				if ($protdegree{$id_a} >= $hubthresh || $protdegree{$id_b} >= $hubthresh) { 
					$num_hub_edges++;
					next; 
				}
				$protdegree{$id_a}++;
				$protdegree{$id_b}++;
				$outputdata{$key} = $data{$key};
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}
			print "Num hub edges found = $num_hub_edges\n";
		} # end $hubremovalmode eq "edge"
		else {
			my $count = 0;
			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
				$outputdata{$key} = $data{$key};
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}			
		}		
		# print
		foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %outputdata) {
			(my $id_a, my $id_b) = split(/\|/, $key);
			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
		}
	} # end if ($mode eq "norm")
	
	
	elsif ($mode eq "noweight") {
		my $count = 0;
		my @shuffled_array = keys %data;		
		fisher_yates_shuffle( \@shuffled_array );
		foreach my $key (@shuffled_array) {
			(my $id_a, my $id_b) = split(/\|/, $key);
			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
			$count++;
			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
				last;
			}
		}
	}
	
#	elsif ($mode eq "max1") {
#		# print
#		my $count = 0;
#		foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			my $score = $data{$key};
#			if ($score > 1) { $score = 1; }
#			print OUTPUTFILE ("$id_a $id_b $score\n");
#			$count++;
#			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#				last;
#			}
#		}
#	}
#	
#	elsif ($mode eq "lnnorm") {
#		# cap at maxscore
#		my $maxscore = 100;
#		foreach my $key (keys %data) {
#			my $score = $data{$key};
#			if ($score > $maxscore) {
#				$data{$key} = $maxscore;
#			}
#		}
#		# apply ln
#		foreach my $key (keys %data) {
#			my $score = $data{$key};
#			$score = log($score+1);
#			$data{$key} = $score;
#		}
#		# normalize
#		my $highest = -9999;
#		foreach my $key (keys %data) {
#			if ($data{$key} > $highest) {
#				$highest = $data{$key};
#			}
#		}
#		foreach my $key (keys %data) {
#			$data{$key} = $data{$key} / $highest;
#		}
#		# print
#		my $count = 0;
#		foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#			$count++;
#			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#				last;
#			}
#		}
#	}
#	elsif ($mode eq "ln") {
#		# cap scores 
#		my $maxscore = 100;
#		foreach my $key (keys %data) {
#			my $score = $data{$key};
#			if ($score > $maxscore) {
#				$data{$key} = $maxscore;
#			}
#		}
#		# scale
#		foreach my $key (keys %data) {
#			my $score = $data{$key};
#			if ($score >= 0.9) {
#				$data{$key} = 0.9 + log($score + 0.1)/log($maxscore+0.1) * 0.1;
#			}
#		}
#		# print
#		my $count = 0;
#		foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#			$count++;
#			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#				last;
#			}
#		}		
#	}
			
	
	else { die; }
	close INPUTFILE;
	close OUTPUTFILE; 
}





# do multiple iterations
else {
	
	my $inputbasename = $argopts{'i'};
	my $outputbasename = $argopts{'o'};
	
	for (my $iter = 0; $iter < $num_iters; $iter++) {
		open (INPUTFILE, ("$inputbasename iter".$iter.".txt")) || die $!;
		open (OUTPUTFILE, (">$outputbasename iter".$iter.".txt")) || die $!;
	
	
		my %data = ();
		foreach my $line (<INPUTFILE>) {
			chomp $line;
			(my $id_a, my $id_b, my $score) = split(/\s/, $line);
			my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
			$data{$key} = $score;
		}
		
		if ($mode eq "none") {
			# create data to be output
			my %outputdata = (); # $outputdata{$key} = score
			my $count = 0;
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
				$outputdata{$key} = $data{$key};
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}			
			# print
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %outputdata) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
			}
		}
		
		elsif ($mode eq "scale") {
			# create data to be output
			my %outputdata = (); # $outputdata{$key} = score
			my $count = 0;
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
				$outputdata{$key} = $scalefactor * $data{$key};
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}			
			# print
			foreach my $key (sort {$outputdata{$b} <=> $outputdata{$a} || $a cmp $b} keys %outputdata) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b $outputdata{$key}\n");
			}
		}
		
		elsif ($mode eq "stretch") {
			# create data to be output
			my %outputdata = (); # $outputdata{$key} = score
			my $count = 0;
			my $minscore = 999;
			my $maxscore = -999;
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
				$outputdata{$key} = $data{$key};
				my $score = $data{$key};
				if ($score < $minscore) {
					$minscore = $score;
				}
				if ($score > $maxscore) {
					$maxscore = $score;
				}
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}	
			my $newminscore = $minscore * $scalefactor;	
			foreach my $key (keys %outputdata) {
				my $score = $outputdata{$key};
				my $newscore = ($score - $minscore) / ($maxscore - $minscore) * ($maxscore - $newminscore) + $newminscore;
				$outputdata{$key} = $newscore;
			}	
			# print
			foreach my $key (sort {$outputdata{$b} <=> $outputdata{$a} || $a cmp $b} keys %outputdata) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b $outputdata{$key}\n");
			}
		}
		
		elsif ($mode eq "topall1") {
			# create data to be output
			my %outputdata = (); # $outputdata{$key} = score
			my $count = 0;
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
				$outputdata{$key} = $data{$key};
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}			
			# print
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %outputdata) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b 1\n");
			}
		}
		
		elsif ($mode eq "norm") {
			# cap at maxscore
			my $maxscore = 100;
			foreach my $key (keys %data) {
				my $score = $data{$key};
				if ($score > $maxscore) {
					$data{$key} = $maxscore;
				}
			}
			# normalize
			my $highest = -9999;
			foreach my $key (keys %data) {
				if ($data{$key} > $highest) {
					$highest = $data{$key};
				}
			}
			foreach my $key (keys %data) {
				$data{$key} = $data{$key} / $highest;
			}
			# create data to be output
			my %outputdata = (); # $outputdata{$key} = score
			if ($hubremovalmode eq "hub") { 
				my $count = 0;
				foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
					$outputdata{$key} = $data{$key};
					$count++;
					if ($num_to_keep ne "all" && $count >= $num_to_keep) {
						last;
					}
				}			
				my %protdegree = (); # $protdegree{$prot} = degree of $prot
				my %hub_edges = ();
				my %hub_prots = ();
				foreach my $key (keys %outputdata) {
					(my $id_a, my $id_b) = split(/\|/, $key);
					$protdegree{$id_a}++;
					$protdegree{$id_b}++;
				}
				foreach my $key (keys %outputdata) {
					(my $id_a, my $id_b) = split(/\|/, $key);
					if ($protdegree{$id_a} >= $hubthresh) {
						$hub_prots{$id_a} = $protdegree{$id_a};
						$hub_edges{$key} = $outputdata{$key};
					}
					if ($protdegree{$id_b} >= $hubthresh) {
						$hub_prots{$id_b} = $protdegree{$id_b};
						$hub_edges{$key} = $outputdata{$key};
					}
				}
				foreach my $key (keys %hub_edges) {
					delete $outputdata{$key};
				}
				print "Num hubs found = ".(scalar keys %hub_prots).", num edges removed = ".(scalar keys %hub_edges)."\n";
				open (HUBSFILE, (">$outputbasename hubs iter".$iter.".txt")) || die $!;
				print HUBSFILE "hubs\n";
				foreach my $prot (keys %hub_prots) {
					print HUBSFILE "$prot\t$hub_prots{$prot}\n";
				}
				print HUBSFILE "edges\n";
				foreach my $key (keys %hub_edges) {
					(my $id_a, my $id_b) = split(/\|/, $key);
					print HUBSFILE "$id_a\t$id_b\t$hub_edges{$key}\n";
				}
				close HUBSFILE;
			}
			elsif ($hubremovalmode eq "edge") {
				my $count = 0;
				my %protdegree = (); # $protdegree{$prot} = degree of $prot
				my $num_hub_edges = 0;
				my %hub_edges = (); # just for debugging
				foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
					(my $id_a, my $id_b) = split(/\|/, $key);
					if ($protdegree{$id_a} >= $hubthresh || $protdegree{$id_b} >= $hubthresh) { 
						$num_hub_edges++;
						$hub_edges{$key} = $data{$key};
						next; 
					}
					$protdegree{$id_a}++;
					$protdegree{$id_b}++;
					$outputdata{$key} = $data{$key};
					$count++;
					if ($num_to_keep ne "all" && $count >= $num_to_keep) {
						last;
					}
				}
				print "Num hub edges found = $num_hub_edges\n";
#				foreach my $hubedge (sort {$hub_edges{$b} <=> $hub_edges{$a}} keys %hub_edges) {
#					print "$hubedge\t$hub_edges{$hubedge}\n";
#				}
			} # end $hubremovalmode eq "edge"
			else {
				my $count = 0;
				foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %data) {
					$outputdata{$key} = $data{$key};
					$count++;
					if ($num_to_keep ne "all" && $count >= $num_to_keep) {
						last;
					}
				}			
			}		
			
			
			# print
			foreach my $key (sort {$data{$b} <=> $data{$a} || $a cmp $b} keys %outputdata) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
			}
		} # end if ($mode eq "norm")
		

		
		elsif ($mode eq "noweight") {
			my $count = 0;
			my @shuffled_array = keys %data;		
			fisher_yates_shuffle( \@shuffled_array );
			foreach my $key (@shuffled_array) {
				(my $id_a, my $id_b) = split(/\|/, $key);
				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
				$count++;
				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
					last;
				}
			}
		}
		
#		elsif ($mode eq "max1") {
#			# print
#			my $count = 0;
#			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				my $score = $data{$key};
#				if ($score > 1) { $score = 1; }
#				print OUTPUTFILE ("$id_a $id_b $score\n");
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}			
#		}
#		
#		elsif ($mode eq "ln") {
#			# cap scores 
#			my $maxscore = 100;
#			foreach my $key (keys %data) {
#				my $score = $data{$key};
#				if ($score > $maxscore) {
#					$data{$key} = $maxscore;
#				}
#			}
#			# scale
#			foreach my $key (keys %data) {
#				my $score = $data{$key};
#				if ($score >= 0.9) {
#					$data{$key} = 0.9 + log($score + 0.1)/log($maxscore+0.1) * 0.1;
#				}
#			}
#			# print
#			my $count = 0;
#			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}
#		}
#				elsif ($mode eq "lnnorm") {
#			# cap at maxscore
#			my $maxscore = 100;
#			foreach my $key (keys %data) {
#				my $score = $data{$key};
#				if ($score > $maxscore) {
#					$data{$key} = $maxscore;
#				}
#			}
#			# apply ln
#			foreach my $key (keys %data) {
#				my $score = $data{$key};
#				$score = log($score+1);
#				$data{$key} = $score;
#			}
#			# normalize
#			my $highest = -9999;
#			foreach my $key (keys %data) {
#				if ($data{$key} > $highest) {
#					$highest = $data{$key};
#				}
#			}
#			foreach my $key (keys %data) {
#				$data{$key} = $data{$key} / $highest;
#			}
#			# print
#			my $count = 0;
#			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}
#		}
		
		else { die; }
		
	
		close INPUTFILE;
		close OUTPUTFILE; 
	}
	
}

#	
#if ($num_iters==0) {
#	open (INPUTFILE, $argopts{'i'}) || die $!;
#	open (OUTPUTFILE, ">$argopts{'o'}") || die $!;
#	
#	# Transform data to [0,1]
#	my %data = ();
#	my $maxscore = 100;
#	my $factor = 1;
#	my $highest = -9999;
#	foreach my $line (<INPUTFILE>) {
#		chomp $line;
#		(my $id_a, my $id_b, my $score) = split(' ', $line);
#		my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
#		my $score = $factor * $score;
#		if ($score > $maxscore) {
#			$score = $maxscore;
#		}
#		if ($doln==1) { 
#			$score = log($score+1); 
#		}
#		if ($score > $highest) {
#			$highest = $score;
#		}
#		$data{$key} = $score;
#	}
#	# normalize
#	foreach my $key (keys %data) {
#		$data{$key} = $data{$key} / $highest;
#	}
#		
#	if ($noweight==1) {		
#		my $count = 0;
#		my @shuffled_array = keys %data;		
#		fisher_yates_shuffle( \@shuffled_array );
#		foreach my $key (@shuffled_array) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#			$count++;
#			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#				last;
#			}
#		}
#	}
#	
#	else {	
#		# filter by number to keep
#		my $count = 0;
#		foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#			(my $id_a, my $id_b) = split(/\|/, $key);
#			print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#			$count++;
#			if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#				last;
#			}
#		}
#	}
#	close INPUTFILE;
#	close OUTPUTFILE; 
#}
#
#
#
#
#
## do multiple iterations
#else {
#	
#	my $inputbasename = $argopts{'i'};
#	my $outputbasename = $argopts{'o'};
#	
#	for (my $iter = 0; $iter < $num_iters; $iter++) {
#		open (INPUTFILE, ("$inputbasename iter".$iter.".txt")) || die $!;
#		open (OUTPUTFILE, (">$outputbasename iter".$iter.".txt")) || die $!;
#	
#		# Transform data to [0,1]
#		my %data = ();
#		my $maxscore = 100;
#		#my $maxscore = 10000;
#		my $highest = -9999;
#		foreach my $line (<INPUTFILE>) {
#			chomp $line;
#			(my $id_a, my $id_b, my $score) = split(' ', $line);
#			my $key = ($id_a lt $id_b)?"$id_a|$id_b":"$id_b|$id_a";
#			if ($score > $maxscore) {
#				$score = $maxscore;
#			}
#			if ($doln==1) { 
#				$score = log($score+1); 
#			}
#			if ($score > $highest) {
#				$highest = $score;
#			}
#			$data{$key} = $score;
#		}
#		
#		
#		# normalize
#		foreach my $key (keys %data) {
#			$data{$key} = $data{$key} / $highest;
#		}
#			
#		if ($noweight==1) {		
#			my $count = 0;
#			my @shuffled_array = keys %data;		
#			fisher_yates_shuffle( \@shuffled_array );
#			foreach my $key (@shuffled_array) {
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}
#		}
#	
#		else {		
#			# filter by number to keep
#			my $count = 0;
#			foreach my $key (sort {$data{$b} <=> $data{$a}} keys %data) {
#				(my $id_a, my $id_b) = split(/\|/, $key);
#				print OUTPUTFILE ("$id_a $id_b $data{$key}\n");
#				$count++;
#				if ($num_to_keep ne "all" && $count >= $num_to_keep) {
#					last;
#				}
#			}
#		}
#		
#		close INPUTFILE;
#		close OUTPUTFILE; 
#	}
#	
#}

	

# ----------- fisher_yates_shuffle( \@array ) ---------------
# generate a random permutation of @array in place
sub fisher_yates_shuffle {
        my $array = shift;
        my $i;
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
}
