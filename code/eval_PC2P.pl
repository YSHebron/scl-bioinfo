

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