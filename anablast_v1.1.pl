#!/usr/bin/perl -w
use strict;
$|++; 

sub parse_fasta {
	################################################################################
	# Obtiene la longitud y la secuencia de las secuencias de un multifasta, crea 
	# un fichero que puede ser utilizado por el script wigtoBigWig y devuelve un 
	# hash con la longitud y la secuencia sin saltos de línea
	# ejemplo: &parse_fasta(<path_to_fasta> <workdir>)
	################################################################################

	my $fasta = shift;
	my $workdir = shift;
	
	print STDERR "\nExtracting Chromosones length\n";

	my %seq;
	my $id = "no_id";

	open FASTA, $fasta or die $!;
	while (my $line = <FASTA>) {
		chomp $line;
		$line =~ s/[\n\r]+$//g;
		next if ($line =~ m/^$/);
		if ($line =~ m/^>([^\s\t\n]+)/) {
			$id = $1;
			next;
		}
		else {
			$line = uc($line);
			$seq{$id}{nn} .= $line;
			$seq{$id}{length} += length($line);
		}
	}
	close FASTA or die $!;

	open SIZE, ">$workdir/chrom.size" or die $!;
	print STDERR "Chromosome sizes:\n";
	for my $key (keys %seq) {
		print STDERR "$key\t$seq{$key}{length}\n";
		print SIZE "$key\t$seq{$key}{length}\n";
	}
	close SIZE or die $!;
	print STDERR "chrom.size created\n";
	return \%seq;
}

sub bitscore_filter_frame_split {
	################################################################################
	# filtra los resultados del blast quedandose solo con los que tienen un bitscore
	# de 30 o superior y elimina las repeticiones
	# ejemplo: &bitscore_filter(<path_to_blast> <workdir>)
	################################################################################

	my ($blast, $workdir, $bit_threshold) = @_;
	print STDERR "\nFiltering Blast by Bitscore >= $bit_threshold and spliting by frame...\n";
	
	open BLAST, $blast or die $!;
	my %path = (
		1 => "$workdir/frame1.blast",
		2 =>"$workdir/frame2.blast",
		3 =>"$workdir/frame3.blast",
		-1 => "$workdir/frame-1.blast",
		-2 => "$workdir/frame-2.blast",
		-3 =>"$workdir/frame-3.blast"
	);
	my %count;
	my %file;
	my @frame = (1, 2, 3, -1, -2, -3);

	foreach (@frame) {
		open $file{$_}, ">$path{$_}" or die $!;
	}

	my $prev_line = "";
	while (my $line = <BLAST>) {
		next if ($prev_line eq $line); # Elimina las líneas repetidas consecutivas
		$prev_line = $line;
		my @colums = split (/\t/, $line);
		if ($colums[5] >= $bit_threshold) {
			foreach (@frame) {			
				if ($colums[6] == $_) {
					$count{$_}++;
					print {$file{$_}} $line;
					last;
				}
			}
		}
	}
	close BLAST or die $!;

	foreach (@frame) {
		close $file{$_} or die $!;
	}

	if (!%count) {
		print STDERR "No blast hit with bitscore higher than 36\n"
	}

	foreach (@frame) {
		if ($count{$_}) {
			print STDERR "$count{$_} lines keeped in frame $_\n";
		}
		else {
			print STDERR "0 lines keeped in frame $_\n";
		}
		
	}

	return %path;
}

sub frame_to_wig_and_bigwig {
	################################################################################
	# Analiza los resultados del blast según el frame que le indiques y te crea un
	# fichero .wig con los picos de anablas para ese frame.
	# necesita pasar el blast por la subrutina bitscore_filter_frame_split previamente
	# para que funcione correctamente y que la subrutina fasta_length cree el fichero
	# chrom.size
	# ejemplo: &parse_frame(<frame>, <blast.result>)
	################################################################################

	my $frame = $_[0];
	my $blast = $_[1];
	my $workdir = $_[2];
	my %chr;
	my $wig_path = "$workdir/frame$frame.wig";

	print STDERR "\nparsing frame: $frame\n";

	open BLAST, $blast or die $!;
	if(-z $blast) {
		print STDERR "Empty blast report\n";
		open OUT, ">$wig_path" or die $!;
		print OUT "track type=wiggle_0 name=\"frame$frame\" description=\"frame$frame\" visibility=full graphType=bar\n";
		print STDERR "Empty wig created. Next\n";
	}
	else {
		while(<BLAST>){
			chomp $_;
			my @colums = split (/\t/, $_);
			if ($colums[2] < $colums[3]) {
				for(my $i = $colums[2]; $i <= $colums[3]; $i++){
					$chr{$colums[1]}[$i]++;
				}
			}
			else {
				for(my $i = $colums[3]; $i <= $colums[2]; $i++){
					$chr{$colums[1]}[$i]++;
				}		
			}
		}
		close BLAST;

		open OUT, ">$wig_path" or die $!;
		print OUT "track type=wiggle_0 name=\"frame$frame\" description=\"frame$frame\" visibility=full graphType=bar\n";
		foreach my $key(sort keys %chr){
			print OUT "variableStep chrom=$key\n";
			for my $i(1 .. $#{$chr{$key}}){
				if(!$chr{$key}[$i]){
					print OUT "$i\t0\n";
					next;
				}
				print OUT "$i\t$chr{$key}[$i]\n";
			}
		}
		close OUT;
		print STDERR "frame$frame.wig created\n";
		print STDERR "Creating BigWig...\n";
		`./wigToBigWig $workdir/frame$frame.wig $workdir/chrom.size $workdir/frame$frame.bw`;
		print STDERR "BigWig Created\n";
	}
	return $wig_path;
}

sub reverse_comp {
	################################################################################
	# Obtiene la reversa complementaria de una secuencia dada
	# ejemplo: &reverse_comp(<secuenca de DNA>)
	################################################################################
	my $seq = shift;
	$seq =~ tr/acgtACGT/tgcaTGCA/;
	$seq = reverse $seq;
	return $seq;
}

sub dna2aa {
	################################################################################
	# Traduce una secuencia de ADN a proteínas, se le debe indicar en que fase se 
	# encuentra (1, 2, 3)
	# ejemplo: &dna2aa(<DNA_seq>, <frame>)
	################################################################################
	#código genético
	my %genetic_code = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '*',    # Stop
		'TAG' => '*',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '*',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
	);
	my $prot;
	my $dna = shift;
	my $frame = shift;
	my $length = length($dna);
	my $res = $length%3;
	my $c = 0;
	for(my $x=$frame-1; $x < $length - $res; $x += 3){
		my $codon = substr($dna,$x,3);
		$codon = uc($codon);
		if (!$genetic_code{$codon}) {
			$prot .= "?";
			$c++;
			#print STDERR "Warning: Extrange nucleotide ($codon) detected in dna2aa subrutine. AA ? created\n";
			next;
		}
		$prot .= $genetic_code{$codon};	
	}
	if ($c != 0) {
		print STDERR "Warning: $c \"?\" Amino acids created\n";
	}
	return $prot;	
}

sub frame_calc {
	################################################################################
	# Calcula el frame de una subsecuencia a partir del frame original de la 
	# secuencia original teniendo en cuenta la reversa complementaria
	# ejemplo: &frame_calc(<frame>, <start>, <end>, <sequence_length>)
	################################################################################

	my $ori_frame = shift;
	my $start = shift;
	my $end = shift;
	my $seq_length = shift;
	
	
	my $frame;
	my $res;

	if ($ori_frame > 0) {
		$res = ($start -1) % 3;
	
		if($res == 0){
			$frame = $ori_frame;
		}
		elsif($res == 1){
			if($ori_frame == 1){
				$frame = 3;	
			}
			elsif($ori_frame == 2){
				$frame = 1;	
			}
			elsif($ori_frame == 3){
				$frame = 2;	
			}	
		}
		elsif($res == 2){
			if($ori_frame == 1){
				$frame = 2;	
			}
			elsif($ori_frame == 2){
				$frame = 3;	
			}
			elsif($ori_frame == 3){
				$frame = 1;	
			}	
		}
	} 
	else {

		$res = ($seq_length - $end) % 3;
	
		if($res == 0){
			if($ori_frame == -1){
				$frame = 1;	
			}
			elsif($ori_frame == -2){
				$frame = 2;	
			}
			elsif($ori_frame == -3){
				$frame = 3;	
			}
		}
		elsif($res == 1){
			if($ori_frame == -1){
				$frame = 3;	
			}
			elsif($ori_frame == -2){
				$frame = 1;	
			}
			elsif($ori_frame == -3){
				$frame = 2;	
			}	
		}
		elsif($res == 2){
			if($ori_frame == -1){
				$frame = 2;	
			}
			elsif($ori_frame == -2){
				$frame = 3;	
			}
			elsif($ori_frame == -3){
				$frame = 1;	
			}	
		}
	}
	# print "\n$ori_frame\t($res) => $frame\n";
	return $frame;
}

sub extract_peaks {
	################################################################################
	# Extrae los datos de los picos de un marco de lectura determinado con unos 
	# determinados valores desde un fichero WIG
	# Se puede modificar la altura mínima del pico y desde que altura se empieza o se
	# deja de condiderar que es un pico
	# ejemplo: &extract_pic(<path_to_wig>, <sequence_reference>)
	################################################################################
	my $wig = shift;
	my $seq_ref = shift;
	my %seq = %{$seq_ref};

	print STDERR "\nExtracting peaks from WIG: $wig\n";
	my $chrom;
	my $frame;
	my $start;
	my $end;
	my $f = 0;
	my $f2 = 0;
	my %pic;
	my $top = 0;
	my $c = 0;
	my $position;
	
	open WIG, "$wig" or die $!;
	
	while(<WIG>){
		chomp;
		if($_ =~ m/variableStep chrom=(.+)/){
			if($f == 1 and $f2 == 1) {
				$f = 0;
				$f2 = 0;
				$c++;
				$end = $position;
				$pic{"frame_$frame"."_peak_$c"}{frame} = $frame;
				$pic{"frame_$frame"."_peak_$c"}{chrom} = $chrom;
				$pic{"frame_$frame"."_peak_$c"}{start} = $start;
				$pic{"frame_$frame"."_peak_$c"}{end} = $end;
				$pic{"frame_$frame"."_peak_$c"}{position} = "$chrom:$start..$end";
				$pic{"frame_$frame"."_peak_$c"}{top} = $top;
				$pic{"frame_$frame"."_peak_$c"}{length} = $end - $start + 1;
				if ($frame >= 1) {
					$pic{"frame_$frame"."_peak_$c"}{strand} = "+";
				}
				else {
					$pic{"frame_$frame"."_peak_$c"}{strand} = "-";
				}
				
				$top = 0;
			}
			$chrom = $1;
			next;
		}
		elsif ($_ =~ m/track type=wiggle_0 name=\"frame(.+)\" description=\"frame/) {
			$frame = $1;
			next;
		}
		next unless($_ =~ m/^(\d+)\t(\d+)$/);
		$position = $1;
		if($2 >= 15 and $f == 0){
			$start = $1;
			$f = 1;
		}
		if($2 >= 20){ #Altura mínima para considerarlo un pico
			$f2 = 1;
			if($2 > $top){ 
				$top = $2;
			} 	
		}
		elsif($2 <= 15 and $f == 1){
			$end = $1;
			$f = 0;
			# print "f$frame\ts$start\tt$top\te$end\n"; ## prueba para ver donde comienzan y terminan los picos
			if($f2 == 1){
				$f2 = 0;
				$c++;
				$pic{"frame_$frame"."_peak_$c"}{frame} = $frame;
				$pic{"frame_$frame"."_peak_$c"}{chrom} = $chrom;
				$pic{"frame_$frame"."_peak_$c"}{start} = $start;
				$pic{"frame_$frame"."_peak_$c"}{end} = $end;
				$pic{"frame_$frame"."_peak_$c"}{position} = "$chrom:$start..$end";
				$pic{"frame_$frame"."_peak_$c"}{top} = $top;
				$pic{"frame_$frame"."_peak_$c"}{length} = $end - $start + 1;
				if ($frame >= 1) {
					$pic{"frame_$frame"."_peak_$c"}{strand} = "+";
				}
				else {
					$pic{"frame_$frame"."_peak_$c"}{strand} = "-";
				}
				
				$top = 0;
			}	
		}
	}

	if ($f == 1 and $f2 == 1) {
		$f = 0;
		$f2 = 0;
		$c++;
		$end = $position;
		$pic{"frame_$frame"."_peak_$c"}{frame} = $frame;
		$pic{"frame_$frame"."_peak_$c"}{chrom} = $chrom;
		$pic{"frame_$frame"."_peak_$c"}{start} = $start;
		$pic{"frame_$frame"."_peak_$c"}{end} = $end;
		$pic{"frame_$frame"."_peak_$c"}{position} = "$chrom:$start..$end";
		$pic{"frame_$frame"."_peak_$c"}{top} = $top;
		$pic{"frame_$frame"."_peak_$c"}{length} = $end - $start + 1;
		if ($frame >= 1) {
			$pic{"frame_$frame"."_peak_$c"}{strand} = "+";
		}
		else {
			$pic{"frame_$frame"."_peak_$c"}{strand} = "-";
		}
	}

	close WIG or die $!;
	print STDERR "Total peaks in $wig = $c\n";

	print STDERR "Extracting $wig peak sequences\n";
	foreach my $key (keys %pic) {
		$pic{$key}{nn} = substr($seq{$pic{$key}{chrom}}{nn}, $pic{$key}{start} -1, $pic{$key}{end} - $pic{$key}{start} +1);
		$pic{$key}{nn} = uc($pic{$key}{nn});
		# print STDERR "$key:$pic{$key}{start}..$pic{$key}{end}\n$pic{$key}{nn}\n"; ## prueba para ver como se ve la secuencia
		if ($pic{$key}{frame} < 0) {
			$pic{$key}{nn} = &reverse_comp($pic{$key}{nn});
		}
		
		$pic{$key}{aa} = &dna2aa($pic{$key}{nn}, &frame_calc($pic{$key}{frame}, $pic{$key}{start}, $pic{$key}{end}, $seq{$pic{$key}{chrom}}{length}));
		my @aa = split ("", $pic{$key}{aa});
		my $stop_codon = 0;
		foreach my $a (@aa) {
			if ($a eq "*") {
				$stop_codon++;
			}
		}
		$pic{$key}{stop_codon} = $stop_codon;
	}
	return \%pic;
}

sub peaks_to_tsv {
	################################################################################
	# Guarda los picos en un fichero de texto en formato tsv
	# ejemplo: &peaks_to_tsv(<ref_peaks>, <path_to_tsv>)
	################################################################################

	print STDERR "\nPrinting peaks in tsv\n";
	my $ref_peaks = shift;
	my %peaks = %{$ref_peaks};
	my $tsv = shift;

	open OUT, ">$tsv" or die $!;
	foreach my $key (keys %peaks) {
		print OUT "id\t";
		foreach my $key2 (sort keys %{$peaks{$key}}) {
			print OUT "$key2\t";
		}
		print OUT "\n";
		last;
	}
	foreach my $key (keys %peaks) {
		print OUT "$key\t";
		foreach my $key2 (sort keys %{$peaks{$key}}) {
			print OUT "$peaks{$key}{$key2}\t";
		}
		print OUT "\n";
	}
	close OUT or die $!;
}

sub tsv_to_peaks {

	my $tsv = shift;
	my %peaks;
	my @header;

	print STDERR "Importing peaks from $tsv\n";

	open TSV, "$tsv" or die $!;
	while (<TSV>) {
		if ($_ =~ m/^id/) {
			@header = split(/\t/, $_);
			next;
		}
		my @data = split(/\t/, $_);

		$peaks{$data[0]}{$header[1]} = $data[1];
		$peaks{$data[0]}{$header[2]} = $data[2];
		$peaks{$data[0]}{$header[3]} = $data[3];
		$peaks{$data[0]}{$header[4]} = $data[4];
		$peaks{$data[0]}{$header[5]} = $data[5];
		$peaks{$data[0]}{$header[6]} = $data[6];
		$peaks{$data[0]}{$header[7]} = $data[7];
		$peaks{$data[0]}{$header[8]} = $data[8];
		$peaks{$data[0]}{$header[9]} = $data[9];
		$peaks{$data[0]}{$header[10]} = $data[10];
		$peaks{$data[0]}{$header[11]} = $data[11];
	}
	return \%peaks;
}

sub peaks_to_gff3 {
	################################################################################
	# Guarda los picos en un fichero en formato gff3
	# ejemplo: &peaks_to_gff3(<ref_peaks>, <path_to_gff>)
	################################################################################

	print STDERR "\nPrinting peaks in GFF3\n";
	my $ref_peaks = shift;
	my %peaks = %{$ref_peaks};
	my $gff = shift;

	open OUT, ">$gff" or die $!;
	foreach my $key (keys %peaks) {
		print OUT "$peaks{$key}{chrom}\tAnABlast\tPeak\t$peaks{$key}{start}\t$peaks{$key}{end}\t$peaks{$key}{top}\t$peaks{$key}{strand}\t.\tID=$key;frame=$peaks{$key}{frame};stop_codons=$peaks{$key}{stop_codon};amino_acids=<textarea readonly=\"\" rows=\"7\" cols=\"62\" class=\"fasta\">$peaks{$key}{aa}</textarea></div>,\n";
	}
	close OUT or die $!;
}

sub peaks_to_fasta {
	################################################################################
	# Guarda la secuencia de aa de los picos en un fichero en formato fasta
	# ejemplo: &peaks_to_fasta(<ref_peaks>, <path_to_fasta>)
	################################################################################
	
	print STDERR "\nPrinting peak seqs in FASTA\n";
		my $ref_peaks = shift;
      my %peaks = %{$ref_peaks};
      my $fasta = shift;

   open OUT, ">$fasta" or die $!;
	foreach my $key (keys %peaks) {
		print OUT ">$key $peaks{$key}{chrom}:$peaks{$key}{start}..$peaks{$key}{end}\n";
		print OUT "$peaks{$key}{aa}\n";
	}
	close OUT or die $!;
}

sub peaks_elong_orf_search {
	################################################################################
	# Busca ORF en los picos
	# ejemplo: &orf_search(<ref_peaks>, <sequence_reference>)
	################################################################################

	my $ref_peaks = shift;
	my %peaks = %{$ref_peaks};
	my $seq_ref = shift;
	my %seq = %{$seq_ref};

	print STDERR "\nSearching ORF in peaks\n";

	my %orf;
	my $c;

	foreach my $key (keys %peaks) {
		
		print "***$key\n";
		my $sub_frame = &frame_calc($peaks{$key}{frame}, $peaks{$key}{start}, $peaks{$key}{end}, $seq{$peaks{$key}{chrom}}{length});
		
		my $nn;
		my $start;
		my $end;
		if ($peaks{$key}{frame} > 0) {
			$start = $peaks{$key}{start} + $sub_frame -1;
			my $rest = ($peaks{$key}{end} - $start + 1) % 3;
			$end = $peaks{$key}{end} - $rest;
			
			while (1) {
				if ($start <= 1) {
					$start = $peaks{$key}{frame};
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $start -1, 3);
				if ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG" or $codon eq "NNN") {
					last;
				}
				else {
					$start -= 3;
				}
			}
			while (1) {
				if ($end >= $seq{$peaks{$key}{chrom}}{length} - 3) {
					$end = $seq{$peaks{$key}{chrom}}{length};
					last
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $end, 3);
				if ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG" or $codon eq "NNN") {
					$end += 3;
					last;
				}
				else {
					$end += 3;
				}
			}
			$nn = substr ($seq{$peaks{$key}{chrom}}{nn}, $start -1, $end - $start +1);
			print "$nn\n";
			my $f = 0;
			for (my $i = 0; $i <= length($nn); $i += 3) {
				my $codon = substr($nn, $i, 3);
				print "$codon | ";
				if ($codon eq "ATG" or $codon eq "ATA" or $codon eq "ATC" or $codon eq "GTG" or $codon eq "TTG" or $codon eq "CTG" or $codon eq "ATT") {
					if ($f == 0) {
						print "<= start | ";
						$c++;
						$f = 1;
						$orf{"$key"."_ORF_$c"}{start} = $start + $i;
						$orf{"$key"."_ORF_$c"}{parent} = $key;
						$orf{"$key"."_ORF_$c"}{chrom} = $peaks{$key}{chrom};
						$orf{"$key"."_ORF_$c"}{strand} = $peaks{$key}{strand};
						$orf{"$key"."_ORF_$c"}{frame} = $peaks{$key}{frame};
					}
				}
				elsif ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG") {
					if ($f == 1) {
						print "<= end\n";
						$f = 0;
						$orf{"$key"."_ORF_$c"}{end} = $start + $i +2;
						$orf{"$key"."_ORF_$c"}{nn} = substr ($seq{$peaks{$key}{chrom}}{nn}, $orf{"$key"."_ORF_$c"}{start} -1, $orf{"$key"."_ORF_$c"}{end} - $orf{"$key"."_ORF_$c"}{start} + 1);
						$orf{"$key"."_ORF_$c"}{aa} = &dna2aa($orf{"$key"."_ORF_$c"}{nn}, 1);
						delete $orf{"$key"."_ORF_$c"} if (length $orf{"$key"."_ORF_$c"}{aa} <= 5);
					}
				}
			}
		}
		else {
			$end = $peaks{$key}{end} - $sub_frame + 1;
			my $rest = ($end - $peaks{$key}{start} + 1) % 3;
			$start = $peaks{$key}{start} + $rest;

			while (1) {
				if ($start <= 1) {
					$start = 1;
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $start -1, 3);
				my $rcodon = &reverse_comp($codon);
				if ($rcodon eq "TGA" or $rcodon eq "TAA" or $rcodon eq "TAG" or $rcodon eq "NNN") {
					last;
				}
				else {
					$start -= 3;
				}
			}
			while (1) {
				if ($end >= $seq{$peaks{$key}{chrom}}{length} - 3) {
					$end = $seq{$peaks{$key}{chrom}}{length} + $peaks{$key}{frame} + 1;
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $end, 3);
				
				my $rcodon = &reverse_comp($codon);
				if ($rcodon eq "TGA" or $rcodon eq "TAA" or $rcodon eq "TAG" or $rcodon eq "NNN") {
					$end += 3;
					last;
				}
				else {
					$end += 3;
				}
			}
			$nn = substr ($seq{$peaks{$key}{chrom}}{nn}, $start -1, $end - $start +1);
			#print "\n$peaks{$key}{chrom}:$start..$end\t$peaks{$key}{frame}\n$nn\n";
			$nn = &reverse_comp($nn);
			
			my $f = 0;

			for (my $i = 0; $i <= length($nn); $i += 3) {
				my $codon = substr($nn, $i, 3);
				if ($codon eq "ATG" or $codon eq "ATA" or $codon eq "ATC" or $codon eq "GTG" or $codon eq "TTG" or $codon eq "CTG" or $codon eq "ATT") {
					if ($f == 0) {
						$c++;
						$f = 1;
						$orf{"$key"."_ORF_$c"}{end} = $end - $i;
						$orf{"$key"."_ORF_$c"}{parent} = $key;
						$orf{"$key"."_ORF_$c"}{chrom} = $peaks{$key}{chrom};
						$orf{"$key"."_ORF_$c"}{strand} = $peaks{$key}{strand};
						$orf{"$key"."_ORF_$c"}{frame} = $peaks{$key}{frame};
					}
				}
				elsif ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG") {
					if ($f == 1) {
						$f = 0;
						$orf{"$key"."_ORF_$c"}{start} = $end - $i - 2;
						$orf{"$key"."_ORF_$c"}{nn} = &reverse_comp(substr ($seq{$peaks{$key}{chrom}}{nn}, $orf{"$key"."_ORF_$c"}{start} -1, $orf{"$key"."_ORF_$c"}{end} - $orf{"$key"."_ORF_$c"}{start} +1));
						$orf{"$key"."_ORF_$c"}{aa} = &dna2aa($orf{"$key"."_ORF_$c"}{nn}, 1);
						delete $orf{"$key"."_ORF_$c"} if (length $orf{"$key"."_ORF_$c"}{aa} <= 5);
					}
				}
			}
		}
	}
	foreach my $key (keys %orf) {
		if (!$orf{$key}{start} or !$orf{$key}{end}) {
			delete $orf{$key};
		}
	}
	return \%orf;
}

sub orf_to_gff3 {
	################################################################################
	# Guarda los ORF en formato GFF3
	# ejemplo: &orf_to_gff3(<ref_orf>, <path_to_gff>)
	################################################################################

	print STDERR "\nPrinting ORF in GFF3\n";
	my $ref_orf = shift;
	my %orf = %{$ref_orf};
	my $gff = shift;
	open OUT, ">$gff" or die $!;
	foreach my $key (keys %orf) {
		print OUT "$orf{$key}{chrom}\tAnABlast\tORF\t$orf{$key}{start}\t$orf{$key}{end}\t.\t$orf{$key}{strand}\t.\tID=$key;frame=$orf{$key}{frame};amino_acids=<textarea readonly=\"\" rows=\"7\" cols=\"62\" class=\"fasta\">$orf{$key}{aa}</textarea></div>,\n"
	}
	close OUT or die $!;
}

sub sma3s {
	################################################################################
	# Ejecuta Sma3s e integra las anotaciones en el gff de picos como subfeatures
	# ejemplo: &sma3s(<path_to_work_dir>, <path_to_uniref90_db>)
	################################################################################
	print STDERR "\nRunning Sm3s\n";

	my $wd = shift;
	my $db = shift;
	
	`/home/apps/Sma3s_v2/sma3s_v2.pl -i $wd/peaks.faa -d $db -num_threads 8 -go -goslim -source -uniprot &> $wd/sma3s.log`;
	
	my %annot;
	open TSV, "$wd/peaks_UniSpr_go_goslim_source.tsv" or die $!;
	while (<TSV>) {
		next if ($_ =~ m/^#/);
		chomp;
		
		my ($name, $gene, $descript, $ec, $p_id, $p_name, $f_id, $f_name, $c_id, $c_name, $kw, $pw, $goslim) = split(/\t/, $_);
		my ($id, $chrom, $start, $end) = split(/ |:|\.\./, $name);
		
		$annot{$chrom}{$id} = "$chrom\tSma3s\tannot\t$start\t$end\t.\t.\t.\tParent=$id;gene=$gene ($descript);ENZYME=$ec;GO (biological process)=";
		
		if ($p_id =~ m/.+/) {
			my @a1 = split(/;/, $p_id);
			my @a2 = split(/;/, $p_name);
			for (my $i = 0; $i <= $#a1; $i++) {
				$annot{$chrom}{$id} .= "$a1[$i] ($a2[$i])<br>";
			}
		}
		
		$annot{$chrom}{$id} .= ";GO (molecular function)=";
		
		if ($f_id =~ m/.+/) {
			my @a1 = split(/;/, $f_id);
			my @a2 = split(/;/, $f_name);
			for (my $i = 0; $i <= $#a1; $i++) {
				$annot{$chrom}{$id} .= "$a1[$i] ($a2[$i])<br>";
			}
		}
		
		$annot{$chrom}{$id} .= ";GO (cellular component)=";
		
		if ($c_id =~ m/.+/) {
			my @a1 = split(/;/, $c_id);
			my @a2 = split(/;/, $c_name);
			for (my $i = 0; $i <= $#a1; $i++) {
				$annot{$chrom}{$id} .= "$a1[$i] ($a2[$i])<br>";
			}
		}
		
		$annot{$chrom}{$id} .= ";Keywords=";
		
		foreach my $key (split(/;/, $kw)) {
			$annot{$chrom}{$id} .= "$key<br>";
		}
		
		$annot{$chrom}{$id} .= ";Pathways=";
		
		foreach my $key (split(/;/, $pw)) {
			$annot{$chrom}{$id} .= "$key<br>";
		}
		
		$annot{$chrom}{$id} .= ";GO slim=";
		
		foreach my $key (split(/;/, $goslim)) {
			$annot{$chrom}{$id} .= "$key<br>";
		}
		
		$annot{$chrom}{$id} .= ";";
	}
	close TSV or die $!;
	
	my %top_blast;
	my $sid = " ";
	open BLAST, "$wd/peaks_UniSpr.blast" or die $!;
	while (<BLAST>) {
		chomp;
		my ($qseqid, $sseqid, $pident, $qcovs, $length, $evalue) = split(/\t/, $_);
		next if ($qseqid eq $sid);
		$top_blast{$qseqid}{sseqid} = $sseqid;
		$top_blast{$qseqid}{pident} = $pident;
		$top_blast{$qseqid}{qcovs} = $qcovs;
		$top_blast{$qseqid}{"length"} = $length;
		$top_blast{$qseqid}{evalue} = $evalue;
		$sid = $qseqid;
	}
	close BLAST or die $!;
	
	my %faa;
	my $c = 0;
	open FAA, "$wd/peaks.faa" or die $!;
	while (<FAA>) {
		next unless ($_ =~ m/^>/);
		chomp;
		$c++;
		my ($odd, $id, $chrom, $start, $end) = split(/>| |:|\.\./, $_);
		$faa{$chrom}{$id}{start} = $start;
		$faa{$chrom}{$id}{end} = $end;
		$faa{$chrom}{$id}{sseqid} = $top_blast{"s".$c}{sseqid};
		$faa{$chrom}{$id}{pident} = $top_blast{"s".$c}{pident};
		$faa{$chrom}{$id}{qcovs} = $top_blast{"s".$c}{qcovs};
		$faa{$chrom}{$id}{"length"} = $top_blast{"s".$c}{"length"};
		$faa{$chrom}{$id}{evalue} = $top_blast{"s".$c}{evalue};
	}
	close FAA or die $!;
	
	open GFF, "$wd/peaks.gff" or die $!;
	open GFF2, ">$wd/temp.gff" or die $!;
	while (<GFF>) {
		print GFF2 $_;
		my ($chrom, $o, $t, $start, $end, $top, $strand, $p, $odd, $id) = split(/\t|=|;/, $_);
		if ($annot{$chrom}{$id}) {
			print GFF2 "$annot{$chrom}{$id}\n";
		}
		if ($faa{$chrom}{$id}{sseqid}) {
			print GFF2 "$chrom\tBlast\tTop BLAST\t$start\t$end\t.\t.\t.\tParent=$id;UniProt ID=$faa{$chrom}{$id}{sseqid};identity (%)=$faa{$chrom}{$id}{pident};coverage (%)=$faa{$chrom}{$id}{qcovs};Alignment length=$faa{$chrom}{$id}{length};E-value=$faa{$chrom}{$id}{evalue};\n";
		}
	}
	close GFF or die $!;
	close GFF2 or die $!;
	`mv $wd/temp.gff $wd/peaks.gff`;
}

sub gff_index {
	################################################################################
	# Parte de un fichero en formato gff3 y lo indexa en memoria para poder trabajar
	# con el devolviendo un hash de hashes
	# ejemplo: &gff_index(<path to gff3>)
	################################################################################
	my ($gff) = @_;
	my %hash;

	print STDERR "\nIndexing gff $gff\n";

	open GFF, $gff or die $!;
	my $c;
	while (<GFF>) {
		chomp;
		$c++;
		next if ($_ =~ m/^#/ or $_ =~ m/^$/);
		my ($chrom, $source, $type, $start, $end, $point, $strand, $point2, $qualifier) = split (/\t/,$_);
		
		my $start_range = 10000 * (int($start/10000) + 1);
		my $end_range = 10000 * (int($end/10000) + 1);

		for (my $i = $start_range; $i <= $end_range; $i += 10000) {
			$hash{$chrom}{$strand}{$i}{$c}{chrom} = $chrom;
			$hash{$chrom}{$strand}{$i}{$c}{source} = $source;
			$hash{$chrom}{$strand}{$i}{$c}{type} = $type;
			$hash{$chrom}{$strand}{$i}{$c}{start} = $start;
			$hash{$chrom}{$strand}{$i}{$c}{end} = $end;
			$hash{$chrom}{$strand}{$i}{$c}{strand} = $strand;
			$hash{$chrom}{$strand}{$i}{$c}{qualifier} = $qualifier;
			$hash{$chrom}{$strand}{$i}{$c}{length} = $end - $start +1;
			$hash{sources}{"$type ($source)"} = 0;
		}
	}
	$hash{total} = $c;
	close GFF or die $!;
	print STDERR "Total item $c in $gff\n";
	return \%hash;
}

sub peaks_vs_gff {
	################################################################################
	# Parte de el gff indexado en memoria y del hash de picos generado por la 
	# subrutina extract_peaks, imprime un fichero tsv con la superposición de los
	# picos con los diferentes item del gff y devuelve la referencia del tsv
	# ejemplo: &peak_vs_gff(<path to gff3> <%peaks> <work dir>)
	################################################################################
	my ($ref_gff, $ref_peaks, $path) = @_;
	my %gff = %{$ref_gff};
	my %peaks = %{$ref_peaks};
	my %tsv;

	print STDERR "\nCreating TSV peaks vs gff\n";

	foreach my $key (keys %peaks) {
		my $start_range = 10000 * (int($peaks{$key}{start}/10000) + 1);
		my $end_range = 10000 * (int($peaks{$key}{end}/10000) + 1);

		my $p_chrom = $peaks{$key}{chrom};
		my $p_strand = $peaks{$key}{strand};
		my $p_start = $peaks{$key}{start};
		my $p_end = $peaks{$key}{end};

		for (my $i = $start_range; $i <= $end_range; $i += 10000) {
			foreach my $key2 (keys %{$gff{$p_chrom}{$p_strand}{$i}}) {

				my $g_start = $gff{$p_chrom}{$p_strand}{$i}{$key2}{start};
				my $g_end = $gff{$p_chrom}{$p_strand}{$i}{$key2}{end};
				my $g_type = $gff{$p_chrom}{$p_strand}{$i}{$key2}{type};
				my $g_source = $gff{$p_chrom}{$p_strand}{$i}{$key2}{source};
				my $g_qualifier = $gff{$p_chrom}{$p_strand}{$i}{$key2}{qualifier};

				next if ($p_start > $g_end or $p_end < $g_start);

				$tsv{$key}{qualifier} .= "$g_qualifier\|\|";
				$tsv{$key}{"$g_type ($g_source)"} .= "$g_start..$g_end;";
			}
		}
	}

	# Calcula los porcentajes de cobertura de los picos con los item del gff
	foreach my $key3 (keys %tsv) {
		foreach my $key4 (keys %{$tsv{$key3}}) {

			next if ($key4 eq "qualifier");
			
			my %nn;
			my $c;
			my $n;
			
			for my $p ($peaks{$key3}{start}..$peaks{$key3}{end}){
				$nn{$p} = 0;
			}

			my @coord = split(";",$tsv{$key3}{$key4});
			
			foreach my $c (@coord){
				my ($g_start, $g_end) = split(/\.\./, $c);
				for my $p ($peaks{$key3}{start}..$peaks{$key3}{end}){
					if($p >= $g_start and $p <= $g_end){
						$nn{$p} = 1; 
					}
				}
			}

			foreach my $key5 (sort keys %nn){
				$n++;
				if($nn{$key5} == 1){
					$c++;
				}	
			}
			$tsv{$key3}{$key4} = $c/$n;
		}
	}

	# Rellena los picos que no tienen coincidencia con ceros
	foreach my $k (keys %peaks) {
		if (!$tsv{$k}{qualifier}) {
			$tsv{$k}{qualifier} = "no_qualifier";
		}		
		foreach my $s (keys %{$gff{sources}}){
			if (!$tsv{$k}{$s}) {
				$tsv{$k}{$s} = 0;
			}
		}
	}

	# Imprime el fichero TSV
	my $file = $path."/peaks_vs_gff.tsv";
	open FILE, ">$file" or die $!;

	foreach my $key6 (keys %peaks) {
		foreach my $key7 (sort keys %{$peaks{$key6}}) {
			print FILE "$key7\t";
		}
		last;
	}

	foreach my $key8 (keys %tsv) {
		foreach my $k (sort keys %{$tsv{$key8}}) {
			print FILE "$k\t";
		}
		last;
	}
	print FILE "\n";

	foreach my $key9 (sort keys %peaks) {
		foreach my $key10 (sort keys %{$peaks{$key9}}) {
			print FILE "$peaks{$key9}{$key10}\t";
		}
		foreach my $key11 (sort keys %{$tsv{$key9}}) {
			if ($key11 eq "qualifier") {
				print FILE "$tsv{$key9}{$key11}\t";
				next;
			}
			printf FILE "%.2f\t", "$tsv{$key9}{$key11}\t";
		}
		print FILE "\n";
	}
	close FILE or die $!;
}

sub sensitivity_specificity {
	################################################################################
	# ejemplo: &sensitivity_specificity(<ref_gff3> <ref_peaks> <min high>)
	################################################################################

	my ($ref_gff, $ref_peaks, $min_top, $seq_ref) = @_;
	print STDERR "\nSensitivity and Specificity (min. top: $min_top)\n";
	
	my %chrom = %{$seq_ref};
	my %peaks = %{$ref_peaks};
	my %gff = %{$ref_gff};

	my %cds;
	my %cds_p;
	my %peak_cov;
	my %peaks_p;
	my %result;
	
	# Extraer los CDS a partir del gff indexado
	foreach my $k (keys %gff) {
		next if ($k eq "sources" or $k eq "total");
		#print "$gff{$k}\n";
		foreach my $k2 (keys %{$gff{$k}}) {
			foreach my $k3 (keys %{$gff{$k}{$k2}}) {
				foreach my $k4 (keys %{$gff{$k}{$k2}{$k3}}) {
					$cds{"$k:$gff{$k}{$k2}{$k3}{$k4}{start}..$gff{$k}{$k2}{$k3}{$k4}{end}"} = 2;
					for my $n ($gff{$k}{$k2}{$k3}{$k4}{start}..$gff{$k}{$k2}{$k3}{$k4}{end}) {
						$cds_p{$k}{$n} = 1;
					}
				}
			}
		}
	}

	$result{total_cds} = keys %cds;
	print STDERR "Total uniq CDS = $result{total_cds}\n";

	# Recorre los picos
	foreach my $key (keys %peaks) {
		#print STDERR "$peaks{$key}{position}\t";

		next if ($peaks{$key}{top} < $min_top);
		$result{total}++;
		
		#Calcula el rango para cada pico
		my $start_range = 10000 * (int($peaks{$key}{start}/10000) + 1);
		my $end_range = 10000 * (int($peaks{$key}{end}/10000) + 1);
		
		my $p_chrom = $peaks{$key}{chrom};
		my $p_strand = $peaks{$key}{strand};
		my $p_start = $peaks{$key}{start};
		my $p_end = $peaks{$key}{end};
		my $p_length = $peaks{$key}{length};
		

		for (my $i = $start_range; $i <= $end_range; $i += 10000) {
			# Recorre los cds de un rango determinado de la cadena fordware del pico
			foreach my $key2 (keys %{$gff{$p_chrom}{$p_strand}{$i}}) {

				my $g_start = $gff{$p_chrom}{$p_strand}{$i}{$key2}{start};
				my $g_end = $gff{$p_chrom}{$p_strand}{$i}{$key2}{end};
				my $g_type = $gff{$p_chrom}{$p_strand}{$i}{$key2}{type};
				my $g_source = $gff{$p_chrom}{$p_strand}{$i}{$key2}{source};
				my $g_qualifier = $gff{$p_chrom}{$p_strand}{$i}{$key2}{qualifier};
				my $g_length = $gff{$p_chrom}{$p_strand}{$i}{$key2}{length};

				next if ($p_start > $g_end or $p_end < $g_start);
				#print STDERR "\t\t$p_chrom:$g_start..$g_end\n";
				my $cover;
				for my $n ($g_start..$g_end) {
					if ($n >= $p_start and $n <= $p_end) {
						$cover++;
					}
				}
				if ($cover/$p_length >= 0.2) {
					$peak_cov{$key} = 1;
					$cds{"$p_chrom:$g_start..$g_end"} = 1;
				}
				elsif ($cover/$g_length >= 0.2) {
					$peak_cov{$key} = 1;
					$cds{"$p_chrom:$g_start..$g_end"} = 1;
				}
			}
			my $rev_strand;
			if ($p_strand eq "+") {
				$rev_strand = "-";
			}
			else {
				$rev_strand = "+";
			}

			# Recorre los cds de un rando determinado en la cadena reversa del pico
			foreach my $key2 (keys %{$gff{$p_chrom}{$rev_strand}{$i}}) {
				my $g_start = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{start};
				my $g_end = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{end};
				my $g_type = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{type};
				my $g_source = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{source};
				my $g_qualifier = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{qualifier};
				my $g_length = $gff{$p_chrom}{$rev_strand}{$i}{$key2}{length};
				
				next if ($p_start > $g_end or $p_end < $g_start);
				#print STDERR "\t\t$p_chrom:$g_start..$g_end\n";

				my $cover;
				for my $n ($g_start..$g_end) {
					if ($n >= $p_start and $n <= $p_end) {
						$cover++;
					}
				}
				if ($cover/$p_length >= 0.2) {
					$peak_cov{$key} = 1;
					$cds{"$p_chrom:$g_start..$g_end"} = 1;
				}
				elsif ($cover/$g_length >= 0.2) {
					$peak_cov{$key} = 1;
					$cds{"$p_chrom:$g_start..$g_end"} = 1;
				}
			}
			if (!$peak_cov{$key}) {
				foreach my $n ($p_start..$p_end) {
					$peaks_p{$p_chrom}{$n} = 1;
				}
			}
		}
	}

	print STDERR "Total peaks >= $min_top = $result{total}\n";

	$result{fp} = 0;
	#open FILE, ">peaks_without_cds_$min_top.txt" or die $!;
	foreach my $key2 (keys %peaks) {
		next if ($peaks{$key2}{top} < $min_top);
		if ($peak_cov{$key2}) {
			#$result{tp_peaks}++;
		}
		else {
			$result{fp}++;
	#		print FILE "$peaks{$key2}{chrom}:$peaks{$key2}{start}..$peaks{$key2}{end}\n";
		}
	}
	#close FILE or die $!;
	
	$result{tp} = 0;
	$result{fn} = 0;
	#open FILE, ">cds_without_peak_$min_top.txt" or die $!;
	foreach my $key3 (keys %cds) {
		if ($cds{$key3} == 1) {
			$result{tp}++;
		}
		else {
	#		print FILE "$key3\n";
			$result{fn}++;
		}
	}
	#close FILE or die $!;

	$result{sensitivity} = $result{tp}/($result{tp} + $result{fn});
	$result{specificity} = $result{tp}/($result{tp} + $result{fp});

	%cds = ();
	%peak_cov = ();
	%peaks = ();
	%gff = ();

	###########################################
	# my %no_cds;
	# foreach my $key (sort keys %chrom) {
	# 	print "Parsing chrom: $key\n";
	# 	for my $n (1..$chrom{$key}{length}) {
	# 		if ($cds_p{$key}{$n}) {
	# 			$no_cds{$key}{$n} = 0;
	# 		}
	# 		else {
	# 			$no_cds{$key}{$n} = 1;
	# 			if ($peaks_p{$key}{$n}) {
	# 				$no_cds{$key}{$n} = 2;
	# 			}
	# 		}
	# 		#print "$key\t$n\t$no_cds{$key}{$n}\n";
	# 	}
	# }

	# my $c;
	# my $f = 0;
	# my %true_n;
	# #open FILE, ">true_netative_$min_top.txt" or die $!;
	# foreach my $key (sort keys %no_cds) {
	# 	for my $n (1..$chrom{$key}{length}) {
	# 		if ($no_cds{$key}{$n} == 1 and $f == 0) {
	# 			$f = 1;
	# 			$c++;
	# 			$true_n{$c} = "$key:$n.."
	# 		}
	# 		if (($no_cds{$key}{$n} != 1 and $f == 1) or ($n == $chrom{$key}{length} and $f == 1)) {
	# 			$f = 0;
	# 			$true_n{$c} .= $n;
	# #			print FILE "$c\t$true_n{$c}\n";
	# 		}
	# 	}
	# }
	# #close FILE or die $!;
	# my $tn = $c;

	# $result{tn} = $tn;
	# $result{acp} = 0.25*($result{tp}/($result{tp}+$result{fn})+$result{tp}/($result{tp}+$result{fp})+$result{tn}/($result{tn}+$result{fp})+$result{tn}/($result{tn}+$result{fn}));
	# $result{ca} = ($result{acp} - 0.5) * 2;
	 
	################################################


	return %result;
}


##############################################################################
##############################################################################
# Ejecutar como: ./anablast_v01.pl <path_to_fasta> <path_to_blast> <workdir>
##############################################################################
print STDERR <<'HEADER';

                          ____  _           _   
     /\             /\   |  _ \| |         | |  
    /  \   _ __    /  \  | |_) | | __ _ ___| |_ 
   / /\ \ | '_ \  / /\ \ |  _ <| |/ _` / __| __|
  / ____ \| | | |/ ____ \| |_) | | (_| \__ \ |_ 
 /_/    \_\_| |_/_/    \_\____/|_|\__,_|___/\__|
HEADER

if (!$ARGV[2]) {
	print "\n\nPlease execute: ./anablast_v01.pl <path_to_fasta> <path_to_blast> <workdir> <gff> <min_bitscore>\n\n";
	exit;
}

my $fasta = $ARGV[0];
my $blast = $ARGV[1];
my $workdir = $ARGV[2];
my $gff = $ARGV[3];
my $min_bitscore = $ARGV[4];


# Genera el fichero chrom.size necesario para generar los BigWig y 
# devuelve un hash con la secuencia y longitud de cada secuencia dentro del
# fasta sin saltos de línea
my $seq_ref = &parse_fasta($fasta, $workdir);

# Divide y filtra el blast inicial según los parámetros establecidos y
# devuelve las rutas a los diferentes blast
#my %blast_path = &bitscore_filter_frame_split($blast, $workdir, $min_bitscore);


# Obtiene los path de los blast parciales de los ficheros que le indiques
# my %blast_path = (
# 	"1" => "$workdir/frame1.blast",
# 	"2" => "$workdir/frame2.blast",
# 	"3" => "$workdir/frame3.blast",
# 	"-1" => "$workdir/frame-1.blast",
# 	"-2" => "$workdir/frame-2.blast",
# 	"-3" => "$workdir/frame-3.blast",
# 	);


# Obtiene los path de los wig de los ficheros que le indiques
# my @wig_path = (
# "$workdir/frame1.wig",
# "$workdir/frame2.wig",
# "$workdir/frame3.wig",
# "$workdir/frame-1.wig",
# "$workdir/frame-2.wig",
# "$workdir/frame-3.wig",
# );

# Recorre los blast convirtiendolos en wig y bigwig y devuelve la ruta a los wig
#my @wig_path;
#foreach my $key (sort keys %blast_path) {
#	push @wig_path, &frame_to_wig_and_bigwig($key, $blast_path{$key}, $workdir);
#}


# Recorre los wig extrayendo los picos de los diferentes frames
#my %peaks;
#foreach (@wig_path) {
#	my $ref_peaks = &extract_peaks ($_, $seq_ref);
#	%peaks = (%peaks, %{$ref_peaks});
#}

# Mete en memoria los picos de un fichero tsv que le indiques
my $ref_peaks = &tsv_to_peaks("$workdir/peaks.tsv");
my %peaks = %{$ref_peaks};


#&peaks_to_tsv(\%peaks, "$workdir/peaks.tsv");
&peaks_to_gff3(\%peaks, "$workdir/peaks.gff");
#&peaks_to_fasta(\%peaks, "$workdir/peaks.faa");


# busca los ORF y los imprime en un fichero gff
#my %orf = %{&peaks_elong_orf_search(\%peaks, $seq_ref)};
#&orf_to_gff3(\%orf, "$workdir/orf.gff");


# Ejecuta el anotador Sma3s
#&sma3s($workdir, "/home/databases/sma3s/swissprot_2017_7/uniprot_sprot.dat");


# Indexa el gff que le indiques
my $ref_gff = &gff_index($gff);


# Compara con los picos para ver cual es la covertura de los picos con 
# respecto a los item del gff
&peaks_vs_gff ($ref_gff, \%peaks, $workdir);


# Genera un fichero con las estadísticas para cada altura de pico con respecto
# a un gff con los CDS oficiales
open STAT, ">anablast_stat.tsv" or die $!;
print STAT "#Min_bitscore\tTotal_peaks\tTotal_CDS\tMin_Top\tTP\tFP\tFN\tSensitivity\tSpecificity\n";
my %stat;
for my $top (20..200) {
	%stat = &sensitivity_specificity($ref_gff, \%peaks, $top, $seq_ref);
	print STAT "$min_bitscore\t$stat{total}\t$stat{total_cds}\t$top\t$stat{tp}\t$stat{fp}\t$stat{fn}\t$stat{sensitivity}\t$stat{specificity}\n";
}
close STAT or die $!;

exit;