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
	# ejemplo: &bitscore_filter(<path_to_blast> <workdir> <bit_threshold> <ref_seq>)
	################################################################################

	my ($blast, $workdir, $bit_threshold, $ref_seq) = @_;
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
		my ($seqid, $chrom, $start, $end, $evalue, $bit, $frame) = split (/\t/, $line);
		my $f = 0;
		for my $key (keys %{$ref_seq}) {
			if ($key eq $chrom) {
				$f = 1;
				last;
			}
		}
		if ($f == 0) {
			print STDERR "Blast warning: no chrom id ($chrom) in seq file\n";
			next;
		}
		if ($bit >= $bit_threshold) {
			foreach (@frame) {			
				if ($frame == $_) {
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
		print STDERR "No blast hit with bitscore higher than $bit_threshold\n"
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

	my ($frame, $blast, $workdir) = @_;
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
	# ejemplo: &extract_pic(<path_to_wig>, <sequence_reference>, <blast_index_reference>)
	################################################################################
	my ($wig, $seq_ref) = @_; # Se ha eliminado el tercer argumento $blat_ref para eliminar el cáculo de pvalor
	my %seq = %{$seq_ref};

	print STDERR "Extracting peaks from WIG: $wig\n";
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
				# $pic{"frame_$frame"."_peak_$c"}{pvalue} = &pvalue_calc($chrom, $frame, $start, $end, $blast_ref); # comentado para eliminar el calculo de pvalor
				# my $p = $pic{"frame_$frame"."_peak_$c"}{pvalue};
				# print "$chrom:$start..$end\t$p\n";
				# exit if ($p == 1);
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
				# $pic{"frame_$frame"."_peak_$c"}{pvalue} = &pvalue_calc($chrom, $frame, $start, $end, $blast_ref); # comentado para eliminar el calculo de pvalor
				# my $p = $pic{"frame_$frame"."_peak_$c"}{pvalue};
				# print "$chrom:$start..$end\t$p\n";
				# exit if ($p == 1);
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
		# $pic{"frame_$frame"."_peak_$c"}{pvalue} = &pvalue_calc($chrom, $frame, $start, $end, $blast_ref); Comentado para eliminar el calculo del pvalor
		# my $p = $pic{"frame_$frame"."_peak_$c"}{pvalue};
		# print "$chrom:$start..$end\t$p\n";
		# exit if ($p == 1);
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

sub pvalue_calc {
	my ($chr, $frame, $start, $end, $ref_blast_index) = @_;
	my %blast_index = %{$ref_blast_index};
	my $pvalue = 1;

	my $start_range = 10000 * (int($start/10000) + 1);
	my $end_range = 10000 * (int($end/10000) + 1);
	my %nn;
	my $length = $end - $start + 1;
	for my $n ($start..$end) {
		$nn{$n} = 1;
	}
	for (my $i = $start_range; $i <= $end_range; $i += 10000) {
		foreach my $k (keys %{$blast_index{$chr}{$frame}{$i}}) {
			my $b_uniref = $blast_index{$chr}{$frame}{$i}{$k}{uniref};
			my $b_evalue = $blast_index{$chr}{$frame}{$i}{$k}{evalue};
			my $b_start = $blast_index{$chr}{$frame}{$i}{$k}{start};
			my $b_end = $blast_index{$chr}{$frame}{$i}{$k}{end};
			my $b_pvalue = $blast_index{$chr}{$frame}{$i}{$k}{pvalue};
			my $b_length = $blast_index{$chr}{$frame}{$i}{$k}{length};
			next if ($start > $b_end or $end < $b_start);
			my $overlap;
			for my $n ($b_start..$b_end) {
				if ($nn{$n}) {
					$overlap++;
				}
			}
			next unless ($overlap/$length >= 0.2 or $overlap/$b_length >= 0.2);
			#print "\t$b_uniref\t$b_start\t$b_end\t$b_evalue\t$b_pvalue\n";
			$pvalue *= $b_pvalue;
		}
	}
	return $pvalue;
}

sub blast_index {
	################################################################################
	# Parte de un fichero en formato blast y lo indexa en memoria para poder trabajar
	# con el devolviendo la referencia de un hash de hashes
	# ejemplo: &gff_index(<path to blast>)
	################################################################################
	my ($path_to_blast) = @_;

	print STDERR "\nIndexing blast $path_to_blast\n";
	my %hash;
	my $c;
	open BLAST, $path_to_blast or die $!;
	while (<BLAST>) {
		chomp;
		$c++;
		next if ($_ =~ m/^#/ or $_ =~ m/^$/);
		my ($uniref, $chr, $start, $end, $evalue, $bit, $frame) = split (/\t/,$_);
		if ($start > $end) {
			my $s_temp = $start;
			my $e_temp = $end;
			$start = $e_temp;
			$end = $s_temp;
		}

		my $e = 2.718281828459045235360;
		my $exp= $e**(-$evalue);
		my $pvalue= 1-$exp;

		my $start_range = 10000 * (int($start/10000) + 1);
		my $end_range = 10000 * (int($end/10000) + 1);

		for (my $i = $start_range; $i <= $end_range; $i += 10000) {
			$hash{$chr}{$frame}{$i}{$c}{evalue} = $evalue;
			$hash{$chr}{$frame}{$i}{$c}{uniref} = $uniref;
			$hash{$chr}{$frame}{$i}{$c}{chr} = $chr;
			$hash{$chr}{$frame}{$i}{$c}{start} = $start;
			$hash{$chr}{$frame}{$i}{$c}{end} = $end;
			$hash{$chr}{$frame}{$i}{$c}{frame} = $frame;
			$hash{$chr}{$frame}{$i}{$c}{pvalue} = $pvalue;
			$hash{$chr}{$frame}{$i}{$c}{length} = $end - $start + 1;
		}
	}
	close BLAST or die $!;
	return \%hash;
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
			chomp @header;
			next;
		}
		my @data = split(/\t/, $_);
		chomp @data;

		foreach my $i (1..$#header) {
			$peaks{$data[0]}{$header[$i]} = $data[$i];
		}
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
			my $f = 0;
			for (my $i = 0; $i <= length($nn); $i += 3) {
				my $codon = substr($nn, $i, 3);
				if ($codon eq "ATG" or $codon eq "ATA" or $codon eq "ATC" or $codon eq "GTG" or $codon eq "TTG" or $codon eq "CTG" or $codon eq "ATT") {
					if ($f == 0) {
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

sub new_orf_search {
	################################################################################
	# Busca ORF en los picos con un núevo método
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
		
		my $nn; # Secuencia de nucleótidos extraida
		my $frame_start; # posición de inicio del frame 
		my $rest; # resto de la posición de inicio cuando lo dividimos entre 3
		my $objetiv_rest; # resto que quiero tener para trabajar con el frame correcto
		my $frame_end; # posición final del frame

		if ($peaks{$key}{frame} > 0) { # Para los frames positivos
			# cálculo de la posición de inicio del frame
			if ($peaks{$key}{frame} == 1) {
				$objetiv_rest = 0;
			}
			elsif ($peaks{$key}{frame} == 2) {
				$objetiv_rest = 1;
			}
			elsif ($peaks{$key}{frame} == 3) {
				$objetiv_rest = 3;
			}

			$frame_start = $peaks{$key}{start}; # iniciamos el frame en la posición de inicio del pico
			$rest = $peaks{$key}{start} % 3;

			for (1..6) { # Corregimos el inicio hasta tener el resto adecuado
				if ($rest != $objetiv_rest) {
					$frame_start++;
				}
				elsif ($rest == $objetiv_rest) {
					last;
				}
			}

			$frame_end = $frame_start + 3; # iniciamos el final del frame en el siguiente codon al inicio del frame
			
			# búsqueda del primer codon de stop previo al inicio del frame
			while (1) {
				if ($frame_start <= 1) {
					$frame_start = $peaks{$key}{frame};
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $frame_start -1, 3);
				if ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG" or $codon eq "NNN") {
					last;
				}
				else {
					$frame_start -= 3;
				}
			}

			# Búsqueda del siguiente codon de stop 
			while (1) {
				if ($frame_end >= $seq{$peaks{$key}{chrom}}{length} - 3) {
					$frame_end = $seq{$peaks{$key}{chrom}}{length};
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $frame_end -1, 3);
				if ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG" or $codon eq "NNN") {
					$frame_end += 3;
					if ($frame_end < $peaks{$key}{end}) { # Si encuentra un codon de stop comprueba que ha superado el final del pico
						next;
					}
					else {
						last;	
					}
				}
				else {
					$frame_end += 3;
				}
			}

			$nn = substr ($seq{$peaks{$key}{chrom}}{nn}, $frame_start -1, $frame_end - $frame_start +1); # Extrae la secuencia de nucleótidos
			my $f = 0;
			for (my $i = 0; $i <= length($nn); $i += 3) {
				my $codon = substr($nn, $i, 3);
				if ($codon eq "ATG" or $codon eq "ATA" or $codon eq "ATC" or $codon eq "GTG" or $codon eq "TTG" or $codon eq "CTG" or $codon eq "ATT") {
					if ($f == 0) {
						$c++;
						$f = 1;
						$orf{"$key"."_ORF_$c"}{start} = $frame_start + $i;
						$orf{"$key"."_ORF_$c"}{parent} = $key;
						$orf{"$key"."_ORF_$c"}{chrom} = $peaks{$key}{chrom};
						$orf{"$key"."_ORF_$c"}{strand} = $peaks{$key}{strand};
						$orf{"$key"."_ORF_$c"}{frame} = $peaks{$key}{frame};
					}
				}
				elsif ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG") {
					if ($f == 1) {
						$f = 0;
						$orf{"$key"."_ORF_$c"}{end} = $frame_start + $i +2;
						$orf{"$key"."_ORF_$c"}{nn} = substr ($seq{$peaks{$key}{chrom}}{nn}, $orf{"$key"."_ORF_$c"}{start} -1, $orf{"$key"."_ORF_$c"}{end} - $orf{"$key"."_ORF_$c"}{start} + 1);
						$orf{"$key"."_ORF_$c"}{aa} = &dna2aa($orf{"$key"."_ORF_$c"}{nn}, 1);
						delete $orf{"$key"."_ORF_$c"} if (length $orf{"$key"."_ORF_$c"}{aa} <= 5);
					}
				}
			}
		}
		else { # para los frames negativos

			# cálculo de la posición de inicio del frame
			if ($peaks{$key}{frame} == -1) {
				$objetiv_rest = 0;
			}
			elsif ($peaks{$key}{frame} == -2) {
				$objetiv_rest = 1;
			}
			elsif ($peaks{$key}{frame} == -3) {
				$objetiv_rest = 3;
			}

			$frame_start = $peaks{$key}{end};
			$rest = ($seq{$peaks{$key}{chrom}}{length} - $frame_start) % 3;

			for (1..6) { # Corregimos el inicio hasta tener el resto adecuado
				if ($rest != $objetiv_rest) {
					$frame_start--;
				}
				elsif ($rest == $objetiv_rest) {
					last;
				}
			}

			$frame_end = $frame_start - 3; # Iniciamos el final del frame en el codon anterior

			# Buscamos el primer codon de fin despues de la posición de inicio
			while (1) {
				if ($frame_start >= $seq{$peaks{$key}{chrom}}{length} - 3) {
					$frame_start = $seq{$peaks{$key}{chrom}}{length} + $peaks{$key}{frame};
					last;
				}
				my $codon = substr ($seq{$peaks{$key}{chrom}}{nn}, $frame_start -1, 3);
				my $rcodon = &reverse_comp($codon);
				if ($rcodon eq "TGA" or $rcodon eq "TAA" or $rcodon eq "TAG" or $rcodon eq "NNN") {
					last;
				}
				else {
					$frame_start -= 3;
				}
			}

			# buscamos el primer codon de inicio despues del final del pico
			while (1) {
				if ($frame_end <= 1) {
					$frame_end = abs($peaks{$key}{frame});
					last;
				}

				my $codon = substr($seq{$peaks{$key}{chrom}}{nn}, $frame_end +1, -3);
				my $rcodon = &reverse_comp($codon);

				if ($rcodon eq "TGA" or $rcodon eq "TAA" or $rcodon eq "TAG" or $rcodon eq "NNN") {
					$frame_end -= 3;
					if ($frame_end > $peaks{$key}{start}) {
						next;
					}
					else {
						last;	
					}
					
				}
				else {
					$frame_end -= 3;
				}
			}
			$nn = substr ($seq{$peaks{$key}{chrom}}{nn}, $frame_start +1, $frame_start - $frame_end +1);
			#print "\n$peaks{$key}{chrom}:$start..$end\t$peaks{$key}{frame}\n$nn\n";
			$nn = &reverse_comp($nn);
			
			my $f = 0;

			for (my $i = 0; $i <= length($nn); $i += 3) {
				my $codon = substr($nn, $i, 3);
				if ($codon eq "ATG" or $codon eq "ATA" or $codon eq "ATC" or $codon eq "GTG" or $codon eq "TTG" or $codon eq "CTG" or $codon eq "ATT") {
					if ($f == 0) {
						$c++;
						$f = 1;
						$orf{"$key"."_ORF_$c"}{end} = $frame_start - $i;
						$orf{"$key"."_ORF_$c"}{parent} = $key;
						$orf{"$key"."_ORF_$c"}{chrom} = $peaks{$key}{chrom};
						$orf{"$key"."_ORF_$c"}{strand} = $peaks{$key}{strand};
						$orf{"$key"."_ORF_$c"}{frame} = $peaks{$key}{frame};
					}
				}
				elsif ($codon eq "TGA" or $codon eq "TAA" or $codon eq "TAG") {
					if ($f == 1) {
						$f = 0;
						$orf{"$key"."_ORF_$c"}{start} = $frame_start - $i - 2;
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
	print STDERR "Total $c item in $gff\n";
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

	my ($ref_gff, $ref_peaks, $min_top, $seq_ref, $print) = @_;
	print STDERR "\nSensitivity and Specificity (min. top: $min_top)\n";
	
	my %chrom = %{$seq_ref};
	my %peaks = %{$ref_peaks};
	my %gff = %{$ref_gff};

	my %cds;
	my %peak_cov;
	my %result;
	
	# Extraer los CDS a partir del gff indexado
	foreach my $k (keys %gff) {
		next if ($k eq "sources" or $k eq "total");
		#print "$gff{$k}\n";
		foreach my $k2 (keys %{$gff{$k}}) {
			foreach my $k3 (keys %{$gff{$k}{$k2}}) {
				foreach my $k4 (keys %{$gff{$k}{$k2}{$k3}}) {
					$cds{"$k:$gff{$k}{$k2}{$k3}{$k4}{start}..$gff{$k}{$k2}{$k3}{$k4}{end}"} = 0;
				}
			}
		}
	}

	$result{total_cds} = keys %cds;
	#print STDERR "Total uniq CDS = $result{total_cds}\n";

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
		
		$peak_cov{$key} = 0; # Asigna 0 a la covertura del pico

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
			next if ($peak_cov{$key} == 1); # si el pico tiene CDS en su cadena no comprueba la complementaria

			my $rev_strand;
			if ($p_strand eq "+") {
				$rev_strand = "-";
			}
			else {
				$rev_strand = "+";
			}

			# Recorre los cds de un rango determinado en la cadena reversa del pico
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
		}
	}

	print STDERR "Total peaks:\t$result{total}\n";

	$result{fp} = 0;
	foreach my $key2 (keys %peaks) {
		next if ($peaks{$key2}{top} < $min_top);
		if ($peak_cov{$key2} == 1) {
			#$result{tp_peaks}++;
		}
		else {
			$result{fp}++;
		}
	}
	
	$result{tp} = 0;
	$result{fn} = 0;
	open FILE, ">cds_$min_top.txt" or die $!;
	foreach my $key3 (keys %cds) {
		if ($cds{$key3} == 1) {
			print FILE "$key3\t$min_top\t$cds{$key3}\n"; #CDS con picos
			$result{tp}++;
		}
		else {
			print FILE "$key3\t$min_top\t$cds{$key3}\n"; #CDS sin picos
			$result{fn}++;
		}
	}
	close FILE or die $!;

	# Imprime en un fichero de texto los picos que solapan o no solapan con CDS y la altura del pico solo si se le indica que lo imprima.
	if ($print eq "print") {
		open FILE, ">peaks_$min_top.txt" or die $!;
		foreach my $key (keys %peaks) {
			next if ($peaks{$key}{top} < $min_top);
			print FILE "$peaks{$key}{chrom}:$peaks{$key}{start}..$peaks{$key}{end}\t$peaks{$key}{top}\t$peak_cov{$key}\n"; #Picos con y sin CDS
		}
		close FILE or die $!;
	}

	$result{sensitivity} = $result{tp}/($result{tp} + $result{fn});
	$result{specificity} = $result{tp}/($result{tp} + $result{fp});

	%cds = ();
	%peak_cov = ();
	%peaks = ();
	%gff = ();

	return %result;
}


##############################################################################
##############################################################################
## Ejecutar como: ./anablast_v01.pl <path_to_fasta> <path_to_blast> <workdir> 
## <gff> <min_bitscore>
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
	print "\n\nPlease execute: ./anablast_v01.pl <path_to_fasta> <path_to_blast> <workdir> <gff> <min_bitscore> <cds_gff>\n\n";
	exit;
}

my ($fasta, $blast, $workdir, $gff, $min_bitscore, $gff_cds) = @ARGV;

print "\n#########\nparameters:\nfasta = $fasta\nblast =$blast\nworkdir = $workdir\nannot gff = $gff\nmin biscore = $min_bitscore\ncds gff = $gff_cds\n##########\n";

## Genera el fichero chrom.size necesario para generar los BigWig y 
## devuelve un hash el identificador, la secuencia y longitud de cada
## secuencia dentro del fasta sin saltos de línea
my $seq_ref = &parse_fasta($fasta, $workdir);

## Divide y filtra el blast inicial según el bit score que le indiques y
## devuelve las rutas a los diferentes blast
# my %blast_path = &bitscore_filter_frame_split($blast, $workdir, $min_bitscore, $seq_ref);

## Obtiene los path de los blast parciales de los ficheros que le indiques
# my %blast_path = (
# 	"1" => "$workdir/frame1.blast",
# 	"2" => "$workdir/frame2.blast",
# 	"3" => "$workdir/frame3.blast",
# 	"-1" => "$workdir/frame-1.blast",
# 	"-2" => "$workdir/frame-2.blast",
# 	"-3" => "$workdir/frame-3.blast",
# 	);

## Recorre los blast parciales convirtiendolos en wig y bigwig y devuelve
## la ruta a los wig
# my @wig_path;
# foreach my $key (sort keys %blast_path) {
# 	push @wig_path, &frame_to_wig_and_bigwig($key, $blast_path{$key}, $workdir);
# }

## Obtiene los path de los wig de los ficheros que le indiques
# my @wig_path = (
# "$workdir/frame1.wig",
# "$workdir/frame2.wig",
# "$workdir/frame3.wig",
# "$workdir/frame-1.wig",
# "$workdir/frame-2.wig",
# "$workdir/frame-3.wig",
# );

## Recorre los wig extrayendo los picos de los diferentes frames.
## Por defecto empieza a considerarlo pico a partir de una altura de 15,
## solo cuando llega hasta 20 lo considera pico y termina el pico cuando 
## vuelve a llegar a 15
# my %peaks;
# foreach (@wig_path) {
# 	my $wig = $_;
# 	my $blast = $wig;
# 	$blast =~ s/\.wig/\.blast/;
# 	# my $blast_ref_index = &blast_index($blast); # comentado para eliminar el calculo de pavalor 
# 	my $ref_peaks = &extract_peaks ($wig, $seq_ref); # Se ha elimindo el tercer argumento $blast_ref_index para eliminar el pvalor
# 	%peaks = (%peaks, %{$ref_peaks});
# }

## Mete en memoria los picos guardados en un fichero tsv que le indiques
my $ref_peaks = &tsv_to_peaks("$workdir/peaks.tsv");
my %peaks = %{$ref_peaks};

## Guarda los picos en un fichero tsv
&peaks_to_tsv(\%peaks, "$workdir/peaks.tsv");

## Guarda los picos en un fichero gff3
&peaks_to_gff3(\%peaks, "$workdir/peaks.gff");

## Guarda las secuencias de los picos en un fichero fasta
#&peaks_to_fasta(\%peaks, "$workdir/peaks.faa");

## Busca los ORF y los guarda en un fichero gff3
my %orf = %{&peaks_elong_orf_search(\%peaks, $seq_ref)};
&orf_to_gff3(\%orf, "$workdir/orf.gff");

## Ejecuta el anotador Sma3s
#&sma3s($workdir, "/home/databases/sma3s/swissprot_2017_7/uniprot_sprot.dat");

## Indexa el gff que le indiques en el argumento 4 ($ARGV[3]). Normalmente
## es un gff de anotaciones del genoma para compararlo con los picos
# my $ref_gff = &gff_index($gff);

## Compara con los picos para ver cual es la covertura de los picos con 
## respecto a los item del gff y guarda los resultados en el fichero 
## peaks_vs_gff.tsv
# &peaks_vs_gff ($ref_gff, \%peaks, $workdir);

## Genera un fichero con las estadísticas para cada altura de pico con
## respecto a un gff con los CDS oficiales
# my $ref_cds_gff = &gff_index("$gff_cds");

# open STAT, ">anablast_stat.tsv" or die $!;
# print STAT "Min_bitscore\tTotal_peaks\tTotal_CDS\tMin_Top\tTP\tFP\tFN\tSensitivity\tSpecificity\n";
# my %stat;
# for my $top (20..200) {
# 	my $print = "no print";
# 	if ($top == 20) {
# 		$print = "print";
# 	}
# 	%stat = &sensitivity_specificity($ref_cds_gff, \%peaks, $top, $seq_ref, $print);
# 	print STAT "$min_bitscore\t$stat{total}\t$stat{total_cds}\t$top\t$stat{tp}\t$stat{fp}\t$stat{fn}\t$stat{sensitivity}\t$stat{specificity}\n";
# }
# close STAT or die $!;

exit;
