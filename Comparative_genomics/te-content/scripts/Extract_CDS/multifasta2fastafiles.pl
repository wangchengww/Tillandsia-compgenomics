#!/usr/bin/perl


#This tool take a fasta file with many sequences as input and generate fasta file with single sequences

$input = $ARGV[0];


open (INPUT, "$input")||die "Pb d\'ouverture";
while (<INPUT>)
  {
    
    if (($header) = (/^>(.*)$/))
      {
	
	($name) = $header =~ /^(\S+)\s*/;
	close (OUT);
	open (OUT, ">$name.fasta")|| die "impossible de creer le fichier $name.fasta";
	#print STDERR "\nCreating $name.fasta";
	print OUT ">$header\n";
	next;
      }
    if (($seq) = (/^(.*)$/))
      {
	print OUT "$seq\n";
	#print STDERR ".";
      }

  }

close (OUT);
close (INPUT);
