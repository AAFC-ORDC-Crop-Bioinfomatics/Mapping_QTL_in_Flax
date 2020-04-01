#!/usr/bin/perl -w
use strict;
# 

########################################################################################################
# Copyright Notice and Disclaimer 
#
# Copyright (c) 2018, Ottawa Research and Development Centre, Agriculture and Agri-Food Canada.
#
# Authors: Dr. Frank M. You
# E mail:  frank.you@canada.ca

########################################################################################################
# This is a pipeline to generate databases for reverse e-PCR
# Input:  (1) genome sequence file (fasta)
#         (2) genome name 
#	  (3) maximum number of sequences for each database (Default: 500)

#  How to generate hash file:
#  (1) Step 1:   famap -tN -b genome.famap fasta_file  
#      (Note: "genome" is any name you can give; the fasta file is 
#      the known reference genome sequence, like flax, rice, Brachypodium)
#  (2) Step 2:   fahash -b genome.hash -w 12 -f3 genome.famap
#      (Note: This command will generate a file called "genome.hash" wich will be used by re-PCR)
# This two steps are implemented in the this script
########################################################################################################

# THE FOLLOWING PACKAGES MUST BE AVAILABLE IN YOUR SERVER MACHINE. PLEASE CHECK YOUR SERVER MACHINE TO
# SEE IF THEY ARE AVAILABLE. IF NOT, YOU NEED INSTALL THEN BEFOR RUN THE APPLICATION

use Getopt::Std;
use vars qw ($opt_i $opt_s $opt_m);
getopts ('i:s:m:');

my $seq_file    = $opt_i;
my $genome_name = $opt_s;
my $max_seqs    = $opt_m || 5000;



if (!$seq_file || !$genome_name) {
	print "Usage: \n";
	print "perl $0 \n";
	print "  -i fasta file name of the reference sequence (* is allowed)\n";
	print "  -s genome name as file name prefix \n";
	print "  -m maximum number of sequences for each database (Default:5000)\n";
    
	exit;
}

# Step 1: split sequence file into several files with maximum number of sequences
my $file_count = &split_seqs($seq_file, $max_seqs);


# Step 2: generate multiple databases

for (my $i = 1; $i <= $file_count; $i++) {
	my $genome_map = $genome_name . "_$i" . ".famap";
	my $genome_hash = $genome_name . "_$i" . ".hash";

	my $cmd = "famap -tN -b $genome_map $seq_file" . "_$i" . ".fas";
	print "$cmd\n";
	system $cmd;

	$cmd = "fahash -b $genome_hash -w 12 -f3 $genome_map";
	print "$cmd\n";
	system $cmd;
}

# end of the main program

sub split_seqs {
    my ($seq_file, $max_seqs) = @_;
    open (SEQ, "<$seq_file") or die ("Can not open the file $seq_file: $!\n");
	
    my $seq = '';
    my $pre_name = '';
    my $name = '';
    my $count = 0;
    my $file_count = 1;
    open (OUT, ">$seq_file" . "_$file_count" . ".fas") or die ("Can not open the file : $!\n");
    while (my $line = <SEQ>) {
	chomp($line);
	$line = &trim($line);
	next if ($line =~ /^\s*$/);
	
	if ($line =~ /^>(\S+)/) {
	    $name = $1;
	    $count++;
	    if ($count > $max_seqs) {
		close OUT;
		$file_count++;
		$count = 1;
		open (OUT, ">$seq_file" . "_$file_count" . ".fas") or die ("Can not open the file : $!\n");
	    }
 	    if ($seq ne '' && $name ne $pre_name && $pre_name ne '') {
		print OUT ">$pre_name\n";
                print OUT &format_seq($seq, 70);
		$seq = '';
	    }
   	    $pre_name = $name;
	} else {
	    $seq .= &trim($line);
	}
    }
    if ($seq ne '') {
	print OUT ">$pre_name\n";
        print OUT &format_seq($seq, 70);
	$seq = '';
    }
    close SEQ;
    close OUT;
    return $file_count;
}


sub format_seq {
    my ($seq, $line_width) = @_;
    my $len = length($seq);
    my $rows = int($len/$line_width);
    $rows++ if ($len % $line_width != 0);
    
    my $seqs = "";
    for (my $i = 0; $i < $rows; $i++) {
        my $offset = $i * $line_width;
        my $width = $line_width;
        $width = ($len - $offset) if ($i == $rows-1); 
        $seqs .= uc(substr($seq, $offset, $line_width)) . "\n";
    }
    
    return $seqs;
}

# Trim function to remove leading and trailing whitespaces
sub trim {
    my $string = shift;
    if ($string =~ /^\s+/) {
	$string =~ s/^\s+//;
    }
    if ($string =~ /\s+$/) {
	$string =~ s/\s+$//;
    }
    return $string;
}