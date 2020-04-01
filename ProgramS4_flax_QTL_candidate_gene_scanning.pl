#!/usr/bin/perl -w
########################################################################################################
# Copyright Notice and Disclaimer 
#
# Copyright (c) 2018, Ottawa Research and Development Centre, Agriculture and Agri-Food Canada.
#
# Authors: Dr. Frank M. You
# E mail:  frank.you@canada.ca

########################################################################################################
# This is a pipeline to find all annotated genes around  QTL within upstream and downstream window sizes.
# For example, a QTL is located at 16920407-18739647 on chr 1. For a given 100kb window size, the program
# will scan all genes between 16920407 - 100 kb to 18739647 + 100 kb on chr 1.

# Input:  (1) QTL file: five columns seperated by a tab key with a header line, for example,
#	Trait	QTL	Chr	Coord_start	Coord_end
#	PM	QPM-crc-LG1 	1	16920407	18739647
#	PM	QPM-crc-LG7	7	3817603	3817863
#	PM	QPM-crc-LG9	9	357191	357510

#         (2) Gene annotation file: two gene annotations files have been provided as supplimentary tables
#             Table S5 and Table S6. Any one of them can be used. 
#	  (3) Upstream or downstream window size (bp) from the QTL location (default: 100000 bp):
#         actural region size for gene scan for a QTL = upstream window size + QTL region (from start
#	  coordinate to end coordinate, sometimes the start coordinate may be equal to the end coordinate)
# 	  + downstream window size 
########################################################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_q $opt_d $opt_w);
getopts ('q:d:w:');

my $qtl_file      = $opt_q;
my $gene_annotation_file = $opt_d;
my $window_size = $opt_w?$opt_w:100000;   # default window size for upstream or downstream: 100 Kb, tataling 200Kb

if (!$qtl_file || !$gene_annotation_file || !$window_size) {
    print "Usage: \n";
    print "perl $0 \n";
    print "     -q QTL file\n";
    print "     -d gene annotation file\n";
    print "     -w upstream or downstream window size (bp) (default: 100000 bp) \n";
    exit;
}

&merge_gene_annotations($gene_annotation_file, $qtl_file);	


############################################## End of the mail program #######################################

#Trim function to remove leading and trailing whitespaces
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


sub find_overlapped_genes {
    my ($start, $end, $genes_str) = @_;

    my $is_overlapped = 0;
    my @genes_info = ();
    my @genes_list = split(/@@@@/, $genes_str);
    for (my $i  = 0; $i < @genes_list; $i++) {
	my @cols = split(/\t/, $genes_list[$i]);
	my $start2 = $cols[2];
	my $end2 = $cols[3];

	if (!$start2 || !$end2 || $start2 eq '' || $end2 eq '') {
	    next;
	}		

	$is_overlapped = &is_overlap($start, $end, $start2, $end2);
	if ($is_overlapped) {
	    push(@genes_info, $genes_list[$i]);
	}
    }
	
    return (\@genes_info);
}

sub is_overlap {
    my ($x1, $x2, $y1, $y2) = @_;
    if ($x1 >= $y1 && $x2 <= $y2 ||
	$x1 <= $y1 && $x2 <= $y2 && $x2 > $y1 ||
	$x1 >= $y1 && $x2 >= $y2 && $x1 < $y2 ||
	$x1 <= $y1 && $x2 >= $y2) {
        return 1;
    } else {
        return 0;
    }
}

sub merge_gene_annotations {
    my ($gene_annotation_file, $qtl_file)= @_;
	
    my $out_file = $qtl_file . "_gene_annotations.txt";
    open (OUT, ">$out_file") or die ("Can not open the file $out_file: $!\n");
    open (TABLE, "<$gene_annotation_file") or die ("Can not open the file $gene_annotation_file: $!\n");
	
    my %gene_hash = ();
    my $count = 0;
    my $header = "";
    while (my $line = <TABLE>) {
        chomp($line);
	&trim($line);
	$count++;
	if ($count == 1){
	    $header = $line;
	    next;
	}
		
        my @cols  = split(/\t/, $line);
	my $chr = $cols[0];
		
	if ($gene_hash{$chr}) {
	    $gene_hash{$chr} .= "@@@@" . $line;			
	} else {
	    $gene_hash{$chr} = $line;
	}
    }
    close TABLE;
	
    open (TABLE, "<$qtl_file") or die ("Can not open the file $qtl_file: $!\n");
    $count = 0;
    
    while (my $line = <TABLE>) {
        chomp($line);
	&trim($line);
	$count++;
	
	
        my @cols  = split(/\t/, $line);

	if ($count == 1) {
	    if (@cols != 5) {
		print "Please check the QTL data file: each row must have 5 columns seperated by a tab key ('\t')\n";
		exit(1);
	    }
	    
	    print OUT "$line\t$header\n";
	    next;
	}
	my $chr = $cols[2];
	my $start = $cols[3];
	my $end = $cols[4];

	# Extend QTL 100 kb in both sides				
	$start -= $window_size;
	$end += $window_size;
		
	if ($start < 1) {
	    $start = 1;
	}
				
	if ($gene_hash{$chr}) {
	    my $gene_info = &find_overlapped_genes($start, $end, $gene_hash{$chr});
	    if (@$gene_info > 0) {
		for (my $i = 0; $i < @$gene_info; $i++) {
		    print OUT "$line\t$gene_info->[$i]\n";
		}
	    } else {
		print OUT "$line\n";
	    }
	}  else {
	    print OUT "$line\n";
	}
    }
    close TABLE;
    close OUT;	
}


