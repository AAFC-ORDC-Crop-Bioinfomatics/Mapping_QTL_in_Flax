#!/usr/bin/perl -w
########################################################################################################
# Copyright Notice and Disclaimer 
#
# Copyright (c) 2018, Ottawa Research and Development Centre, Agriculture and Agri-Food Canada.
#
# Authors: Dr. Frank M. You
# E mail:  frank.you@canada.ca

########################################################################################################
# This is a program to convert coordinates of markers on scaffolds based reference sequences to the
# recently released flax chromosome-scale pseudomolecules (You et al. 2018, Chromosome-scale pseudomolecules
# refined by optical, physical and genetic maps in flax. The Plant Journal, 95,371–384)

# Input:
#       (1) a lookup file for scaffold to pseudomolecule mapping: Table S4 must be used.

#       (2) A file with marker coordinates on flax scaffolds based reference sequences, 
#         which must have three columns seperated by a tab key - marker names and marker coodinates \n";
#	Marker	Scaffold_ID	Coordinate_on_scaffold
#	scaffold112_114241	scaffold112	114241
#	scaffold1491_318496	scaffold1491	318496
#	scaffold31_1800846	scaffold31	1800846
#	scaffold344_309662	scaffold344	309662
#	scaffold51_1349321	scaffold51	1349321
#	scaffold59_572553	scaffold59	572553
#	scaffold156_641874	scaffold156	641874
#	scaffold147_367986	scaffold147	367986
#	scaffold859_123972	scaffold859	123972
#	scaffold297_275113	scaffold297	275113
#	scaffold361_14957	scaffold361	14957
#	scaffold273_68457	scaffold273	68457

# Output:
#	Marker	Scaffold_ID	Coordinate_on_scaffold	Chr	New_Chr_Coord
#	scaffold112_114241	scaffold112	114241	1	18444086
#	scaffold1491_318496	scaffold1491	318496	6	14006651
#	scaffold31_1800846	scaffold31	1800846	3	3929932
#	scaffold344_309662	scaffold344	309662	1	11008279
#	scaffold51_1349321	scaffold51	1349321	4	10532424
#	scaffold59_572553	scaffold59	572553	1	10051709
#	scaffold156_641874	scaffold156	641874	3	5906791
#	scaffold147_367986	scaffold147	367986	5	11288517
#	scaffold859_123972	scaffold859	123972	15	1939372
#	scaffold297_275113	scaffold297	275113	1	16435852
#	scaffold361_14957	scaffold361	14957	1	16726904
#	scaffold273_68457	scaffold273	68457	8	585113

########################################################################################################

use strict;
use Getopt::Std;
use vars qw ($opt_m $opt_d);
getopts ('m:d:');


my $mapping_file                  = $opt_m;
my $scaffold_coordinate_file      = $opt_d;

if (!$mapping_file || !$scaffold_coordinate_file) {
    print "Usage: \n";
    print "perl $0 \n";
    print "     -m scaffold to pseudomolecule mapping file. Table S4 must be used.\n";
    print "     -d scaffold coordinate data file which must have three column: marker name, scaffold IDs and coordinates \n";
    exit;
}

# read mapping file
my $scaffold_coords = &read_mapping_data($mapping_file);

# convert SNP coordinates from scaffolds to chromosome/linkage group
&convert_scaffold_coordinates($scaffold_coordinate_file, $scaffold_coords);

############################################## End of the mail program #######################################

sub read_mapping_data {
    my $mapping_file = shift;
	print "Read mapping data ...\n";
	
    my %scaffold_coords = ();
    open (TABLE, "<$mapping_file") or die ("Can not open the file $mapping_file: $!\n");
    
    my $count = 0;
    while (my $line = <TABLE>) {
	$count++;
        next if ($count == 1);  # skip header

        chomp($line);
        my @cols  = split(/\t/, $line);
        my $scaffold_name               = $cols[0]; 
        my $new_scaffold_name           = $cols[1]; 
        my $scaffold_start              = $cols[2]; 
        my $scaffold_end                = $cols[3]; 
        my $chr                         = $cols[4];
        my $chr_start                   = $cols[5];
        my $chr_end                     = $cols[6];
        my $scaffold_orientation        = $cols[7];
        

        my $coords = $scaffold_start . "&" . $scaffold_end . "&" . $chr_start . "&" . $chr_end . "&" . $chr . "&" . $new_scaffold_name  . "&" . $scaffold_orientation;
        if (exists($scaffold_coords{$scaffold_name})) {
            $scaffold_coords{$scaffold_name} .= "/" . $coords;
        } else {
            $scaffold_coords{$scaffold_name} = $coords;
        }
    }    
    close TABLE;
    
    return \%scaffold_coords;
 }


sub convert_scaffold_coordinates {
    my ($snp_data_file, $scaffold_coords) = @_;
 
    print "Convert SNP coordinates ...\n";
 
    open (TABLE, "<$snp_data_file") or die ("Can not open the file $snp_data_file: $!\n");

    my $out_file = $snp_data_file . ".converted.txt";
    open (OUT, ">$out_file") or die ("Can not open the file $out_file: $!\n");


    my $count = 0;
    my $found_count = 0;
    my $removed_count = 0;
    my $low_maf_count = 0;
    my $row_count = 0;
    my $converted_count = 0;
    my $total_count = 0;
    while (my $line = <TABLE>) {
        chomp($line);
	
	$row_count++;
	if ($row_count == 1) {
	    $line =~ s/\-/_/g;
	    print OUT "$line\tChr\tNew_Chr_Coord\n";   
	    next;
	}
	next if ($line =~ /^\s*$/);
        $total_count++;
	
	my @cols  = split(/\t/, $line);
        my $marker_name             = $cols[0]; 
        my $scaffold_name           = $cols[1]; 
        my $scaffold_coord          = $cols[2];
        if (!$scaffold_coords->{$scaffold_name}) {
            print OUT "$line\n";
            $count++;
            next;
        }
        
        my @coord_arr = split(/\//, $scaffold_coords->{$scaffold_name});
        my $found = 0;
        for (my $i = 0; $i < @coord_arr; $i++) {
            my @coord_arr2 = split(/&/, $coord_arr[$i]);
            
            my $scaffold_start      = $coord_arr2[0];
            my $scaffold_end        = $coord_arr2[1];
            my $chr_start           = $coord_arr2[2];
            my $chr_end             = $coord_arr2[3];
            my $chr                 = $coord_arr2[4];
            my $new_scaffold_name   = $coord_arr2[5];
	    my $orientation         = $coord_arr2[6];

            
            if ($scaffold_coord >= $scaffold_start && $scaffold_coord <= $scaffold_end) {

                my $new_coord = &get_new_coord($scaffold_coord, $scaffold_start, $chr_start, $chr_end, $orientation);  # start from 1
		print OUT "$marker_name\t$scaffold_name\t$scaffold_coord\t$chr\t$new_coord\n";
		$converted_count++;
		$found = 1;
                last;
            }
        }
        
        if (!$found) {
	    print OUT "$line\n";
        }
    }    
    close TABLE;
    close OUT;
    
    print "A total of $converted_count from $total_count markers have been sucessfully converted. The results is saved in the file \"$out_file\".\n";
}

sub get_new_coord {
    my ($scaffold_snp_coord, $scaffold_start, $chr_start, $chr_end, $orientation) = @_;
    my $new_snp_coord = 0;
    
    my $c1 = $scaffold_snp_coord - $scaffold_start;
    
    if ($orientation eq '+') {
        $new_snp_coord = $chr_start + $c1;
    } else {
	my $tmp = $chr_start;
	$chr_start = $chr_end;
	$chr_end = $tmp;
        $new_snp_coord = $chr_start - $c1;
    }
    
    return $new_snp_coord;
}

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

