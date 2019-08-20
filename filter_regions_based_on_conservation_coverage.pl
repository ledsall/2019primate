#!/usr/bin/perl

# ################################################################################################################################################################################################################
# filter_regions_based_on_conservation_coverage.pl
#
# CC-BY Lee Edsall
#	email: edsall57@gmail.com
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
#
# Remove regions that aren't conserved across the chimp, gorilla, orangutan, and macaque genomes
#
# Input parameters
#	1. name of input Galaxy file
#	2. coverage threshold (as decimal)
#	3. name for output bed file
#	4. name for output text file
#
# Galaxy file format
#	column 1: chromosome
#	column 2: start
#	column 3: end
#	column 4: genome
#	column 5: number of nucleotides
#	column 6: number of gaps
#
# Output files
#	1. bed file containing regions that passed the filter
#	2. text file containing all of the regions and the coverage for each genome
#
# ################################################################################################################################################################################################################

use strict;
use warnings;

my $coverage_threshold;
my $input_filename;
my $output_filename_bed_file;
my $output_filename_text_file;

my %human_hash;
my %chimp_hash;
my %gorilla_hash;
my %orangutan_hash;
my %macaque_hash;

my $line;
my @line_array;
my $hash_key;
my $checksum;

my $chrom;
my $start;
my $end;
my $minimum_coverage;

my $region_counter_total = 0;
my $region_counter_passed = 0;
my $region_counter_failed = 0;

unless (scalar @ARGV == 4)
{
        print STDERR "Usage: $0: \n\t parameter 1: name of Galaxy file \n\t parameter 2: coverage threshold (as decimal) \n\t parameter 3: name for output bed file \n\t parameter 4: name for output text file \n";
        exit(0);
}

$input_filename = $ARGV[0];
$coverage_threshold = $ARGV[1];
$output_filename_bed_file = $ARGV[2];
$output_filename_text_file = $ARGV[3];

# ###################################################################################################################################################################################################
# Read the galaxy file and set up the hashes
#	The human hash will contain the number of nucleotides
#	The other hashes will be set to zero and updated in a later loop to contain the coverage percentage (as a decimal)
# ###################################################################################################################################################################################################
print STDOUT "\nReading Galaxy file and setting up hashes ... ";

open (INPUT_FILE, "<$input_filename") or die "couldn't open input file $input_filename\n";

while ($line = <INPUT_FILE>)
{
	if ($line =~ /hg19/)
	{
		$region_counter_total++;
		chomp($line);
		@line_array = split(/\t/,$line);
		$hash_key = $line_array[0] . "~" . $line_array[1] . "~" . $line_array[2];
		$human_hash{$hash_key} = $line_array[4];
		$chimp_hash{$hash_key} = 0;
		$gorilla_hash{$hash_key} = 0;
		$orangutan_hash{$hash_key} = 0;
		$macaque_hash{$hash_key} = 0;
	}
}

print STDOUT "done\n";

close INPUT_FILE;

# ###################################################################################################################################################################################################
# Read the galaxy file and update the non-human hashes to contain the coverage percentage
# ###################################################################################################################################################################################################
print STDOUT "Reading Galaxy file and calculating coverage percentage ... ";

open (INPUT_FILE, "<$input_filename") or die "couldn't open input file $input_filename\n";

while ($line = <INPUT_FILE>)
{

	if ($line =~ /hg19/)
	{
		next;
	}

	chomp($line);
	@line_array = split(/\t/,$line);
	$hash_key = $line_array[0] . "~" . $line_array[1] . "~" . $line_array[2];

	if ($line_array[3] =~ /panTro4/)
	{
		$chimp_hash{$hash_key} = sprintf("%.2f", $line_array[4] / $human_hash{$hash_key});
		next;
	}

	if ($line_array[3] =~ /gorGor3/)
	{
		$gorilla_hash{$hash_key} = sprintf("%.2f", $line_array[4] / $human_hash{$hash_key});
		next;
	}

	if ($line_array[3] =~ /ponAbe2/)
	{
		$orangutan_hash{$hash_key} = sprintf("%.2f", $line_array[4] / $human_hash{$hash_key});
		next;
	}

	if ($line_array[3] =~ /rheMac3/)
	{
		$macaque_hash{$hash_key} = sprintf("%.2f", $line_array[4] / $human_hash{$hash_key});
		next;
	}

}

print STDOUT "done\n";

close INPUT_FILE;

# ###################################################################################################################################################################################################
# Go through the hashes and process the regions
#	1. Compare the coverage to the threshold
#	2. Print the coverage information to the text file
#	3. Print regions that passed the threshold to the bed file
# ###################################################################################################################################################################################################

print STDOUT "Filtering the regions and writing the output files (bed file will not be sorted) ... ";

open (OUTPUT_FILE_BED_FILE, ">$output_filename_bed_file") or die "couldn't open output file $output_filename_bed_file\n";
open (OUTPUT_FILE_TEXT_FILE, ">$output_filename_text_file") or die "couldn't open output file $output_filename_text_file\n";

print OUTPUT_FILE_TEXT_FILE "chr" . "\t" . "start" . "\t" . "end" . "\t" . "length" . "\t" . "min_cov" . "\t" . "status" . "\t" . "chimp" . "\t" . "gorilla" . "\t" . "orang" . "\t" . "macaque" . "\n";

foreach $hash_key (keys(%human_hash))
{

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the minimum value
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($chimp_hash{$hash_key} <= $gorilla_hash{$hash_key})
	{
		$minimum_coverage = $chimp_hash{$hash_key};
	}
	else
	{
		$minimum_coverage = $gorilla_hash{$hash_key};
	}

	if ($orangutan_hash{$hash_key} <= $minimum_coverage)
	{
		$minimum_coverage = $orangutan_hash{$hash_key};
	}

	if ($macaque_hash{$hash_key} <= $minimum_coverage)
	{
		$minimum_coverage = $macaque_hash{$hash_key};
	}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Print the first 4 columns to the text file (chromosome, start, end, length)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	($chrom, $start, $end) = split(/~/,$hash_key);

	print OUTPUT_FILE_TEXT_FILE $chrom . "\t" . $start . "\t" . $end . "\t" . $human_hash{$hash_key} . "\t" . $minimum_coverage . "\t";

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Compare the minimum coverage to the threshold
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($minimum_coverage >= $coverage_threshold)
	{
		$region_counter_passed++;
		print OUTPUT_FILE_BED_FILE $chrom . "\t" . $start . "\t" . $end . "\n";
		print OUTPUT_FILE_TEXT_FILE "pass\t";
	}
	else
	{
		$region_counter_failed++;
		print OUTPUT_FILE_TEXT_FILE "fail\t";
	}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Print the remaining columns to the text file (coverage percentages for chimp, gorilla, orangutan, macaque)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	print OUTPUT_FILE_TEXT_FILE $chimp_hash{$hash_key} . "\t" . $gorilla_hash{$hash_key} . "\t" . $orangutan_hash{$hash_key} . "\t" . $macaque_hash{$hash_key} . "\n";

} # END:  foreach $hash_key (keys(%human_hash))

print STDOUT "done\n\n";

close OUTPUT_FILE_BED_FILE;
close OUTPUT_FILE_TEXT_FILE;


# ###################################################################################################################################################################################################
# FINISH UP
# ###################################################################################################################################################################################################

$checksum = $region_counter_total - $region_counter_passed - $region_counter_failed;

print STDOUT "Checksum is $checksum. Checksum should be 0\n\n";

print STDOUT "Parameters\n";
print STDOUT "\tInput Galaxy file\t\t" . $input_filename . "\n";
print STDOUT "\tCoverage threshold\t\t" . $coverage_threshold . "\n";
print STDOUT "\tOutput bed file\t\t\t" . $output_filename_bed_file . "\n";
print STDOUT  "\tOutput text file\t\t" . $output_filename_text_file . "\n\n";

print STDOUT "Region counts\n";
print STDOUT "\tPassed filter\t" . $region_counter_passed . "\n";
print STDOUT "\tFailed filter\t" . $region_counter_failed . "\n";
print STDOUT "\tTotal\t\t" . $region_counter_total . "\n\n";
