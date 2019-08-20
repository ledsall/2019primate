#!/usr/bin/perl

# ################################################################################################################################################################################################################
# zero_count_filter_HCGOM_min_2_replicates.pl
#
# CC-BY Lee Edsall
#	email: edsall57@gmail.com
#	Twitter: @LeeEdsall
#
# This script was used to analyze data for Edsall et al. 2019 which compared DNase-seq data from 5 primates (human, chimpanzee, gorilla, orangutan, macaque)
#
# Remove entries from a text file that have a score of zero for 2 or 3 replicates in a species
# Hard coded for 3 replicates per species in this order: human, chimp, gorilla, orangutan, macaque
#
# Columns are:
#	column 1	chromosone
#	column 2	start
#	column 3	end
#	column 4	human b1t2
#	column 5	human b2t1
#	column 6	human b3t1
#	column 7	chimp b1t2
#	column 8	chimp b2t2
#	column 9	chimp b3t1
#	column 10	gorilla b1t1
#	column 11	gorilla b2t1
#	column 12	gorilla b3t1
#	column 13	orangutan b1t1
#	column 14	orangutan b2t1
#	column 15	orangutan b3t1
#	column 16	macaque b1ta
#	column 17	macaque b2ta
#	column 18	macaque b3t1
#
# print stats to STDOUT and produce 2 output files
#	1. filtered text file (all non-zero entries)
#	2. text file containing entries that were removed
#
# 2 input parameters
#	1. input file
#	2. base name for output files
# ################################################################################################################################################################################################################

use strict;
use warnings;

my $input_filename;
my $output_filename_base;
my $output_filename_kept;
my $output_filename_removed;

my $line;
my @line_array;

my $checksum1;
my $checksum2;
my $checksum3;
my $checksum4;
my $checksum5;
my $checksum6;

my $human_flag;
my $chimp_flag;
my $gorilla_flag;
my $orangutan_flag;
my $macaque_flag;
my $flag_counter;

my $entry_counter_input_file = 0;
my $entry_counter_kept = 0;

my $entry_counter_removed_H = 0;
my $entry_counter_removed_C = 0;
my $entry_counter_removed_G = 0;
my $entry_counter_removed_O = 0;
my $entry_counter_removed_M = 0;

my $entry_counter_removed_HC = 0;
my $entry_counter_removed_HG = 0;
my $entry_counter_removed_HO = 0;
my $entry_counter_removed_HM = 0;
my $entry_counter_removed_CG = 0;
my $entry_counter_removed_CO = 0;
my $entry_counter_removed_CM = 0;
my $entry_counter_removed_GO = 0;
my $entry_counter_removed_GM = 0;
my $entry_counter_removed_OM = 0;

my $entry_counter_removed_HCG = 0;
my $entry_counter_removed_HCO = 0;
my $entry_counter_removed_HCM = 0;
my $entry_counter_removed_HGO = 0;
my $entry_counter_removed_HGM = 0;
my $entry_counter_removed_HOM = 0;
my $entry_counter_removed_CGO = 0;
my $entry_counter_removed_CGM = 0;
my $entry_counter_removed_COM = 0;
my $entry_counter_removed_GOM = 0;

my $entry_counter_removed_HCGO = 0;	# All except macaque
my $entry_counter_removed_HCGM = 0;	# All except orangutan
my $entry_counter_removed_HCOM = 0;	# All except gorilla
my $entry_counter_removed_HGOM = 0;	# All except chimp
my $entry_counter_removed_CGOM = 0;	# All except human

my $entry_counter_removed_1_species = 0;
my $entry_counter_removed_2_species = 0;
my $entry_counter_removed_3_species = 0;
my $entry_counter_removed_4_species = 0;
my $entry_counter_removed_total = 0;

# =================================================================================================================================================================
# Make sure the necessary parameters are present
# =================================================================================================================================================================
unless (scalar @ARGV == 2)
{
        print STDERR "Usage: $0: \n\tparameter 1: input file containing scores for 3 replicates per species in this order: human, chimp, gorilla, orangutan, macaque\n\tparameter 2: base name for output files\n";
        exit(56);
}

# =================================================================================================================================================================
# Open files
# =================================================================================================================================================================

$input_filename = $ARGV[0];
$output_filename_base = $ARGV[1];

$output_filename_kept = $output_filename_base . ".zero_filtered.txt";
$output_filename_removed = $output_filename_base . ".zero_entries.txt";

open (INPUT_FILE, "<$input_filename") or die "couldn't open input file $input_filename\n";
open (OUTPUT_FILE_KEPT, ">$output_filename_kept") or die "couldn't open output file $output_filename_kept\n";
open (OUTPUT_FILE_REMOVED, ">$output_filename_removed") or die "couldn't open output file $output_filename_removed\n";

# =================================================================================================================================================================
# Main loop
# =================================================================================================================================================================

print STDOUT "\nReading input file ($input_filename) and removing DHS sites that have zero counts for 2 or more replicates ... ";

while ($line = <INPUT_FILE>)
{
	$entry_counter_input_file++;
	chomp($line);
	@line_array = split(/\t/,$line);

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Set the "remove entries?" flag to zero
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	$human_flag = 0;
	$chimp_flag = 0;
	$gorilla_flag = 0;
	$orangutan_flag = 0;
	$macaque_flag = 0;

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Check the human replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( ($line_array[3] == 0 && $line_array[4] == 0) || ($line_array[3] == 0 && $line_array[5] == 0) || ($line_array[4] == 0 && $line_array[5] == 0) )
	{
		$human_flag = 1;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Check the chimp replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( ($line_array[6] == 0 && $line_array[7] == 0) || ($line_array[6] == 0 && $line_array[8] == 0) || ($line_array[7] == 0 && $line_array[8] == 0) )
	{
		$chimp_flag = 1;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Check the gorilla replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( ($line_array[9] == 0 && $line_array[10] == 0) || ($line_array[9] == 0 && $line_array[11] == 0) || ($line_array[10] == 0 && $line_array[11] == 0) )
	{
		$gorilla_flag = 1;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Check the orangutan replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( ($line_array[12] == 0 && $line_array[13] == 0) || ($line_array[12] == 0 && $line_array[14] == 0) || ($line_array[13] == 0 && $line_array[14] == 0) )
	{
		$orangutan_flag = 1;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Check the macaque replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( ($line_array[15] == 0 && $line_array[16] == 0) || ($line_array[15] == 0 && $line_array[17] == 0) || ($line_array[16] == 0 && $line_array[17] == 0) )
	{
		$macaque_flag = 1;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Calculate the total number of flags
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	$flag_counter = $human_flag + $chimp_flag + $gorilla_flag + $orangutan_flag + $macaque_flag;

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# All 5 species are flagged, which means something is wrong. Print an error message and exit
	# 	The DHS sites are created by taking the union set of 2+ peaks so at least one species must have reads in 2 replicates
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 5)
	{
		print STDERR "ERROR! No species had reads in 2 or more replicates for " . $line_array[0] . ":" . $line_array[1] . "-" . $line_array[2] . "\n";
		print STDOUT "ERROR! No species had reads in 2 or more replicates for " . $line_array[0] . ":" . $line_array[1] . "-" . $line_array[2] . "\n";
		close INPUT_FILE;
		close OUTPUT_FILE_KEPT;
		close OUTPUT_FILE_REMOVED;
	        exit(56);
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# No species are flagged so print the entry to the "kept" output file and move to the next line
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 0)
	{
		# All species have at least 2 non-zero replicates
		$entry_counter_kept++;
		print OUTPUT_FILE_KEPT $line . "\n";
		next;
	}

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# Between 1 and 4 species are flagged. Print the entry to the "removed" output file and increment the "total" counter
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# At least one of the species has zero entries for 2 or 3 replicates
	$entry_counter_removed_total++;
	print OUTPUT_FILE_REMOVED $line . "\n";

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# 1 species is flagged
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 1)
	{
		$entry_counter_removed_1_species++;

		if ($human_flag == 1)
		{
			$entry_counter_removed_H++;
			next;
		}

		if ($chimp_flag == 1)
		{
			$entry_counter_removed_C++;
			next;
		}

		if ($gorilla_flag == 1)
		{
			$entry_counter_removed_G++;
			next;
		}

		if ($orangutan_flag == 1)
		{
			$entry_counter_removed_O++;
			next;
		}

		if ($macaque_flag == 1)
		{
			$entry_counter_removed_M++;
			next;
		}
	}  # END: if ($flag_counter == 1)

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# 2 species are flagged
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 2)
	{
		$entry_counter_removed_2_species++;

		if ($human_flag == 1 && $chimp_flag == 1)
		{
			$entry_counter_removed_HC++;
			next;
		}

		if ($human_flag == 1 && $gorilla_flag == 1)
		{
			$entry_counter_removed_HG++;
			next;
		}

		if ($human_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_HO++;
			next;
		}

		if ($human_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HM++;
			next;
		}

		if ($chimp_flag == 1 && $gorilla_flag == 1)
		{
			$entry_counter_removed_CG++;
			next;
		}

		if ($chimp_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_CO++;
			next;
		}

		if ($chimp_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_CM++;
			next;
		}

		if ($gorilla_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_GO++;
			next;
		}

		if ($gorilla_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_GM++;
			next;
		}

		if ($orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_OM++;
			next;
		}

	} # END: if ($flag_counter == 2)

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# 3 species are flagged
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 3)
	{
		$entry_counter_removed_3_species++;

		if ($human_flag == 1 && $chimp_flag == 1 && $gorilla_flag == 1)
		{
			$entry_counter_removed_HCG++;
			next;
		}

		if ($human_flag == 1 && $chimp_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_HCO++;
			next;
		}

		if ($human_flag == 1 && $chimp_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HCM++;
			next;
		}

		if ($human_flag == 1 && $gorilla_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_HGO++;
			next;
		}

		if ($human_flag == 1 && $gorilla_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HGM++;
			next;
		}

		if ($human_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HOM++;
			next;
		}

		if ($chimp_flag == 1 && $gorilla_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_CGO++;
			next;
		}

		if ($chimp_flag == 1 && $gorilla_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_CGM++;
			next;
		}

		if ($chimp_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_COM++;
			next;
		}

		if ($gorilla_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_GOM++;
			next;
		}

	}  # END: if ($flag_counter == 3)

	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	# 4 species are flagged
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if ($flag_counter == 4)
	{
		$entry_counter_removed_4_species++;

		# All except macaque
		if ($human_flag == 1 && $chimp_flag == 1 && $gorilla_flag == 1 && $orangutan_flag == 1)
		{
			$entry_counter_removed_HCGO++;
			next;
		}

		# All except orangutan
		if ($human_flag == 1 && $chimp_flag == 1 && $gorilla_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HCGM++;
			next;
		}

		# All except gorilla
		if ($human_flag == 1 && $chimp_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HCOM++;
			next;
		}

		# All except chimp
		if ($human_flag == 1 && $gorilla_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_HGOM++;
			next;
		}

		# All except human
		if ($chimp_flag == 1 && $gorilla_flag == 1 && $orangutan_flag == 1 && $macaque_flag == 1)
		{
			$entry_counter_removed_CGOM++;
			next;
		}

	}  # END: if ($flag_counter == 4)

}  # END: while ($line = <INPUT_FILE>)

print STDOUT "done.\n\n";

# =================================================================================================================================================================
# Close files
# =================================================================================================================================================================
close INPUT_FILE;
close OUTPUT_FILE_KEPT;
close OUTPUT_FILE_REMOVED;

# =================================================================================================================================================================
# Calculate checksums
# =================================================================================================================================================================
$checksum1 = $entry_counter_input_file - $entry_counter_kept - $entry_counter_removed_total;
$checksum2 = $entry_counter_removed_total - $entry_counter_removed_1_species - $entry_counter_removed_2_species - $entry_counter_removed_3_species - $entry_counter_removed_4_species;
$checksum3 = $entry_counter_removed_1_species - $entry_counter_removed_H - $entry_counter_removed_C - $entry_counter_removed_G - $entry_counter_removed_O - $entry_counter_removed_M;
$checksum4 = $entry_counter_removed_2_species - $entry_counter_removed_HC - $entry_counter_removed_HG - $entry_counter_removed_HO - $entry_counter_removed_HM - $entry_counter_removed_CG - $entry_counter_removed_CO - $entry_counter_removed_CM - $entry_counter_removed_GO - $entry_counter_removed_GM - $entry_counter_removed_OM;
$checksum5 = $entry_counter_removed_3_species - $entry_counter_removed_HCG - $entry_counter_removed_HCO - $entry_counter_removed_HCM - $entry_counter_removed_HGO - $entry_counter_removed_HGM - $entry_counter_removed_HOM - $entry_counter_removed_CGO - $entry_counter_removed_CGM - $entry_counter_removed_COM - $entry_counter_removed_GOM;
$checksum6 = $entry_counter_removed_4_species - $entry_counter_removed_HCGO - $entry_counter_removed_HCGM - $entry_counter_removed_HCOM - $entry_counter_removed_HGOM - $entry_counter_removed_CGOM;

# =================================================================================================================================================================
# Print parameters to STDOUT
# =================================================================================================================================================================
print STDOUT "=================================================================================================================================================\n";
print STDOUT "Parameters\n";
print STDOUT "=================================================================================================================================================\n";
print STDOUT "Input file\t\t\t\t" . $input_filename . "\n";
print STDOUT "Output file for DHS sites removed\t" . $output_filename_removed . "\n";
print STDOUT "Output file for DHS sites kept\t\t" . $output_filename_kept . "\n\n";

# =================================================================================================================================================================
# Print checksums to STDOUT
# =================================================================================================================================================================
print STDOUT "=================================================================================================================================================\n";
print STDOUT "Checksums (All should be zero)\n";
print STDOUT "=================================================================================================================================================\n";
print STDOUT $checksum1 . "\tChecksum 1: \"total DHS sites\" - \"DHS sites kept\" - \"DHS sites removed\"\n";
print STDOUT $checksum2 . "\tChecksum 2: \"DHS sites removed\" - \"1 species flagged\" - \"2 species flagged\" - \"3 species flagged\" - \"4 species flagged\"\n";
print STDOUT $checksum3 . "\tChecksum 3: \"1 species flagged\" - counts for each species\n";
print STDOUT $checksum4 . "\tChecksum 4: \"2 species flagged\" - counts for each 2-species combination\n";
print STDOUT $checksum5 . "\tChecksum 5: \"3 species flagged\" - counts for each 3-species combination\n";
print STDOUT $checksum6 . "\tChecksum 6: \"4 species flagged\" - counts for each 4-species combination\n";

# =================================================================================================================================================================
# Print counts to STDOUT
# =================================================================================================================================================================

print STDOUT "=================================================================================================================================================\n";
print STDOUT "Counts\n";
print STDOUT "=================================================================================================================================================\n";
print STDOUT $entry_counter_input_file . "\tTotal DHS sites\n";
print STDOUT "\t" . $entry_counter_kept . "\tDHS sites kept\n";
print STDOUT "\t" . $entry_counter_removed_total . "\tDHS sites removed\n";
print STDOUT "\t\t" . $entry_counter_removed_1_species . "\t1 species flagged\n";
print STDOUT "\t\t" . $entry_counter_removed_2_species . "\t2 species flagged\n";
print STDOUT "\t\t" . $entry_counter_removed_3_species . "\t3 species flagged\n";
print STDOUT "\t\t" . $entry_counter_removed_4_species . "\t4 species flagged\n";

print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT "1 species flagged\n";
print STDOUT "\t" . $entry_counter_removed_1_species . "\tTotal\n";
print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT $entry_counter_removed_H . "\tHuman\n";
print STDOUT $entry_counter_removed_C . "\tChimp\n";
print STDOUT $entry_counter_removed_G . "\tGorilla\n";
print STDOUT $entry_counter_removed_O . "\tOrangutan\n";
print STDOUT $entry_counter_removed_M . "\tMacaque\n\n";

print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT "2 species flagged\n";
print STDOUT "\t" . $entry_counter_removed_2_species . "\tTotal\n";
print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT $entry_counter_removed_HC . "\tHuman and chimp\n";
print STDOUT $entry_counter_removed_HG . "\tHuman and gorilla\n";
print STDOUT $entry_counter_removed_HO . "\tHuman and orangutan\n";
print STDOUT $entry_counter_removed_HM . "\tHuman and macaque\n";
print STDOUT $entry_counter_removed_CG . "\tChimp and gorilla\n";
print STDOUT $entry_counter_removed_CO . "\tChimp and orangutan\n";
print STDOUT $entry_counter_removed_CM . "\tChimp and macaque\n";
print STDOUT $entry_counter_removed_GO . "\tGorilla and orangutan\n";
print STDOUT $entry_counter_removed_GM . "\tGorila and macaque\n";
print STDOUT $entry_counter_removed_OM . "\tOrangutan and macaque\n\n";

print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT "3 species flagged\n";
print STDOUT "\t" . $entry_counter_removed_3_species . "\tTotal\n";
print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT $entry_counter_removed_HCG . "\tHuman and chimp and gorilla\n";
print STDOUT $entry_counter_removed_HCO . "\tHuman and chimp and orangutan\n";
print STDOUT $entry_counter_removed_HCM . "\tHuman and chimp and macaque\n";
print STDOUT $entry_counter_removed_HGO . "\tHuman and gorilla and orangutan\n";
print STDOUT $entry_counter_removed_HGM . "\tHuman and gorilla and macaque\n";
print STDOUT $entry_counter_removed_HOM . "\tHuman and orangutan and macaque\n";
print STDOUT $entry_counter_removed_CGO . "\tChimp and gorilla and orangutan\n";
print STDOUT $entry_counter_removed_CGM . "\tChimp and gorilla and macaque\n";
print STDOUT $entry_counter_removed_COM . "\tChimp and orangutan and macaque\n";
print STDOUT $entry_counter_removed_GOM . "\tGorilla and orangutan and macaque\n\n";

print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT "4 species flagged\n";
print STDOUT "\t" . $entry_counter_removed_4_species . "\tTotal\n";
print STDOUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
print STDOUT $entry_counter_removed_HCGO . "\tAll except macaque\n";
print STDOUT $entry_counter_removed_HCGM . "\tAll except orangutan\n";
print STDOUT $entry_counter_removed_HCOM . "\tAll except gorilla\n";
print STDOUT $entry_counter_removed_HGOM . "\tAll except chimp\n";
print STDOUT $entry_counter_removed_CGOM . "\tAll except human\n\n";



