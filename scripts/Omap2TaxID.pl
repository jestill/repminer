#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# OMAP TO TAXON ID                                          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 01/23/2007                                       |
# UPDATED: 01/23/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Takes OMAP seq data that I have modified for use         |
#  in running the RepMiner process. Parses out the          |
#  unqiue ID number and the short name for the taxon. The   |
#  output file is a *.NA file suitable for use with         |
#  Cytoscape.                                               |
#                                                           |
# USAGE:                                                    |
#  Omap2TaxID.pl -i Infile.fasta -o OutFile.NA [-q]         |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#  -DBI                                                     |
#  -MySQL                                                   |
#                                                           |
#-----------------------------------------------------------+

print "The program has started\n";
=head1 PROGRAM VARIABLES
=cut

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use Getopt::Std;               # Allows options flags at command line

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+
my $quiet;
my $SeqNum = 0;                # Var to keep track of nummber of seqs 

=head1 CMD LINE OPTIONS
=cut
#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
# Declare the Options hash to hold the options from the command line
# Can provide the option for a variable input format here
# q flag is for running in quiet mode if desired.
my %Options;
getopts('i:o:q', \%Options);
my $Usage = "Omap2TaxID.pl -i InputFilePath -o OutputFilePath";
my $SeqFile = $Options{i} || 
    die "You must provide an input file path.\n$Usage\n";
my $NAOut = $Options{o} ||
    die "You must provide an output file path for the NA file.\n$Usage\n";
my $quiet = $Options{q};

#-----------------------------+
# FILE I/O HANDLES            |
#-----------------------------+
open (SEQIN, "<".$SeqFile);
open (NAOUT, ">".$NAOut);
print NAOUT "TaxonID\n";
# WRITE IO HANDLE

=head1 MAIN BODY OF PROGRAM
The main work that the program does is here.
=cut

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+
my $TooMany = 100;

while (<SEQIN>)
{
    chomp;                     # Remove newline character
    if (m/\>.*/)
    {
	$SeqNum++;

	# BAIL OUT FOR DEBUG
	#if ($SeqNum > $TooMany){exit;}

	# The following split does work
	#my ($SeqID, $Desc) = split(/\|/);

	# The following regexp match works well
	#m/\>(.*)\|(.*)/;

	# The following regexp reduces the output
	# to just what I need
	m/\>(.*)\|(.{6})/;
	my $SeqID = $1;
	my $TaxID = $2;
	print NAOUT $SeqID."=".$TaxID."\n";
	
        # PRINT SEQ NUM INFO FOR DEBUG AND TRACKING OUTPUT
	print $SeqNum."\n";
	

    }

} # END OF FOR EVERY LINE IN THE INPUT SEQUENCE FILE

exit;
#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 01/23/2007
#  - Started the program, wrote main body of the program
#  - Played around with regexp to get just the info I needed
#    from the fasta > line.
