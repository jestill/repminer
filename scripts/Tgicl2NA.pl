#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# TGICL CLUSTERS TO NODE ATTRIBUTES                         |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 01/23/2007                                       |
# UPDATED: 01/23/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  The purpose of this program is to convert the *_cluster  |
#  output from the tgicl tclust program to a format that    |
#  can be drawn on top of the networks that are being       |
#  visualized in cytoscape.                                 |
#                                                           |
# USAGE:                                                    |
#  Tgicl2NA.pl -i Infile.fasta -o OutFile.NA [-q]           |
#                                                           |
# REQUIREMENTS:                                             |
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
my $ClusNum = 0;                # Var to keep track of nummber of seqs 
my $NumSeqs;                    # The number of sequences in the cluster
my $ClustId;                    # Cluster ID
my $SeqId;                      # Unique sequence ID for working with
                                # output in cytoscape.
my $SeqDesc;                    # The full sequence description
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
my $ClustFile = $Options{i} || 
    die "You must provide an input file path.\n$Usage\n";
my $NAOut = $Options{o} ||
    die "You must provide an output file path for the NA file.\n$Usage\n";
my $quiet = $Options{q};

#-----------------------------+
# FILE I/O HANDLES            |
#-----------------------------+
open (CLUSIN, "<".$ClustFile);
open (NAOUT, ">".$NAOut);
print NAOUT "ClustID\n";
# WRITE IO HANDLE

=head1 MAIN BODY OF PROGRAM
The main work that the program does is here.
=cut

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+
my $TooMany = 100;

while (<CLUSIN>)
{
    chomp;                     # Remove newline character


    # Get cluster ID
    if (m/\>.*/)
    {
	$ClusNum++;
	($ClustId, $NumClust) = split;

	print $ClustId."\n";
    }else{
	my @SeqList = split;
	
	# Need to iterate across the array
	foreach (@SeqList)
	{
	    # The $_ string will hold the individual element in SeqList
	    print "\t".$_."\n";
	    # since I am splitting the defaul $_ no need to give
	    # string value as a variable
	    my @SeqParse = split(/\|/);
	    # The unique sequence id is the first element in the array
	    # that is split by the pipe character
	    print "\t".$SeqParse[0]."\n";
	    print NAOUT $SeqParse[0]."=".$ClustId."\n";
	}

	#print "\t".$SeqList[0]."\n";

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
