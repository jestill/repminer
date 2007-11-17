#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# Fasta2NA : FASTA to Node Attribute File                   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 03/18/2007                                       |
# UPDATED: 03/18/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a RepMiner formatted FASTA file converts the       |
#  sequence records into a Noded Attribute file for use in  |
#  Cytoscape.                                               |
#                                                           |
# USAGE:                                                    |
#  Fasta2Na.pl -i InFastaFile.fasta -o OutNodeAtt.NA        |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#                                                           |
#-----------------------------------------------------------+

print "Fasta2NA has started\n.";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Getopt::Std;               # Allows options flags at command line

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
my $Usage = "Fasta2Na.pl -i InputFilePath -o OutputFilePath\n";

my %Options;
getopts('i:o:q', \%Options);


my $SeqFile = $Options{i} || 
    die "You must provide an input file path\n$Usage\n";

my $OutPath = $Options{o} || 
    die "You must provide an output path\n$OutPath\n$Usage\n";
# VARIABLES NOT REQUIRED AT THE COMMAND LINE
my $quiet = $Options{q};

#-----------------------------+
# FILE IO                     |
#-----------------------------+

my $inseq = Bio::SeqIO->new(-file   => "<$SeqFile",
			    -format => 'fasta' ) ||
    die "Can not open the intput file:\n$inseq\n";

open (OUT, ">$OutPath") ||
    die "Can not open output file:\n$OutPath\n";
print OUT "FullSeqDat\n";

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+
while (my $seq = $inseq->next_seq) 
{
    # Print to screen to show process

    # Print to output file

    my @IdSplit = split (/\|/ ,$seq->primary_id());

#    print "Processing: ".$seq->primary_id()."\n";
    print "Processing: ".$IdSplit[0]."\n";
    print OUT $IdSplit[0];
    print OUT "=";
    print OUT $seq->seq();
    print OUT "\n";
} # End of while $seq input loop

    
#-----------------------------+
# CLOSE OUTPUT FILE AND       |
# EXIT PROGRAM                |
#-----------------------------+
close OUT;
exit;

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 03/18/2007
# - Wrote body of the program with input and output flags
