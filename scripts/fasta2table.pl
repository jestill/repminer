#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta2table.pl - Convert fatsta file to tab delim text    |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 03/11/2008                                       |
# UPDATED: 03/11/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Add a number prefix to the id line from a fasta file     |
#  to allow for easy parsing of fasta sequences.            |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Getopt::Long;               # Allows options flags at command line

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $infile;                    # Path to the infile
my $outfile;                   # Path to the outfile
my $feature_name;              # Name of the sequence feature

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "f|feat=s"    => \$feature_name,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);



my $inseq = Bio::SeqIO->new(-file   => "<$infile",
			    -format => 'fasta' );

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+

while (my $seq = $inseq->next_seq) 
{

#    print STDOUT $seq->primary_id()."\n";

    my ($node_id,$cruft) = split (/\|/,$seq->primary_id());
    print STDOUT "$node_id\t".$seq->seq()."\t$feature_name\t\n";
    
} # END OF THE FOR EVERY SEQUENCE RECORD 

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 03/11/2008
# - Started basic program
