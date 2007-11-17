#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# FASTA ADD NUMBER                                          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/12/2006                                       |
# UPDATED: 12/18/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Add a number prefix to the id line from a fasta file     |
#  to allow for easy parsing of fasta sequences. The process|
#  also provides each sequence with a unique incremental    |
#  id that can be used for reference.                       |
#                                                           |
# USAGE:                                                    |
#  ParseASGR.pl -i infilepath -o outfilepath                |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#  -DBI                                                     |
#  -MySQL                                                   |
#                                                           |
#-----------------------------------------------------------+
#
# NOTE:
# It may just as easy to do this without any references to
# the bioperl modules. Just read through the input file and
# append the incremental number to the beginning of the 
# sequence. An advantage to using the Bioperl path is that
# a different sequence input format could be used.
#

print "The FastaAddNum program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use DBI;                       # Allows connection to MySQL database
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Text::Wrap;                # Allows printing of wrapped text
use Getopt::Std;               # Allows options flags at command line

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
# USAGE STATEMENT
my $Usage = "FastaAddNum.pl -i InFilePath -o OutFilePath\n";
#my $InSeqFile = "/home/jestill/projects/liriodendron/Lirio_data.fasta";
#my $OutSeqPath = "/home/jestill/projects/liriodendron/Lirio_data_num";
my $SeqFormat = "fasta";       # The input sequence format
my $SeqNum = 0;                # Var to keep track of nummber of seqs 
my $SeqUniqueId;               # Unique ID attributed to the seq read record
my $NewId;

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
# Declare the Options hash to hold the options from the command line
# Can provide the option for a variable input format here
# q flag is for running in quiet mode if desired.
my %Options;
getopts('i:o:q', \%Options);

my $InSeqPath = $Options{i} || 
    die "You must provide an input file\n$Usage\n";
my $OutSeqPath = $Options{o} || 
    die "You must provide an output file\n$Usage\n";
my $quiet = $Options{q};

#-----------------------------+
# FASTA FILE IO               |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$InSeqPath",
			    -format => $SeqFormat ) ||
    die "ERROR Can not open infile:\n $InSeqPath\n";

my $outAll = Bio::SeqIO->new(-file   => ">$OutSeqPath",
			     -format => "fasta") ||
    die "ERROR Can not open outfile:\n $OutSeqPath\n";

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+

while (my $seq = $inseq->next_seq) 
{
    # This process may not work with all FASTA files, some
    # of the information will be left out of the new FASTA file
    # that is produced. It would be good to just capture the entire
    # FASTA header and append the incremented number.
    $SeqNum++;
    $SeqUniqueId = $seq->primary_id;
    $NewID = $SeqNum."|".$SeqUniqueId ;


    # THE FOLLOWING MAY NOT BE NEEDED 
    # WITH THE LIRIODENDRON DATA
    #$seq->desc('');
    $seq->primary_id( $NewID );
    $seq->display_id( $NewID );

    #-----------------------------+
    # PRINT THE NEW ID            |
    #-----------------------------+
    if (! $quiet)   # If not quiet
    {
	print "\tNEW ID: ".$NewID."\n";
    }

    #-----------------------------+
    # WRITE SEQUENCE RECORDS OUT  |
    # TO A NEW FILE               |
    #-----------------------------+
    $outAll->write_seq($seq);
	

} # END OF THE FOR EVERY SEQUENCE RECORD 

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/12/2006
# - Started program
#
# 12/18/2006
# - Changed Name to FastaAddNum.pl
# - Added input variable flags
