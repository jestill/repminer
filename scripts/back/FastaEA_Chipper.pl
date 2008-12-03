#!/usr/bin/perl -w
#
# MODIFIED TO SELECT SEQUENCES FROM 
#
#-----------------------------------------------------------+
#                                                           |
# FastaChipper                                              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 05/25/2006                                       |   
# UPDATED: 01/12/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
# Takes a sequence input file in fasta format and then      |
# randomly selects [s] sequences of length [l].             |
# This is an extremely slow way to do this process. Better  |
# is to load the sequences in a database with an            |
# incremental row number. Then select the sequence string   |
# from the database using the row number as the selection   |
# criteria. Alternativley, the sequence strings could be    |
# into an array that is used to fetch sequence data. This   |
# latter method will be memory intensive, but probably the  |
# fastest.                                                  |
#                                                           |
# Note that in this selection model, the total length of the|
# sequence database is not used in the selection proces.    |
# In other words all sequence records have an equal         |
# likelihood of being selected regardless of length.        | 
#                                                           |
# Example Usage:                                            |
#                                                           |
# ./RandomSeq.pl -i InFasta.fasta -o OutFile.fasta          |
#                -l LengthOfSeqToReturn                     |
#                -n Number of sequences to return           |
#-----------------------------------------------------------+

#----------------------------+
# INCLUDED MODULES           |
#----------------------------+
use Getopt::Std;              # Get options from the command line
use Bio::SeqIO;               # Bioperl Sequence interface

#-----------------------------+
# SET LOCAL VARIABLE SCOPE    |
#-----------------------------+
my $i;                         # Iterative variable for array refs
my $OutCount = "0";            # The count of sequences sent to the
                               # output file.
my @SeqData = ();              # Array to store the sequence strings
my @SeqRecord = ();            # Single row worth of data to push to
                               # the SeqData array
#my $TotalNumberSeqs = 0;     # The total number of seq in the input file
                               # This will be determined when loading
                               # the sequence strings into the array.
my $RandRow;                   # A random row number to select
my $RandSeq;                   # A random sequence string
my $RanSeqLen;                 # The length of the random sequence string
my $quiet;                     # Var to keep output quiet if desired
my $stealth;                   # Var to run program without any output
                               # to DISPLAY.
my $TotOutCount;               # Count of the total number of output sequences

#-----------------------------+
# GET COMMAND LINE VARIABLES  |
#-----------------------------+
my %Options;
getopts('i:o:n:l:m:qs', \%Options);

my $Usage = "./RandomSeq.pl -i InFasta.fasta -o OutFile.fasta\n".
    "-l LengthOfSeqToReturn\n".
    "-n Number of sequences to return\n";
# Num seq from each BAC
#my $AllByAll = $Options{i} || 
#    die "You must provide an input FASTA file\n$Usage\n";
my $InFile = $Options{i} || 
    die "You must provide a FASTA file as input:\n$Usage\n";
my $OutFile = $Options{o} ||
    die "You must provide a path for the FASTA file output\n$Usage\n";
my $MaxNum = $Options{n} ||
    die "You must provide the number of sequences to return\n$Usage\n";
#
# FOLLOWING VARS HAVE DEFAULT VALUES
#
# Default value for length is 700
my $SeqLen = $Options{l} || "700";
my $MinInputLen = $Options{m} || "700";
$stealth = $Options{s};
$quiet = $Options{q};
# If we are stealth then we are also running in quiet mode.

#----------------------------+
# FILEHANDLES I/O            |
#----------------------------+
my $seqio_object = Bio::SeqIO->new( '-file' => $InFile,
                                    '-format' => "fasta" ) ||
    die "ERROR: Can not open the sequence file:\n$InFile\n";

open (FILEOUT, ">$OutFile") ||
    die ("ERROR: Can not open the output file:\n$OutFile");

#-----------------------------+
# LOAD THE SEQUENCE STRINGS   |
# AND INFO INTO AN ARRAY      |
#-----------------------------+
# This creates a 2d array in the following form
# where N = $TotalNumOfSeqs
# 
# +-----------+-----------+
# | PrimaryId | SeqString |
# +-----------+-----------+
# |  [0,0]    |   [0,1]   |
# |    :      |     :     |
# |    :      |     :     |
# | [N-1,0]   |  [N-1,1]  |
# +-----------+-----------+
#
print "Loading sequences into the array\n";

# ARRAY INDICES START AT ZERO
# TO TAKE INTO COUNT RANGE OF
# RANDOM NUMBER GENERATOR
my $TotalNumSeqs = '0';
while (my $inseq = $seqio_object->next_seq)
{
    print "\tLoading: $TotalNumSeqs\n";
    
    # LOAD SEQUENCE DATA TO THE ARRAY
    $SeqData[$TotalNumSeqs][0] = $inseq->primary_id();
    $SeqData[$TotalNumSeqs][1] = $inseq->seq();

    $TotalNumSeqs++;
}

# If too many short sequences get loaded into
# the array. May want to add some code to
# ignore sequences that are too short to be 
# subdivided: 
# $inseq->len << $SeqLen 

#-----------------------------+
# WORKING ON SELECTION A      |
# RANDOM SEQUENCE FROM THE DB | 
#-----------------------------+
# Select a random number between 0(or 1?)
# and the total number. Using the int converts
# the random number to an integer.
my $IdiotCount = "0";
my $IdiotNum = "30000";
my $BACSeqNum = 0;

my $MaxIndexNum = $TotalNumSeqs - 1;


#for (my $SeqNum=0;$SeqNum<=$TotalNumSeqs ;$SeqNum++)
for (my $SeqNum=0;$SeqNum<=$MaxIndexNum ;$SeqNum++)
{

    print "SeqNum $SeqNum\n";


    while ($OutCount < $MaxNum)
    {
	$IdiotCount++;
	$TotOutCount++;
	
	# I will continue to use the RandRow Variable name even
	# though this is not really random any more
	$RandRow = $SeqNum;
#	$RandRow = int (rand($TotalNumSeqs));
#	print "The Random Number is $RandRow \t";
	
	$RandSeq = $SeqData[$RandRow][1];
	$RandId = $SeqData[$RandRow][0];
	$RandSeqLen = length($RandSeq);
	
	print "Ran Seq Len: $RandSeqLen\n";
	
	# Need code here to make sure that good sequnece
	# was produced as output. If true then increment
	# outcout, otherwise continue.
	my $OutProduced = &SelSubSeq($RandId, $RandSeq, $RandSeqLen);
	
	if ($OutProduced =~ "true"){$OutCount++;}
	print "\tGood Output:\t$OutProduced\n";
	
	# Will eventually need to pass additional information to the output file
	# such as the original sequence source so that I can trace these
	# back to the cluster that was expected.
	
	#-----------------------------+
	# IDIOT SWITCH TO KEEP FROM   |
	# SPINNING OFF AND INFITITE   |
	# LOOP                        |
	#-----------------------------+
	if ($IdiotCount >> $IdiotNum)
	{
	    print "HA HA. You appear to have tripped the idiot switch\n";
	    exit;
	}
	
    } # End of while $OutCount < Max Numstatement
    
    # Reset OutCount to zero
    $OutCount = 0;
} # End for SeqNum


#my $NumOut = $OutCount*
print "\n\n";
print "The FastaChipper program has completed.\n";
print "The total number of input sequences: $TotalNumSeqs\n";
print "$TotOutCount Sequences were placed in the output file:\n\t$OutFile\n";
print "\n\n";
exit;

#
sub SelSubSeq 
{
#-----------------------------+
# Given a string sequence and |
# other info this will        |
# subselect a sequence string |
# of length l. Sequences      |
# shorter then l are ignored. |                
#-----------------------------+
    my $InSeqId = $_[0];       # The primary ID of the input sequence
    my $InSeqStr = $_[1];      # The sequence string
    my $InSeqLen = $_[2];      # The length of 
    
    my $GoodOut = "";         # Was good output produced
    my $FakeCloneNum = $OutCount + 1;
    
    # Determine the Max end length
    my $MaxEndLen = $InSeqLen - $SeqLen;

    if ($MaxEndLen <= 0)
    {
	$GoodOut = "false";
    }else{
	
	my $RandStartPos = int(rand($MaxEndLen));
	my $RandSeqEndPos = $RandStartPos + $SeqLen;
	my $FakeClone = substr $InSeqStr, $RandStartPos,$SeqLen;
	#print "\tRand Start Is:".$RanNum."\n";
	# Output string includes 
	#  - Unique Name
	#  - Input sequence primary ID
	#  - Original clone length
	#  - Start position of the substring in the source sequence
	#  - End position of the substring in the source sequence
	print FILEOUT ">FakeClone".$FakeCloneNum."|".$InSeqId."|".
	    $InSeqLen."|".$RandStartPos."|".$RandSeqEndPos."\n";
	print FILEOUT $FakeClone."\n";
	$GoodOut = "true";
    }
    
    return $GoodOut;
}


#Close the output file
close FILEOUT;

# 
# TO DO
# INCLUDE INFORMATION ON THE SOURCE SEQUENCE SO THE EXPECTED
# CLUSTER CAN BE GENERATED FROM THE OUTPUT SEQUENCE FILE
# SHOULD ALSO INCLUDE THE START AND STOP POSITION FROM THE
# ORIGINAL TO SEE IF SOME REGIONS OF THE SEQUENCE ARE MORE
# LIKELY TO BE INCORRECTLY ASSIGNED.
