#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# CALCULATE PERCENT IDENTITY                                |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 03/04/2007                                       |
# UPDATED: 03/04/2007                                       |
# DESCRIPTION:                                              |
#  Given a fasta file that represent the an aligned set of  |
#  sequences. Determine the percent identity using a simple |
#  count algorithm.                                         |
#                                                           |
# NOTES:                                                    |
#  - Array indices start at zero                            |
#-----------------------------------------------------------+

# TODO
# Add a threshold for reporting
print "The CalcPID program has started.\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Getopt::Std;               # Allows options flags at command line


#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
#my $SeqFile;                   # Input fasta file
#my $OutFile;                   # Output pairise PID file
my @SeqStrings;                # Sequence Strings
my @SeqIds;                    # Sequence IDS
my $TotComp;                   # The total number of pairwise comparisions

=head1 CMD LINE OPTIONS
=cut
#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
# q flag is for running in quiet mode if desired.
my %Options;
getopts('i:o:t:q', \%Options);
my $Usage = "CalPid.pl -i InFilePath.fasta -o OutFilePath.txt".
    " (-t threshold)\n";
my $SeqFile = $Options{i} || 
    die "You must provide an input file path\n$Usage\n";
my $OutFile = $Options{o} || 
    die "You must provide an output file path\n$Usage\n";
my $Thresh = $Options{t} || 0.75;

=head1 FILE I/O
Open files for input
=cut
#-----------------------------+
# FILE/DB IO                  |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$SeqFile",
			    -format => 'fasta' ) ||
    die "Can not open input file:\n$SeqFile\n";

# The SIF output file for cytoscape
open (SIFOUT, ">$OutFile") ||
    die "Can not open output file:\n$OutFile\n";

=head1 MAIN BODY
Main body of the program.
=cut
#-----------------------------+
# LOAD ALL SEQUENCE STRINGS   |
# TO AN ARRAY SEQ ID IN SECOND|
# ID ARRAY                    |
#-----------------------------+
print "Loading arrays.\n"; 

while (my $seq = $inseq->next_seq) 
{
    # SEQUENCE STRING
    push (@SeqStrings, $seq->seq);
    

    # SEQUENCE 
    push (@SeqIds, $seq->primary_id());

    print "PUSHING: ".$seq->primary_id()."\n\t".$seq->seq."\n";

};

my $NumSeqs = @SeqIds;
my $MaxNum = $NumSeqs - 1;
print $NumSeqs." will be compared\n";

#-----------------------------+
# COMPARE SEQ STRINGS         |
#-----------------------------+
# The name used to refer to the sequence IDS will 
# start at one even though the arrays are at zero.
#
# I will just do this for i<j for now
print "Determing all pairwise PIDS.\n";

for ( my $i=0 ; $i<=$MaxNum; $i++)
{
    # DEBUG PRINT FOLLOWS
    #print $SeqIds[$i]."\n";
    for ( my $j=0 ; $j<=$MaxNum; $j++)
    {
	if ($i < $j)
	{
	    # DEBUG PRINT FOLLOWS
	    #print $i."\t".$j."\n";
	    
	    my $MatchCount = 0;                   # Inititalize the MatchCount
	    my $PID = 0;
	    my @SeqI = split(//,$SeqStrings[$i]); # Load Seq[i] to array
	    my @SeqJ = split(//,$SeqStrings[$j]); # Load Seq[j] to array
	    my $SeqILen = @SeqI;                  # Length of the SeqIArray
	    # The following for debug
	    #print $SeqILen."\n";
	    #exit;
	    my $SeqMaxNum = $SeqILen - 1;
	    #$SeqILen = $SeqILen - 1;              #
	    
	    # $x is position along the two sequences

	    #-----------------------------+
	    # COUNT MATCHING BASES        |
	    #-----------------------------+
	    for (my $x=0; $x<$SeqMaxNum; $x++)
	    {
		# If the two residues match increment $MatchCount
		if ($SeqI[$x] =~ $SeqJ[$x])
		{
		    $MatchCount++;
		}else{
		    
		}
	    }
	    
	    #-----------------------------+
	    # CALCULATE PERCENT IDENTITY  |
	    #-----------------------------+
	    # Need to determine if $Count is greter then zero
	    # to avoid any divide by zero problems. Report
	    # Null for divide by zero.
	    if ($SeqILen > 0)
	    { 
		$PID = $MatchCount/$SeqILen;
	    } else {
		$PID = "Null";
	    }

	    #-----------------------------+
	    # PRINT RESULT TO SCREEN      |
	    #-----------------------------+
#	    print $SeqIds[$i]."\t".$SeqIds[$j]."\n";
#	    print $SeqIds[$i]."\t".$MatchCount."\t".$SeqIds[$j]."\n";

	    # ONLY PRINT FOR PID GREATER THEN THRESHOLD VALUE
	    if ($PID > $Thresh)
	    {

		print $SeqIds[$i]."\t".$PID."\t".$SeqIds[$j]."\n";

		# XCrd and YCrd are numeric IDs of the
		# sequences. Thess must be incremented
		# by one since IDs start at 1 and arrays
		# are all indexed starting at zero
		my $XCrd = $i + 1;
		my $YCrd = $j + 1;
		print SIFOUT $XCrd."\tbl\t".$YCrd."\n";

	    } # END of IF PID > THRESH

	    #exit; # Temp exit

	} # END OF IF i < j
	
    } # END OF FOR J

} # END OF FOR I

#-----------------------------+
# CLOSE THE OUTPUT FILE       |
# HANDLE AND EXIT PROGRAM     | 
#-----------------------------+
print "Closing output file handle\n";
close SIFOUT;
exit;

# The following just for debug right now
#print "The number of sequences is ".$NumSeqs."\n";
#print "The max array reference is ".$MaxNum."\n";

#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 03/04/2007
# - Program started
