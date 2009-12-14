#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# tsd_check.pl - Target Site Duplication Checker            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 06/08/2006                                       |
# UPDATED: 06/08/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Extract target site duplication sequences from a genomic |
#  sequence given the start and stop positions of an LTR    |
#  retrotransposon.                                         |
#                                                           |
# INPUT VARIABLES:                                          |
#   [f] FASTA file containing ALL contigs of interest       |
#       This can be a single chrom, multiple chrom,         |
#       multiple BACs etc.                                  | 
#   [a] Annotation file in GFF Format                       |
#   [o] Output file path                                    |
#   [l] Length of LTR seq to return                         |
#       Default value is 3                                  |
#   [t] Length of putative TDS to return                    |
#       Default value is 5                                  |
#   [q] Quiet mode. Do not print program status to screen   | 
#                                                           |
# REQUIRES:                                                 |
#   - bioperl                                               |
#     Makes use of the BioPERL seqIO function               |
#                                                           |
# USAGE:                                                    |
#  tsd_sel.pl -f SeqFile.fasta -a LTRAnn.gff -o OutFile.txt |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;                # Allows for fasta sequence input/output
use Getopt::Long;              # Get input from the command line

#-----------------------------+
# GLOBAL VARIABLES            |
#-----------------------------+
# These variables need to be local to the program
# but need to be globally avaiable to the 
# subfunctions. The user will not need to 
# modify these variables.
my @StartGeneExons;            # Array to hold information for start exon
my @Annot;                     # Array to hold all annotation info for BAC
my @StartGeneList;             # Array to hold the list of gene names
                               # to use as the start genes.
my $StartPos;                  # Random start position on the BAC for
                               # one end of a simulated clone
my $SearchPos;                 # Search position to look for gene on other
                               # end of the fake clone.
my $CumSum = "0";              # Var to hold cumulative sum
my $ans;                       # Var to hold answers to questions
my $NumExonStartGene;          # Number of exons in the gene model
my $UsableLen;                 # Usable length of the exon
my $SeqId;                     # Sequence ID of contig sequence record to use
my $RandCumSum;                # Random cumulative sum
my $SeqLength;                 # Length of the sequence file
my $NumTSDMatch;               # Number of matches in TSD
my $MatchString;               # String to visualize match
my $MaxNum;                    # Max position of LTR elements when
                               # stored in the array

# Vars with default values
my $TSDLen = "5";              # Length of the target site duplication
my $LTRLen = "3";              # Length of LTR shown in verbose output
my $gff_out;                   # Output in the gff2 format
my $gff3_out;                  # Output in the gff3 format
                               # For now this assumes input in gff2 format

# Booleans
my $quiet = 0;
my $verbose = 0;
my $test= 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0;
my $show_usage = 0;
my $print_color = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "l|ltr-len=s"    => \$LTRLen,
                    "t|tsd-len=s"    => \$TSDLen,
		    "o|outfile=s"    => \$OutFile,    # Text file output
		    "gff-out=s"      => \$gff_out,    # Output in gff2 
		    "gff3-out=s"     => \$gff3_out,   # Output in gff3 format
		    "f|fasta=s"      => \$FastaFile,
		    "a|annot|gff=s"  => \$AnnotFile,
		    "q|quiet"        => \$quiet,
		    "verbose"        => \$verbose,
		    "color"          => \$print_color,
		    # ADDITIONAL INFORMATION
		    "usage"          => \$show_usage,
		    "test"           => \$test,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,);

# Load ANSICOLOR if printing in color
if ($print_color) {
    use Term::ANSIColor;           # Allows for pretty print of ANSI output
}

#-----------------------------+
# FILEHANDLES I/O             |
#-----------------------------+
if ($OutFile) {
    open (OUTFILE, ">$OutFile") ||
	die ("Can not open outfile:\n$OutFile\n");
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not print to STDOUT\n";
}

# GFF2 Format output
if ($gff_out) {
    open (GFFOUT, ">$gff_out") ||
	die ("Could not open gff outfile:\n $gff_out\n");
}

# GFF 3 Format output if desired
if ($gff3_out) {
    open (GFF3OUT, ">$gff3_out") ||
	die ("Could not open gff3 outfile:\n $gff3_out\n");
}

# Get the input sequence file in the bioperl SeqIO format
my $seqio_object = Bio::SeqIO->new( '-file' => $FastaFile,
				    '-format' => "fasta" ) ||
    die ("Can not open fasta file:\n$FastaFile");

#-----------------------------+
# LOAD ANNOTATION DATA        |
#-----------------------------+
&LoadAnnot($AnnotFile);        # Load full annotation file to 2d Array
my $LenAnnot = @Annot; 
print STDERR "$LenAnnot Annotation Objects\n" if $verbose;


# PRINT METADATA TO HEADER
print OUTFILE "# FASTA-IN:\t$FastaFile\n";
print OUTFILE "# GFF-IN:\t$AnnotFile\n";
if ($OutFile) {
    print OUTFILE "# OUT-FILE:\t$OutFile\n";
}
else {
    print OUTFILE "# OUT-FILE:\tSTDOUT\n";
}
if ($gff_out) {
    print OUTFILE "# GFF-OUT:\t $gff_out\n";
}
print OUTFILE "# LTR LENGTH:\t$LTRLen\n";
print OUTFILE "# TSD LENGTH:\t$TSDLen\n";

#-----------------------------+
# FOR EACH SEQUENCE RECORD    |
# IN THE FASTA FILE SEARCH FOR|
# LTRS TO EXTRACT             |
#-----------------------------+
while (my $inseq = $seqio_object->next_seq) {

    $SeqId = $inseq->primary_id;
    
    if ($verbose) {
	print STDERR color 'bold';
	print STDERR "SEQUENCE: ".$SeqId."\n";
	print STDERR color 'reset';
    }

    #-----------------------------+
    # FOR EACH ANNOTATED LTR      |
    # IN THE GFF FILE             |
    #-----------------------------+
    $MaxNum = $LenAnnot - 1;
    for ($i=0; $i<=$MaxNum; $i++) {
	if ($SeqId =~ $Annot[$i][0]) {
	    my $LTR_Seq = $Annot[$i][0];
	    my $LTR_ID = $Annot[$i][8];
	    my $LTR_Start = int( $Annot[$i][3] );
	    my $LTR_End = int( $Annot[$i][4] );

	    if ($verbose) {
		print STDERR "\tFeature:".$LTR_ID." : ".
		    $LTR_Start."-".$LTR_End."\n";
	    }


	    # This first pulls out TSD and portion of the LTR
	    # and then just selects the TSD for the comparison
	    my $LeftEnd = $inseq->subseq( $LTR_Start - $TSDLen, 
					  $LTR_Start + $LTRLen-1);

	    my $RightEnd = $inseq->subseq( $LTR_End - $LTRLen+1,
					   $LTR_End + $TSDLen );


	    my $LeftMatch = substr($LeftEnd, 0, $TSDLen);
	    my $RightMatch = substr($RightEnd, $LTRLen);
	    my $MatchLen = length($LeftMatch);

	    #-----------------------------+
	    # GENERATE A STRING SHOWING   |
	    # MATCHING RESIDUES           |
	    #-----------------------------+
	    # DOING THIS MATCH ON A BASE BY BASE COMPARISION
	    # SO THAT N charcaters can be avoided
	    $NumTSDMatch = "0";
	    $MatchString = "";

	    for ($j=0; $j<=$MatchLen; $j++) {

		my $L = substr($LeftMatch,$j,1);
		my $R = substr($RightMatch,$j,1);
		
		# DO NOT COUNT N BELOW
		
		if ($L =~ $R && 
		    uc($L) ne "N" && 
		    uc ($R) ne "N") {
		    $NumTSDMatch++;
		    $MatchString = $MatchString."|";  
		}
		else {
		    $MatchString = $MatchString." "; 
		}

	    }
	    

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    # IF POSITIVE FOR TSD MATCH   |
	    #-----------------------------+
	    if ($NumTSDMatch == $TSDLen) {
		
		$CumSum++;
		
		# PRINT OUTPUT IN GFF2 FORMAT
		# ASSUMES INPUT IN GFF2
		if ($gff_out) {
		    print GFFOUT $Annot[$i][0]."\t".
			$Annot[$i][1]."\t".
			$Annot[$i][2]."\t".
			$Annot[$i][3]."\t".
			$Annot[$i][4]."\t".
			$Annot[$i][5]."\t".
			$Annot[$i][6]."\t".
			$Annot[$i][7]."\t".
			$Annot[$i][8]."\n";
		}


		# PRINT OUTPUT IF GFF3 FORMAT
		if ($gff3_out) {
		    
		    # Feature ID
		    my $gff3_id = $Annot[$i][8];
		    $gff3_id = "ID=".$gff3_id."_".$CumSum;

		    # Feature Name
		    my $gff3_name = "Name=".$Annot[$i][8];

		    # Features Note, separated by pipes
		    my $gff3_note = "Note=tsd_length|".$TSDLen."|".
			"tsd_seq|".$LeftMatch;
		    
		    # Feature Attribute
		    my $gff3_attribute = $gff3_id.";".
			$gff3_name.";".
			$gff3_note;
		    
		    print GFF3OUT $Annot[$i][0]."\t".
			$Annot[$i][1]."\t".
			$Annot[$i][2]."\t".
			$Annot[$i][3]."\t".
			$Annot[$i][4]."\t".
			$Annot[$i][5]."\t".
			$Annot[$i][6]."\t".
			$Annot[$i][7]."\t".
			$gff3_attribute."\n";
		}


	    }  # End of output for TSD matches

	    print OUTFILE $Annot[$i][0]."\t".
		$Annot[$i][1]."\t".
		$Annot[$i][2]."\t".
		$Annot[$i][3]."\t".
		$Annot[$i][4]."\t".
		$Annot[$i][5]."\t".
		$Annot[$i][6]."\t".
		$Annot[$i][7]."\t".
		$Annot[$i][8]."\t".
		$NumTSDMatch."\t".
		#$LeftEnd."\t".
		#$RightEnd."\t".
		#$MatchString."\t".
		$LeftMatch."\t".
		$RightMatch."\n";
	
	    #-----------------------------+
	    # PRETTY PRINT OUTPUT TO THE  |
	    # SCREEN                      |
	    # THIS WILL BE SORTED BY THE  |
	    # BAC OCCURRENCE IN THE FASTA |
	    # INPUT FILE                  |
	    #-----------------------------+
	    if ($verbose) {
		if ($NumTSDMatch == $TSDLen) {
		    print STDERR color 'green';
		    print STDERR "\t\tMATCH:".$NumTSDMatch."\n";
		    print STDERR "\t\t";
		    print STDERR " " x $LTRLen;
		    print STDERR $LeftEnd."\n";
		    print STDERR "\t\t";
		    print STDERR " " x $LTRLen;
		    print STDERR $MatchString."\n";
		    print STDERR "\t\t".$RightEnd."\n";
		    print STDERR color 'reset';
		}
		else {
		    print STDERR color 'red';
		    print STDERR "\t\tMATCH:".$NumTSDMatch."\n";
		    print STDERR "\t\t";
		    print STDERR " " x $LTRLen;
		    print STDERR $LeftEnd."\n";
		    print STDERR "\t\t";
		    print STDERR " " x $LTRLen;
		    print STDERR $MatchString."\n";
		    print STDERR "\t\t".$RightEnd."\n";
		    print STDERR color 'reset';
		}
	    }
	    
	    #-----------------------------+
	    # SAVE RESULTS TO ANNOTATION  |
	    # ARRAY                       |
	    #-----------------------------+
	    $Annot[$i][9] = $NumTSDMatch;
	    $Annot[$i][10] = $LeftEnd;
	    $Annot[$i][11] = $RightEnd;
	    $Annot[$i][12] = $MatchString;
	    $Annot[$i][13] = $LeftMatch;
	    $Annot[$i][14] = $RightMatch;



	} # End of found seq/annotation match

    } # End of for each annotation in the gff file

} # End of for each sequence in the output file


if ($verbose) {

    print STDERR "\n\n";
    
    print STDERR color 'bold';    
    print STDERR color 'black';
    print STDERR "CONTIG\tELEMENT\tNUM_MATCH\tLEFT\tRIGHT\t\n";
    print STDERR color 'reset';
    
    for ($i=0; $i<=$MaxNum; $i++) {
	print STDERR $Annot[$i][0]."\t";
	print STDERR $Annot[$i][8]."\t";
	print STDERR $Annot[$i][9]."\t";
	print STDERR $Annot[$i][10]."\t";
	print STDERR $Annot[$i][11]."\t";
	print STDERR $Annot[$i][13]."\t";
	print STDERR $Annot[$i][14]."\t";
	print STDERR "\n";
    }

}

#-----------------------------+
# Close output files
#-----------------------------+

if ($gff3_out) {
    close (GFF3OUT);
}


if ($gff_out){
    close (GFFOUT);
}

close (OUTFILE);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

#-----------------------------+
# LOAD BAC ANNOTATION DATA    |
#-----------------------------+
# This should be in tab delim gff format.
sub LoadAnnot  {
    my $GffFile = $_[0];
    my (@line, $SeqName, $Source, $Feature, $Start, $End);
    my ($Score, $Strand, $Frame);

    open (INFILE, $GffFile) || 
	die ("Can not load file $GffFile \n");
    while (<INFILE>) {

	# Get in from tab delim text file 
	my @line = split(/\t/, $_); 
	my $Seqname = $line[0];
	my $Source = $line[1];
	my $Feature = $line[2];
	my $Start = $line[3];
	my $End = $line[4];
	my $Score = $line[5];
	my $Strand = $line[6];
	my $Frame = $line[7];
	my $Attribute = $line[8]; # Gene name column
	chomp $Attribute;         # Get rid of newline character

	# Can get rid of comment below to show the loading
	# of the individual features
	#print "Load: $Attribute $Start\-$End\n";

	# Load the information into the @Annot array
	push (@Annot, [$Seqname, $Source, $Feature, $Start, $End,
		       $Score, $Strand, $Frame, $Attribute]);

    }
    close (INFILE);

    my $TestAryLen = @Annot;

    if ($TestAryLen < 1) {
	print color 'red';
	print "ERROR: I did not find any annotation data\n";
	print color 'reset';
	exit;
    }

}

#-----------------------------+
# USER VARIABLE CHECK         |
#-----------------------------+
# Let's the user check the input variables

sub UserVarCheck {
    print STDERR "\nYOU HAVE SELECTED THE FOLLOWING VARIBLES:\n";
    #-----------------------------+
    # FASTA FILE                  |
    #-----------------------------+
    print STDERR "FASTA File:\n";
    if (-e $FastaFile) {
	print STDERR "\t$FastaFile\n";
    }
    else {
	print STDERR "\tWARNING: $FastaFile \n\tDoes Not Exist\n";
    }
    
    #-----------------------------+
    # ANNOTATION FILE             |
    #-----------------------------+
    print "ANNOTATION FILE:\n";
    if (-e $AnnotFile) {
	print STDERR "\t$AnnotFile\n";
    }
    else {
	print STDERR "\tWARNING: $AnnotFile \n\tDoes Not Exist\n";
    }
    
    #-----------------------------+
    # OUTPUT FILE                 |
    #-----------------------------+
    print "OUTPUT FILE:\n";
    if (-e $OutFile) {
	print "\t$OutFile already exists\n";
	print "\tThe existing file will be overwritten\n";
    }
    else {
	print "\t$OutFile\n";
    }
    
    #-----------------------------+
    # LTR LENGTH                  |
    #-----------------------------+
    print "LTR LENGTH:\n\t$LTRLen\n";

    #-----------------------------+
    # TSD LENGTH                  |
    #-----------------------------+
    print "TSD LENGTH\n\t$TSDLen\n";
    
    #-----------------------------+
    # QUIET MODE                  |
    #-----------------------------+
    # If this is in quiet mode the answer to
    # the question is Y and the user will not
    # need to supply feedback at the command line.
    if (! $quiet){
	$ans = &UserFeedback("\nDo you want to continue (y/n)?");
    }else{
	$ans = "Y";
    }
    
    if ($ans =~"N" )
    { 
	print "Goodbye\n";
	exit;
    }elsif ($ans =~"Y"){
	print "Starting the process...\n";	
    }    
}



#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/08/2006
# - Begin modification of DblHit to do a dts_sel for clem
#
# 05/15/2009
# - Modifiying to long intput format
# - Removing user feedback to facilitate use on r cluster
# - Added output of fetures with tsds to a separate file
#
# 05/18/2009
# - Added count of cumulative sum
# - Added option to provide output in GFF format
#     - This assumes input is GFF2 format
#     - This will append cum_sum to input ID to insure unique id
#     - This will give information on the TSD in the Note field 
#       of the GFF3 formatted attributes file. This note 
#       sep by pipe characters
# 
#-----------------------------------------------------------+
# MODEL                                                     |
#-----------------------------------------------------------+
#
# Given a genomic sequence and Start/End Positions of the
# LTR Retrotransposon for a given genomic sequence:
#
# GENOME
# |=========================================================================|
#                    
#                  START                           END
#                   ||                             ||
# LTR               \/                             \/
#                    [=5'LTR=]-------------[3'=LTR=]
#                                      
# RETURNED          
# SEQS              =:                              :=
#                  5'Seq                             3'Seq
#
# GIVEN: 
# START is position of start of 5' LTR
# END is position of end of 3' LTR
