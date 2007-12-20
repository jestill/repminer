#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# RepMiner : PAN ANALYSIS 
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/17/2006                                       |
# UPDATED: 08/03/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the total pipeline of analyses for PANS that were    |
#  selected from Cytoscape. This will first transform the   | 
#  output text file from Cytoscape to a fasta file and then |
#  cluster the PANs with TGICL. The TGICL clusters are next |
#  split into a separate dir for each cluster. The clusters |
#  are identified by tBLASTx to a series of repeat databases|
#  and are then subjected to HMMER analysis to identify the | 
#  potential MITES and MULES.          
#                                                           |
# USAGE:                                                    |
#  BLParser.pl                                              |
#                                                           |
# DEPENDENCIES:                                             |
#  -BioPerl                                                 |
#                                                           |
#-----------------------------------------------------------+

=head1 INCLUDED CONTENT
PERL Modules required by the program
=cut
#-----------------------------+ 
# INCLUDES                    |
#-----------------------------+
use Text::Wrap;                # Text wrap for output of seq strings
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Bio::SearchIO;             # Parse BLAST output
use DBI();                     # Connect to the database
use Term::ANSIColor;           # Allows to print error messages in red
use Getopt::Std;               # Get options from the command line
use File::Copy;                # To copy files
use Bio::Tools::HMMER::Results;# To work with the HMMER output

#-----------------------------+
# GET START TIME              |
# AND SCRIPT WORKING DIR      |
#-----------------------------+
my $StartTime = time;

=head1 VARIABLES
The program variables
=cut
#-----------------------------------------------------------+
# ASSEMBLY RELATED OPTIONS                                  |
#-----------------------------------------------------------+

#-----------------------------+
# REPEAT BLAST VARIABLES      |
#-----------------------------+
my $MaxE = '0.0001';          # Max E value for repblast
my $MinQryLen = '50';         # Minimum query length
my $MinScore = '50';          # Minimum BitScore

#-----------------------------+
# TGICL CLUSTERING OPTIONS    |
#-----------------------------+
# Set options to
#  - 75% identity
#  - 21 base overlap (The min allowable with tgicl)
#my $Topt = "-O '-p 75 -o 21'";   # 21 Pans -> 44 Contigs
#my $Topt = "-O '-p 80 -o 30'";   # 21 Pans -> 44 Contigs
my $Topt = "-O '-p 90 -o 30'";    # 21 Pans -> 45 Contigs


#-----------------------------------------------------------+
# TAXON SPECIFIC OPTIONS                                    |
#-----------------------------------------------------------+
my $Species = "cen";
#my $Species = "pen";

# SET SCOPE OF VARS
my $BaseDir;
my $SourceDb;
my $SourceDbName;
my @Name;

##-----------------------------+
## CENCHRUS VARIABLES          |
##-----------------------------+
## The species that was the source
#my $Species = "cen";               # Species to look up id number for
#

if ($Species =~ "cen")
{
## THE BASE DIR THAT ALL PAN FILES EXIST IN
$BaseDir = "/home/jestill/projects/asgr/Asgr_C_Net20060802/";
#
## BLAST FORMATTED SOURCE FILE DIR AND DB NAME
$SourceDb = "/home/jestill/projects/asgr/blast/asgr_c/asgr_c";
$SourceDbName = "asgr_c";
#
## NAME OF PAN FILES
## Expects data in the format:
##  /BaseDir/Name/Name.txt"
@Name = ( "CPAN001",
	  "CPAN002",
	  "CPAN003",
	  "CPAN004",
	  "CPAN005",
	  "CPAN006",
	  "CPAN007",
	  "CPAN008",
	  "CPAN009",
	  "CPAN010",
	  "CPAN011",
	  "CPAN012",
	  "CPAN013",
	  "CPAN014",
	  "CPAN015",
	  "CPAN016",
	  "CPAN017",
	  "CPAN018",
	  "CPAN019",
	  "CPAN020",
	  "CPAN021",
	  "CPAN022",
	  "CPAN023",
	  "CPAN024",
	  "CPAN025",
	  "CPAN026",
	  "CPAN027",
	  "CPAN028",
	  "CPAN029",
	  "CPAN030",
	  "CPAN031",
	  "CPAN032",
	  "CPAN033",
	  "CPAN034",
	  "CPAN035",
	  "CPAN036",
	  "CPAN037",
	  "CPAN038",
	  "CPAN039",
	  "CPAN040",
	  "CPAN041",
	  "CPAN042",
	  "CPAN043",
	  "CPAN044",
	  "CPAN045",
	  "CPAN046",
	  "CPAN047",
	  "CPAN048"
	  );
}elsif ($Species =~ "pen")
{
    
#-----------------------------+
# PEN VARIABLES               |
#-----------------------------+
    $BaseDir = "/home/jestill/projects/asgr/Asgr_p_Net/";
    $SourceDb = "/home/jestill/projects/asgr/blast/asgr_p/asgr_p";
    $SourceDbName = "asgr_p";

    @Name = ( "PPAN001",
	      "PPAN002",
	      "PPAN003",
	      "PPAN004",
	      "PPAN005",
	      "PPAN006",
	      "PPAN007",
	      "PPAN008",
	      "PPAN009",
	      "PPAN010",
	      "PPAN011",
	      "PPAN012",
	      "PPAN013",
	      "PPAN014",
	      "PPAN015",
	      "PPAN016",
	      "PPAN017",
	      "PPAN018",
	      "PPAN019",
	      "PPAN020",
	      "PPAN021",
	      "PPAN022",
	      "PPAN023",
	      "PPAN024",
	      "PPAN025",
	      "PPAN026",
	      "PPAN027",
	      "PPAN028",
	      "PPAN029",
	      "PPAN030",
	      "PPAN031",
	      "PPAN032",
	      "PPAN033",
	      "PPAN034",
	      "PPAN035",
	      "PPAN036",
	      "PPAN037",
	      "PPAN038",
	      "PPAN039",
	      "PPAN040",
	      "PPAN041",
	      "PPAN042",
	      "PPAN043",
	      "PPAN044",
	      "PPAN045",
	      "PPAN046",
	      "PPAN047",
	      "PPAN048",
	      "PPAN049",
	      "PPAN050",
	      "PPAN051",
	      "PPAN052",
	      "PPAN053",
	      "PPAN054",
	      "PPAN055",
	      "PPAN056",
	      "PPAN057",
	      "PPAN058",
	      "PPAN059",
	      "PPAN060"
	      );
}else{
    print "\nI DO NOT RECOGNIZE THE SPECIES\n";
    exit;
}


=head1 TEMP 
This is an attempt to read in files from a folder
to use as the qry dataset.
=cut
#-----------------------------------------------------------+
# THE FOLLOWING IS AN AT WORK IDEA
#
#my $PanPath = $BaseDir."PAN/"; # Source of all PANS from Cytoscape
#
## Read in the files from the base PanPath
#opendir( DIR, $PanPath );
#  my @fileList = grep !/^\.\.?$/, readdir DIR ;     # Read the directory ignoring the . and ..
#closedir( DIR );                                    # Close the directory
# 
##
## make a directory for the output file unless already exists
#mkdir $output, 0777 unless (-e $output); # set permissions
## DO THIS FOR EACH FILE IN THE DIR
#
#
# END OF AT WORK IDEA
#-----------------------------------------------------------+



#-----------------------------------------------------------+
# GLOBALS                                                   |
#-----------------------------------------------------------+
#-----------------------------+
# OUTPUT FILE PATHS           |
#-----------------------------+
my $MetaPath = $BaseDir."MetaData.html";
my $ConIDPath = $BaseDir."ContigId.NA";
my $AllFastaPath = $BaseDir."AllPAN.fasta";

# THE BASE DIR THAT THE SCRIPT IS WORKING IN
my $ScriptDir = "/home/jestill/projects/asgr/scripts/";
# The padded length to use 
# ie (3 results in CONTIG001,CONTIG002 etc.
my $pad_len = 3;

#-----------------------------+
# GET AND PRINT INFORMATION   |
# ON PAN NAME LIST            |
#-----------------------------+
my $LenName = @Name;
#show the number of PANs to analyze
print color 'blue bold';
print $LenName." PANs to analyze.\n";
print color 'reset';

#-----------------------------+
# INITIALIZE VARIABLES AND    |
# SET SCOPE                   |
#-----------------------------+
my $KillCount = '0';
my $TotConCount = '0';         #  Total number of contigs

#-----------------------------+
# VARS FOR META DATA TABLE    |
#-----------------------------+
# These vars need to be global in scope so 
# that the subfunctions can write to vars
my $AsmLen;                    # Length of the assembly
my $NumAsmSeqs;                # Number of seqs in the assembly
my $NumAsm;                    # Number of assemblies in a PAN

=head1 COMMAND LINE
Get variables from the command line.
=cut

#-----------------------------+
# GET VARS FROM COMMAND LINE  |
#-----------------------------+
my %Options;
getopts('qab:', \%Options);
my $quiet = $Options{q};       # Run in quiet mode
my $all = $Options{a};         # Show all hits from the repeat blast
my $desc = $Options{d};        # Descriptive blast information
my $BlProg = $Options{b} ||    # Type of blast program to run 
    "tblastx";                 # default is tblastx 

#-----------------------------+
# CHECK THAT A VALID BLAST    |
# PROGRAM WAS SELECTED        |
#-----------------------------+
unless ($BlProg =~ "tblastx" || $BlProg =~ "blastn")
{
    print("\007");
    print color 'red';
    print "A valid BLAST program was not selected.\n";
    print "Valid values include:".
	"\n\t*blastn\n\t*tblastx\n";
    print color 'reset';
    exit;
}

=head1 DATABASE CONNECT
Connect to the datbase. This is needed
to get the BAC id for the ASGR seqs.
=cut
#-----------------------------------------------------------+
# DATABASE CONNECTIONS                                      |
#-----------------------------------------------------------+

#-----------------------------+
# DB VARIABLES                |
#-----------------------------+
my $DbName = "dbASGR";             # Database name to connect to
my $tblSrcTable = "tblSeqData";    # Table with src sequence data
my $tblDatCat = "tblDatCat";       # Table with information about data
                                   # categories. Can include color
                                   # and symbol for drawing.
my $DbUserName = "jestill";        # User name for the db connection

#-----------------------------+
# GET DB USER PASSWORD        |
#-----------------------------+
print "\nPassword for $DbUserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
$DbUserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $DbUserPassword;

#-----------------------------+
# FILE I/O DB CONNECT         |
#-----------------------------+
my $dbh = DBI->connect("DBI:mysql:database=$DbName;host=localhost",
		       $DbUserName, $DbUserPassword,
		       {'RaiseError' => 1});


$RepDbName = "dbRep";
my $RepDB = DBI->connect("DBI:mysql:database=$RepDbName;host=localhost",
			   $DbUserName, $DbUserPassword,
			   {'RaiseError' => 1});

=head1 OPEN FILES
Open files for Input/Output
=cut
#-----------------------------------------------------------+
# OPEN FILES FOR I/O                                        |
#-----------------------------------------------------------+
#-----------------------------+
# ALL PANCON IN ONE FASTA FILE|
#-----------------------------+
#open (ALLFAST, ">".$AllFastaPath);
my $all_seq_out = Bio::SeqIO->new (
				   '-format' => 'fasta',
				   '-file' => ">".$AllFastaPath ) ||
    die "Count not open $AllFastaPath\n";

#-----------------------------+
# META DATA HTML              |
#-----------------------------+
# HTML formatted meta data that connects
# to the files that are created and gives
# an overview of the analysis
open (METAOUT, ">".$MetaPath);
print METAOUT "<HTML>\n";
print METAOUT "<HEAD>\n";
print METAOUT "<TITLE>JABA BLAST METADATA</TITLE>".
    "</HEAD>\n<BODY>\n";
print METAOUT "<H1>".$Species." PAN Analysis</H1>\n";
print METAOUT "<HR>\n";
print METAOUT "<P><B>Number of PANS: </B>".$LenName."</P>";
print METAOUT "<P><B>TGICL Options: </B>".$Topt."</P>\n";
print METAOUT "<P><B>REPEAT ID BLAST VARIABLES:</B></P>\n";
print METAOUT "<UL>".
    "<LI><B>PROGRAM: </B>".$BlProg."</LI>".
    "<LI><B>MAX-E: </B>".$MaxE."</LI>".
    "<LI><B>MIN-SCORE: </B>".$MinScore."</LI>".
    "<LI><B>MIN-QRY LEN: </B>".$MinQryLen."</LI>".
    "</UL>";
print METAOUT "<HR>\n";

#-----------------------------+
# CLUSTER DATA NETWORK ATTR   |
#-----------------------------+
# The cluster data network attribute file
# This indicates the contig ID that the 
# cluster belongs to.
open ( CONOUT ,">".$ConIDPath);
print CONOUT "ContigID\n";

=head1 MAIN BODY
Main Body of the program
=cut

for ( $i=0; $i<$LenName; $i++)
{


    # ONLY DO FIRST TWO RECORDS FOR DEBUG/TEST
    #if ($i == 3){exit;}

    # Had to add a Kill switch because
    # the loop would spiral out of control
    #$KillCount++;
    #if ($KillCount > 100){exit;}

    #-----------------------------+
    # Show the analysis process   | 
    # number and cpan id          |
    #-----------------------------+
    my $ProcNum = $i + 1;          # Process Number
    my $ProcPan = $Name[$i];       # PAN that is being processed
    print color 'bold blue';
    print "\n========================================\n";
    print " PROCESS $ProcNum of $LenName\n";
    print " PAN $ProcPan\n";
    print "========================================\n\n";
    print color 'reset';

    print METAOUT "<H3>$ProcPan</H3>\n";

    my $In = $BaseDir.$Name[$i]."/".$Name[$i].".txt";
    my $Out = $BaseDir.$Name[$i]."/".$Name[$i].".fasta";
    #print "\t".$In."\n\t->".$Out."\n";

    my $Img = $BaseDir.$Name[$i]."/".$Name[$i].".png";
    #print "\t".$In."\n\t->".$Out."\n";

    #-----------------------------+
    # PRINT IMAGE OF THE PAN IF   |
    # ON EXISTS                   |
    #-----------------------------+
    if (-e $Img)
    {
	# Percent did not render well in Mozilla Firefox so I just
	# use a straight width with and make browser keep aspect
	# Link to the original image is provided
	#print METAOUT "<P><IMG BORDER=2 WIDTH=\"100\" SRC=".
	#    $Img."></P>\n";
	print METAOUT "<P><A HREF=".$Img."><IMG BORDER=2 WIDTH=\"100\" SRC=".
	    $Img."></A></P>\n";
    }

    # INCLUDE AN IF THEN STATEMENT TO FIRST CHECK
    # IF THE FILE ACTUALLY EXIST BEFORE TRYING 
    # TO DO THE TRANSFORMATION
    if (-e $In)
    {
	$PanWorkDir =  $BaseDir.$Name[$i]."/";
	# Change the working dir to where I would like
	# the tgicl information stored
	chdir ( $PanWorkDir );

	#-----------------------------+
	# CONVERT CYTOSCAPE TO FASTA  |
	# FORMAT                      |
	#-----------------------------+
	print color 'bold blue';
        print "\n$Name[$i]: Converting to FASTA.\n";
	print color 'reset';
	&Cyto2Fasta( $In, $Out);
	
	#-----------------------------+
	# RUN TGICL                   |
	#-----------------------------+
	print color 'bold blue';
	print "\n$Name[$i]: Clustering with tgicl\n";
	print color 'reset';
	&RunTGICL ( $Out );

	#-----------------------------+
	# FOR EACH FASTA CONTIG IN THE|
	# asm_1/contig FILE CREATE A  |
	# SUBDIR IN THE PARENT DIR TO |
	# HOLD THE ANALYSIS AND RUN   |
	# THE ASSEMBLY ANALYSIS       |
	#-----------------------------+
	print color 'bold blue';
	print "\n$Name[$i]: Analyzing TGICL Assemblies\n";
	print color 'reset';
	&AsmAnalysis( $PanWorkDir , $Name[$i] );

    }else{
	# Print error message out to terminal
	print("\007");             # Sound the bell for file not exist errors
	print color 'bold red';
	print "The file does not exist.\n$In\n";
	print color 'reset';
	
	# Print error message out to metadata file
	print METAOUT "<P><FONT COLOR=#CD0000>File could not be found:<BR>\n";
	print METAOUT $In."</FONT></P>\n";
    }
}

=head1 CLOSE FILES
Close open files and write footer information.
=cut
#-----------------------------+
# GET END TIME AND TOTAL TIME |
#-----------------------------+
my $EndTime = time;
my $TotTime = $EndTime - $StartTime;


#-----------------------------+
# CLOSE METADATA HTML FILE    |
#-----------------------------+
print METAOUT "<P>\n";
print METAOUT "<HR>\n";
print METAOUT "<P><P>";
print METAOUT "<H3>Summary</H3>\n";
print METAOUT "<P>Total Contigs: ".$TotConCount."</P>\n";
print METAOUT "<P>Total Time: ".$TotTime." Seconds</P>\n";
$MyFinTime = &MyTimeStamp();
print METAOUT "<P><I>Analysis Complete:<BR>".$MyFinTime."</I></P>";
print METAOUT "</BODY>\n</HTML>";
close METAOUT;

#-----------------------------+
# CLOSE THE CONTIG ID FILE    |
#-----------------------------+
close CONOUT;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub MyTimeStamp
{
#-----------------------------+
# GET A TIMESTAMP USING       |
# THE LOCALTIME FUNCTION      |
#-----------------------------+
# Modified from:
# http://perl.about.com/od/perltutorials/a/perllocaltime_2.htm
    my ($second, $minute, $hour, $dayOfMonth, $month,
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings);
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    ($second, $minute, $hour, $dayOfMonth, $month, 
     $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year"; 
    
   return $theTime;
}

sub RunTGICL
{
#-----------------------------+ 
# RUN TGICL SUBFUNCTION       |
#-----------------------------+
    my $TFile = $_[0];          # Fasta file to cluster with tgicl
    
    my $cmd = "tgicl ". $TFile." ".$Topt;
    system ($cmd);   

    # Print TGICL options to the metadata file
    #print METAOUT "<P>TGICL Options: ".$Topt."</P>\n";

}

sub Blast2Gff
{
#-----------------------------+
# CONVERT BLAST OUTPUT TO A   |
# GFF FILE FORMAT FOR         |
# VISUALIZATION IN APOLLO     |
#-----------------------------+
# 
# This requires a global scale connection to
# a dbh database object.

    # OUTPUT IN FORM
    #LTR/gypsy	RepeatMasker: AllRep	LTR/gypsy	27676	27772	268

    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my $CurCon = $_[0];
    my $BlastFile = $_[1];
    my $GFFOutDir = $_[2];

    my $GFFOut = $GFFOutDir."/".$CurCon.".gff";

    #print color 'red';
    #print $GffOut."\n";
    #print color 'reset';

    #-----------------------------+
    # BLAST RELATED VARS          |
    #-----------------------------+
    my $MaxE = 0.00001;
    my $MinQryLen = 100;

    #-----------------------------+ 
    # SET LOCAL SCOPE             |
    #-----------------------------+ 
    my $HitName;  
    my $BlastResult;
    my $BlastHit;
    my $BlastHSP;

    #-----------------------------+
    # OPEN BLAST REPORT           |
    #-----------------------------+
    my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
					  '-file'   => $BlastFile,
					  '-signif' => $MaxE,
					  '-min_query_len' => $MinQryLen) 
	|| die "Could not open BLAST input file:\n$BlastFile.\n";
    
    #-----------------------------+
    # OPEN THE GFF OUTPUT FILE    |
    #-----------------------------+
    open (GFFOUT, ">$GFFOut") 
	|| die "Can not open file:\n $GFFOut\n";
    
    #-----------------------------+
    # PRINT LINE TO THE GFF FILE  |
    # FOR EVERY HSP IN THE BLAST  |
    # OUTPUT FILE                 |
    #-----------------------------+
    while ($BlastResult = $BlastReport->next_result())
    {
    	while ( $BlastHit = $BlastResult->next_hit())
	{
	    while ($BlastHSP = $BlastHit->next_hsp())
	    {
		
		$HitName = $BlastHit->name();
		my @TmpSplit = split( /\|/, $HitName);
		my $HitIdNum = $TmpSplit[0] || die "Split did not work";
		#my $Bac = &GetBacName( $HitName, $Species);
		my $Bac = &GetBacName( $HitIdNum, $Species);
		
		print GFFOUT $BlastHit->name()."\t".
		    "$Bac\t".
		    "$Bac\t".
		    $BlastHSP->start('query')."\t".
		    $BlastHSP->end('query')."\t".
		    $BlastHSP->score()."\t".
		    ".\t.\t".
		    "$Bac\n";
		
	    } # End of while next hsp
	    
	} # End of while next hit
    } # End of while next result

    close GFFOUT;

    # Return the full path of the GFF file when this works
    return $GFFOut;

} # END OF Blast2Gff Subfunction


sub GetBacName
{
#-----------------------------+
# GET THE BAC ID FOR THE      |
# ID NUMBER AND SPECIES NAME  |
# PASSED TO THE SUBFUNCTION   |
#-----------------------------+

    my $SpecNum = $_[0];
    my $SpecName = $_[1];
    my $result;

    my $SelectSQL = "SELECT bac FROM ".$tblSrcTable.
	" WHERE spec_num = '".$SpecNum."'".
	" AND species = '".$SpecName."'";
    
    #if (! $quiet)
    #{
    #	print "SQL STATEMENT\n$SelectSQL\n";
    #}

    $dbh->do($SelectSQL);
    $cur = $dbh->prepare($SelectSQL);
    $cur->execute();
    @row=$cur->fetchrow;
    $result=$row[0] || "SQL ERROR";
    $cur->finish();
    return $result;

}

sub AsmAnalysis
{
#-----------------------------+
# ASSEMBLY ANALYSIS           |
# SUBFUNCTION                 |
#                             |
#-----------------------------+
# FOR A CLUSTER FILE PRODUCDED|
# BY TGICL, CREATE A DIR FOR  |
# EACH CLUSTER AND RUN BLAST  |
# AGAINST THE INITIAL DB      |
# AND THE SET OF REPEAT DBS   |
#-----------------------------+

    # This assumes that all output is in the asm_1 directory
    my $PanDir = $_[0];        # Directory that holds the PAN Data
    my $PanId = $_[1];         # The unique ID of the PAN

    my $ConFile = $PanDir."asm_1/contigs";  #
    my $AceFile = $PanDir."asm_1/ACE";      # The ACE file produced by tgicl
    my $ClustFile = $PanDir.$PanId.".fasta_clusters";
    my $NewAce = $PanDir.$PanId.".ace";     # New version of ACE file with .ace
                                            # extension for opening with clview.
    my @AceMetaInfo;           # Array to hold meta data parsed from 
                               # the ace file to pass when working with contig file

    #-----------------------------------------------------------+
    # ACE FILE                                                  |
    #-----------------------------------------------------------+
    if (-e $AceFile)
    {

	# This uses the perl File::Copy module
	copy ($AceFile, $NewAce);
	# Provide a link to the Ace file in the metadata html file
	print METAOUT "<A HREF=".$NewAce.">".$PanId.".ace</A>\n";

	# WORK THROUGH THE ACE FILE TO GET INFORMATION 
	# ABOUT THE CLUSTERS
	
	my $ConNum = 0; 
	my ($seq_in, $seq_out);
	my $ConPad; 
	my $AceConId;
	my @AceSplit;
	my $AceLen;
	my $AceNumSeq;

	open (ACEIN, $AceFile) || 
	    die "Could not open\n$AceFile\n";

	while (<ACEIN>)
	{	    
	    my $FlagStr = substr($_,0,2);
	    
	    # Commented out 07/27/2006
	    #if ($FlagStr =~ "AS")
	    #{
		#$ConNum++;                 # The contig number
		## Pad the contig number with leading zeros
		#$ConPad = sprintf("%0${pad_len}d", $ConNum);
		#$AceConId = $PanId."CON".$ConPad;
	    #}

	    #Now can try to get the contig information such as
	    #length and number of seqs from the ace file itself
	    if ($FlagStr =~ "CO")
	    {
		
		# Added these here 07/27/2006
		$ConNum++;                 # The contig number
		$ConPad = sprintf("%0${pad_len}d", $ConNum);
		$AceConId = $PanId."CON".$ConPad;

		@AceSplit = split ( / /, $_ );
		my $AceLen = $AceSplit[2];
		my $AceNumSeq = $AceSplit[3];
		$AceMetaInfo[0][$ConNum] = $AceLen;
		$AceMetaInfo[1][$ConNum] = $AceNumSeq;
		$AceMetaInfo[2][$ConNum] = $AceConId;
	    }

	    if ($FlagStr =~ "AF")
	    {
		@AceSplit = split( / /, $_ );  # Split by white space
		$AceSeqId = $AceSplit[1] || "ERR";
		print CONOUT $AceSeqId."=".$AceConId."\n";
	    }
	}

	close ACEIN;       # Close the ACE file
    }

    #-----------------------------------------------------------+
    # CONTIG FILE                                               |
    #-----------------------------------------------------------+
    my $ConNum = 0;             # Counter for number of clusters within the pan
    my ($seq_in, $seq_out);     # Set scope of seq objects
    my $ConPad;                 # Set scope for padded contig number

    if (-e $ConFile)
    {
	print color 'green';
	print "Found contig file:\n$ConFile\n";
	print color 'reset';

	$seq_in  = Bio::SeqIO->new (
				    '-format' => 'fasta',
				    '-file' => "<".$ConFile );


	# START TABLE IN METADATA HTML FILE
	print METAOUT "\n<TABLE border=1>\n";
	print METAOUT "<TR>".
	    "<TD align=center><B>Contig</B><BR>".
	    "Database</TD>".
	    "<TD align=center><B>FastLen</B><BR>".
	    "Name</TD>".                      #Name
	    "<TD align=center><B>AceLen</B><BR>".
	    "Class</TD>".                     #Class
	    "<TD align=center><B>NumSeq</B><BR>".
	    "Hits</TD>".                      #Hits
	    "<TD align=center><B>FASTA</B><BR>".
	    "Raw</TD>".
	    "<TD align=center><B>Apollo</B><BR>".
	    "EVal</TD>".
	    "<TD align=center><B>NA</B><BR>".
	    "Length</TD>".
	    "</TR>\n";

	#-----------------------------+
	# FOR EVERY SEQEUNCE IN THE   |
	# FASTA FILE                  |
	#-----------------------------+
	while( my $seqobj = $seq_in->next_seq() )   
	{
	    $TotConCount++;            # Total num of contigs in the analysis
	    $ConNum++;                 # The contig number
	    # Pad the contig number with leading zeros
	    $ConPad = sprintf("%0${pad_len}d", $ConNum);

	    # Output path of the individual contig
	    my $OutPath = $PanDir.$PanId."CON".$ConPad."/".
		$PanId."CON".$ConPad.".fasta";

	    #Output Dir for individual Contig
	    my $OutDir =  $PanDir.$PanId."CON".$ConPad."/";
	    mkdir $OutDir, 0777 unless (-e $OutDir); # set permissions

	    # Show the contig number and output path
	    my $CurCon = $PanId."CON".$ConPad;
	    print color 'bold';
	    print "\t$CurCon\n";
	    print color 'reset';

	    #-----------------------------+
	    # OUTPUT RECORD AS FASTA FILE |
	    #-----------------------------+
	    $seq_out =  Bio::SeqIO->new (
					 '-format' => 'fasta',
					 '-file' => ">".$OutPath ) ||
					 die "Count not open $OutPath\n";

	    # Reset name to the CPAN.CON name
	    my $AsmLen = $seqobj->length;  # Length of the Assembly
	    #print METAOUT "<TR>\n".
		#"<TD>".$CurCon."</TD>\n".
		#"<TD><A HREF=".$OutPath.">FASTA</A></TD>\n".
		#"<TD>".$AsmLen."</TD>\n".
		#"</TR>\n";
	    $seqobj->display_id($CurCon);
	    $seq_out->write_seq($seqobj);  

	    
	    #-----------------------------+
	    # ADD RECORD TO THE FASTA FILE|
	    # OF ALL PANS                 |
	    #-----------------------------+
	    $all_seq_out->write_seq($seqobj);

	    #-----------------------------+
	    # RUN THE ANALYSIS ON THE     |
	    # ASSEMBLY FASTA FILE         |
	    #-----------------------------+
	    if (-e $OutPath)
	    {

		# BLAST ASSEMBLY AGAINST THE INDIVIDUAL READS
		# THAT WERE THE SOURCE SEQUENCES FOR THE
		print color 'blue';
		print "\tBLASTING $CurCon against source seqs\n";
		print color 'reset';
		#&SrcBlast ( $CurCon, $OutPath, $OutDir, $SourceDbName, $SourceDb );
		# I made the SrcBlast subfunction return the name
		# of the resulting blast outoput file so that
		# this can be passed to the Blast2Gff subfun
		my $SrcBlastOut = &SrcBlast ( $CurCon, $OutPath, $OutDir, 
					      $SourceDbName, $SourceDb );

		# MODIFY THE SOURCE BLAST TO THE GFF FORMAT
		# NEEDED FOR APOLLO
		print color 'blue';
		print "\tConverting BLAST Output to GFF.\n";
		print color 'reset';
		my $GffOut = &Blast2Gff ( $CurCon, $SrcBlastOut, $OutDir );
		
		# USE THE APOLLO FUNCTIONS TO OPEN THE GFF FILE WITH
		# THE FASTA FILE IN APOLLO AND MAKE AN APOLLO BACKUP FILE
		print color 'blue';
		print "\tConverting GFF to APOLLO FORMAT.\n";
		print color 'reset';
		my $ApolloOutPath = $OutDir.$CurCon.".apollo";
		&ApolloConvert ( $GffOut, 'gff', $ApolloOutPath, 'backup', 
				 $OutPath, 'na' );


		# Move down to here to include the apollo format file
		print METAOUT "<TR>\n".
		    "<TD bgcolor=#CCCCCC>".$CurCon."</TD>\n".
		    "<TD bgcolor=#CCCCCC align=center>".$AsmLen."</TD>\n".
		    "<TD bgcolor=#CCCCCC align=center>".$AceMetaInfo[0][$ConNum]."</TD>\n". # Length
		    "<TD bgcolor=#CCCCCC align=center>".$AceMetaInfo[1][$ConNum]."</TD>\n". # NumSeqs
		    "<TD bordercolor=black bgcolor=#CCCCCC><A HREF=".$OutPath.">FASTA</A></TD>\n".
		    "<TD bordercolor=black bgcolor=#CCCCCC><A HREF=".
		    $ApolloOutPath.">Apollo</A></TD>".
		    #"<TD bgcolor=#CCCCCC>".$AceMetaInfo[2][$ConNum]."</TD>\n". # ConID
		    "<TD bgcolor=#CCCCCC align=center>\&nbsp</TD>\n".
		    "</TR>\n";		
		
		# BLAST ASSEMBLY AGAINST THE REPEAT DATABASES
		print color 'blue';
		print "\tBLASTING $CurCon against repeats\n";
		print color 'reset';
		&RepBlast ( $CurCon, $OutPath, $OutDir );
		
		#-----------------------------+
		# RUN HMM MODELS              |
		# - added 08/08/2006          |
		#-----------------------------+
		&RepHmmerRun ( $CurCon, $OutPath, $OutDir, "MITE" );
		&RepHmmerRun ( $CurCon, $OutPath, $OutDir, "MULE" );
		&RepHmmerRun ( $CurCon, $OutPath, $OutDir, "TPASE" ); # No hits in Cen
		
	    } # End of if the fasta file $OutPath exists

	} # End of while seq_obj next sequence

	print METAOUT "</TABLE>\n";
	

    }else{
	print color 'bold red';
	print "Could not find contig file:\n$ConFile\n";
	print color 'reset';

	print METAOUT "<P><FONT COLOR=#CD0000>Could not find contig file:<BR>\n";
	print METAOUT  "$ConFile</FONT></P>"
    } # End of if contig file exists



}


sub SrcBlast
{

#-----------------------------+
# BLAST ASSEMBLY AGAINST THE  |
# SOURCE SEQUENCES            |
#-----------------------------+
# Can be used to blast against an db at a single path
# Given a fasta file path
# and working dir, blast the fasta file
# against a database of 
    #-----------------------------+
    # GET VARS PASSED TO SUBFUN   |
    #-----------------------------+
    my $SrcIndQry = $_[0];          # The name of the fasta file 
                                    # This is the contig id
    my $SrcQryPath = $_[1];         # FASTA file with qry sequence
    my $SrcWorkDir = $_[2];         # Working dir to place BLAST output
    my $SrcDbName = $_[3];          # The source db name
                                    # This var is used to name the out file
    my $SrcDbPath = $_[4];          # The blast db that contains the src seqs
    
    #-----------------------------+
    # VAR SET IN SUBFUN           |
    #-----------------------------+
    my $BlSuf = "-e 0.00001 -a 2";                      # Blast suffix
    
    #-----------------------------+
    # SET VARIABLE SCOPE TO LOCAL |
    # AND INITIALIZE COUNTERS     |
    #-----------------------------+
    my ( $IndDb , $BlastCmd );
    my $SrcOutPath = $SrcWorkDir.$SrcIndQry."_".$SrcDbName.".blo";
    $BlastCmd = "blastall -p blastn -i $SrcQryPath -d".
	    " $SrcDbPath -o $SrcOutPath $BlSuf";
    print "\t".$SrcIndQry."_".$SrcDbName."\n";
    system ($BlastCmd); # Do the blast command 
    
    return $SrcOutPath;           # Return the source output path
}

sub RepBlast
{

#-----------------------------+
# REPEAT BLAST SUBFUNCTIONS   |
#-----------------------------+
# Given a fasta file path
# and working dir, blast the fasta file
# against a database of 
    #-----------------------------+
    # GET VARS PASSED TO SUBFUN   |
    #-----------------------------+
    my $RBIndQry = $_[0];          # The name of the fasta file 
                                   # This is the contig id
    my $RBQryPath = $_[1];         # FASTA file with qry sequence
    my $RBWorkDir = $_[2];         # Working dir to place BLAST output
    
    #-----------------------------+
    # VAR SET IN SUBFUN           |
    #-----------------------------+
    my $RBDbDir = "/home/jestill/blast/repeats/";    # Base dir for db sequences
    my $BlSuf = "-e 0.001 -a 2";                     # Blast suffix
    
    #-----------------------------+
    # SET VARIABLE SCOPE TO LOCAL |
    # AND INITIALIZE COUNTERS     |
    #-----------------------------+
    my @RepDb;                     # List of the names of the repeat databases
    my ( $IndDb , $BlastCmd );
    my $RBProcNum = 0;

    #-----------------------------+
    # BLAST REPEAT QRY DATABASES  |
    #-----------------------------+
    @RepDb = ( "gram_rep",         # Gramineae v3.1 repeat database from TIGR
	       "os_rep",           # Oryza sativa repeat database from TIGR
	       "RB_pln",           # All plants repeat database from RepBase
	       "TREP_8",           # TREP v.8 database
	       "zm_rep",           # Zea mays repeat database from TIGR
	       "SanMiguel",        # Phillip SanMiguel DB
	       "Wessler"           # Wessler lab MAGI database
	       );
    
    my $RBLenDb =  @RepDb;
    my $RBNumProc = $RBLenDb;
    
    # RUN THE BLAST FOR EVERY QUERY DATABASE 
    # IN THE QUERY DB FILE
    for $IndDb (@RepDb)
    {
	$RBProcNum++;
	print "\t$RBProcNum of $RBNumProc: ";
	print $RBIndQry."_".$IndDb."\n";
	my $RBDbPath = $RBDbDir.$IndDb;
	my $TestFile = $RBDbPath.".nhr";
	my $RBOutPath = $RBWorkDir.$RBIndQry."_".$IndDb.".blo";
	my $RBImgPath = $RBWorkDir.$RBIndQry."_".$IndDb.".gif";
	my $RBTblPath = $RBWorkDir.$RBIndQry."_".$IndDb.".tbl";
	#$BlastCmd = "blastall -p blastn -i $RBQryPath -d".
	#    " $RBDbPath -o $RBOutPath $BlSuf";

	$BlastCmd = "blastall -p ".$BlProg." -i $RBQryPath -d".
	    " $RBDbPath -o $RBOutPath $BlSuf";
	
	# DO THE BLAST
	system ($BlastCmd); # Do the blast command 

	# GENERATE A GIF IMAGE FROM THE BLAST OUTPUT
	# The blast-imager.pl program is currently not working properly
	#my $Bl2TblCmd = $ScriptDir."blast2table.pl ".$RBOutPath." > ".
	#    $RBTblPath;
	#my $BlImgCmd = $ScriptDir."blast-imager.pl ".$RBTblPath."> ".$RBImgPath;
	#system ($Bl2TblCmd); # Convert blast output to table format
	#system ($BlImgCmd);  # Convert blast output to GIF image
	
	&GetRepClass ($RBOutPath);            # Get blast information

	# ONCE THE BLAST IS COMPLETED, IT MAY BE USEFUL TO PARSE THE
	# BLAST OUTPUT AND DETERMINE THE REPEAT CLASS OF THE BEST HIT

    }

    
}


sub Cyto2Fasta
{
#-----------------------------+ 
# CONVERTS CYTOSCAPE TEXT     |
# WITH SEQUENCE FILES TO A    |
# FASTA FORMAT                |
#-----------------------------+

    #-----------------------------+
    # VARIABLES                   |
    #-----------------------------+
    my $InFilePath = $_[0];
    my $OutFilePath = $_[1];
    my $LineNum = 0;
    
    #-----------------------------+
    # FILE I/O CONNECTIONS        |
    #-----------------------------+
    open (INFILE, "<$InFilePath") ||
	die "Can not open INFILE:\n$InFilePath\n";
    
    open (OUTFILE, ">$OutFilePath") ||
	die "Can not open OUTFILE:\n$OutFilePath";

    #-----------------------------+
    # WORK WITH THE BLAST RESULT  |
    #-----------------------------+
    while (<INFILE>)
    {
	chomp;                    # Get rid of newline character
	
	$LineNum++;
	#print $LineNum/"\n";
	
	#WILL NEED TO IGNORE THE FIRST LINE (HEADER LINE)
	if ($LineNum > 1)
	{
	    my @Split = split(/\t/, $_ );
	    my $SeqId = $Split[0]  || "ERROR";
	    my $SeqStr = $Split[1]    || "ERROR";
	    
	    # Print Record ID TO Scree
	    my $RecNum = $LineNum - 1;
	    #print "Record: ".$RecNum."\tSeq: ".$SeqId."\n";
	    
	    # PRINT OUTPUT TO OUTPUT FILE
	    print OUTFILE ">".$SeqId."\n";
	    print OUTFILE $SeqStr."\n";

	}
	
    }
    
    #-----------------------------+
    # CLOSE THE FILES             |
    #-----------------------------+ 
    close INFILE;
    close OUTFILE;
}

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ApolloConvert
{
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my $InFile = $_[0];        # Input file path
    my $InForm = $_[1];        # Output file format:
                               # game|gff|gb|chadoxml|backup
    my $OutFile = $_[2];       # Output file path
    my $OutForm = $_[3];       # Ouput file foramt
                               # chadoDB|game|chadoxml|genbank|gff|backup
    my $SeqFile = $_[4];       # The path of the sequence file
                               # This is only required for GFF foramt files
                               # When not required this can be passed as na
    my $DbPass = $_[5];        # Database password for logging on to the 
                               # chado database for reading or writing.
    my ( $ApPath, $ApCmd );

    $ApPath = "/home/jestill/Apps/Apollo/";

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ApPath."Apollo -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}


sub GetWesClass
{
#-----------------------------+
# GETS THE TERMINAL           |
# CLASSIFICATION NAME FOR HITS|
# TO THE WESSLER REPEAT       |
# DATABASE FROM MAGI SITE     |
#-----------------------------+
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    my %WesClass;

    %WesClass = (
		 # RETROTRANSPOSONS
		 'LTR/copia' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/Copia' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/COPIA' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/copia?' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/gypsy' => 'Ty3-gypsy LTR-Retrotransposon',
		 'LTR/Gypsy' => 'Ty3-gypsy LTR-Retrotransposon',
		 'LTR/?' => 'LTR-Retrotransposon',
		 'LTR/unknown' => 'LTR-Retrotransposon',
		 'LTR/Barley' => 'LTR-Retrotransposon',
		 'LTR/Maize' => 'LTR-Retrotransposon',
		 'SINE/12-15bp?' => 'SINE Retrotransposon',
		 'LINE' => 'LINE Retrotransposon',
		 'LTR/Tto1_NTlike' => 'LINE Retrotransposon',
		 'LTR/Zea' => 'LTR-Retrotransposon',
		 # TRANSPOSONS
		 'DNA/8bp' => 'Transposon',
		 'DNA/hAT' => 'hAT Transposon',
		 'DNA/CACTA' => 'CACTA, En/Spm Transposon',
		 'DNA/spm' => 'CACTA, En/Spm Transposon',
		 'DNA/CACTG' => 'CACTA, En/Spm Transposon',
		 'DNA/TAA' => 'ping/pong/SNOOPY Transposon',
		 # FOLDBACK ELEMENTS
		 "MITE/NNN" => 'MITE',
		 'MITE/TA' => 'MITE',
		 'MITE/TA(stw)' => 'MITE',
		 # OTHER
		 'vector' => 'Vector',
		 'Vector/Bao' => 'Vector',
		 'unknown' => 'Unknown',
		 'Novel/?' => 'Unknown',
		 'unknown/Bao' => 'Unknown'
		 );

    $cl = $WesClass{$qry} || "UNK";
    return $cl;   
}

sub GetTREPClass
{
#-----------------------------+
# GET THE CLASSIFICATION OF   |
# AN ELEMENT GIVEN THE TREP   |
# ACCESSION ID                |
#-----------------------------+
# This will need to use tblTREPMeta table in the database
    
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned

    my $QryTbl = "tblTREP";

    my $SelQry = "SELECT j_class FROM ".$QryTbl.
	" WHERE accession = '".$qry."'";

    $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    @row=$cur->fetchrow;
    $cl = $row[0] || "TREP CLASS UNK";
    $cur->finish();

    return $cl;

}


sub GetTREPName
{
    
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    
    my $QryTbl = "tblTREP";

    my $SelQry = "SELECT name FROM ".$QryTbl.
	" WHERE accession = '".$qry."'";

    $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    @row=$cur->fetchrow;
    $cl = $row[0] || "TREP NAME UNK";
    $cur->finish();

    #$cl = "Currently Unknown";
    return $cl;

}

sub GetTIGRClass
{

    #-----------------------------+
    # PARSE TIGR REPEAT CODES     |
    # Parse the codes used by TIGR|
    # for Repeat names. Returns   |
    # SuperClass-Class-SubClass   |
    #-----------------------------+


    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    my %TIGRRep;            # Hash to translate from code to class
    my $code = substr ($qry,4,7);
    
    #-----------------------------+
    # HASH TO TRANSLATE FROM TIGR |
    # CODES TO REPEAT CLASSES     |
    #-----------------------------+
    %TIGRRep = (
		# RETROTRANSPOSONS
		TERT001 => 'Ty1-copia LTR-Retrotransposon',
		TERT002 => 'Ty3-gypsy LTR-Retrotransposon',
		TERT003 => 'LINE Retrotransposon',
		TERT004 => 'SINE Retrotransposon',
		TERTOOT => 'Unclassified Retrotransposon',
		# DNA TRANSPOSONS
		TETN001 => 'Ac/Ds Transposon',
		TETN002 => 'CACTA, En/Spm Transposon',
		TETN003 => 'Mutator (MULE) Transposon',
		TETN004 => 'Mariner (MLE) Transposon',
		TETN005 => 'ping/pong/SNOOPY Transposon',
		TETNOOT => 'Unclassified Transposon',
		# MITES
		TEMT001 => 'Tourist MITE',
		TEMT002 => 'Stowaway MITE',
		TEMT003 => 'Crackle MITE',
		TEMT004 => 'Explorer MITE',
		TEMT005 => 'Gaijin/Gaigin MITE',
		TEMT006 => 'Castaway MITE',
		TEMT007 => 'Snap MITE',
		TEMT008 => 'Amy/LTP MITE',
		TEMT009 => 'Ditto MITE',
		TEMT010 => 'Wanderer MITE',
		TEMT011 => 'p-SINE1 MITE',
		TEMT012 => 'Pop MITE',
		TEMT013 => 'Krispie MITE',
		TEMT014 => 'Snabo MITE',
		TEMT015 => 'MITE-adh, typeA MITE',
		TEMT016 => 'MITE-adh, typeB MITE',
		TEMT017 => 'MITE-adh, typeD MITE',
		TEMT018 => 'MITE-adh, typeG MITE',
		TEMT019 => 'MITE-adh, typeH MITE',
		TEMT020 => 'MITE-adh, typeI MITE',
		TEMT021 => 'MITE-adh, typeJ MITE',
		TEMT022 => 'MITE-adh, typeK MITE',
		TEMT023 => 'MITE-adh, typeL MITE',
		TEMT024 => 'MITE-adh, typeM MITE',
		TEMT025 => 'MITE-adh, tympN MITE',
		TEMT026 => 'MITE-adh-1',
		TEMT027 => 'MITE-adh-2',
		TEMT028 => 'MITE-adh-3',
		TEMT029 => 'MITE-adh-4',
		TEMT030 => 'MITE-adh-5',
		TEMT031 => 'MITE-adh-6',
		TEMT032 => 'MITE-adh-7',
		TEMT033 => 'MITE-adh-8',
		TEMT034 => 'MITE-adh-9',
		TEMT035 => 'MITE-adh-10',
		TEMT036 => 'MITE-adh-11',
		TEMT037 => 'MITE-adh-12',
		TEMT038 => 'Buhui MITE',
		TEMT039 => 'Casin MITE',
		TEMT040 => 'Centre MITE',
		TEMT041 => 'Delay MITE',
		TEMT042 => 'ECR MITE',
		TEMT043 => 'Helia MITE',
		TEMT044 => 'ID-2 MITE',
		TEMT045 => 'ID-3 MITE',
		TEMT046 => 'ID-4 MITE',
		TEMT047 => 'Lier MITE',
		TEMT048 => 'Stola MITE',
		TEMT049 => 'Stone MITE',
		TEMT050 => 'Susu MITE',
		TEMT051 => 'Wuji MITE',
		TEMT052 => 'Youren MITE',
		TEMT053 => 'Micron MITE',
		TEMT054 => 'Truncator MITE',
		TEMT055 => 'Heart Breaker MITE',
		TEMT056 => 'Frequent Flyer MITE',
		TEMT057 => 'Heart Healer MITE',
		TEMT058 => 'Acrobat MITE',
		TEMT059 => 'mPIF MITE',
		TEMT060 => 'Pangrangja MITE',
		TEMT061 => 'Kiddo MITE',
		TEMT062 => 'mPing/miniSNOOPY MITE',
		TEMT063 => 'MDM MITE',
		TEMTOOT => 'Unclassified MITE',
		# CENTROMERE RELATED SEQUENCES
		CMCM001 => 'Centromere Specific Retrotransposons',
		CMCM002 => 'Centromeric Satellite Repeats',
		CMCMOOT => 'Unclassified Centromere Sequence',
		# TELOMERE RELATED SEQUENCES
		TRTM000 => 'Telomere',
		TRTA000 => 'Telomere Associated',
		# RIBOSOMAL RNA GENES
		RGRR000 => '45s rDNA',
		RGRR005 => '5s rDNA',
		# UNCLASSIFIED
		OTOT000 => 'Unclassified'
		);

    $cl = $TIGRRep{$code} ||
	"UNKNOWN TIGR CODE";
    return $cl;

}

sub GetRBClass
{
#-----------------------------+
# GET REPEAT BASE CLASS       |
#-----------------------------+

# A standardized ontology for the class is returned given the name  
# used by RepBase for the repeat class.

    my $qry = $_[0];           # The qry code sent to the subfunction
    my $cl;                    # The classification that is returned
    my %RBClass;               # Hash to translate from code to class

    %RBClass = (
		#-----------------------------+
		# INTERSPERSED REPEATS        |
		#-----------------------------+
		Interspersed => 'Interspersed Repeat',
		#DNA TRANSPOSONS
		Mariner => 'Mariner (MLE) Transposon',
		hAT => 'DNA Transposon',
		MuDR => 'DNA Transposon',
		EnSpm => 'CACTA, En/Spm Transposon',
		piggyBac => 'DNA Transposon',
		P => 'DNA Transposon',
		Merlin => 'DNA Transposon',
		Harbinger => 'DNA Transposon',
		Transib => 'DNA Transposon',
		Novosib => 'DNA Transposon',
		Mirage => 'DNA Transposon',
		Helitron => 'Helitron',
		Polinton => 'DNA Transposon',
		Rehavkus => 'DNA Transposon',
		DNA => 'DNA Transposon',
		#LTR Retrotransposon
		Gypsy => 'Ty3-gypsy LTR-Retrotansposon',
		Copia => 'Ty1-copia LTR-Retrotransposon',
		LTR => 'LTR Retrotransposon',
		BEL => 'LTR Retrotransposon',
		DIRS => 'LTR Retrotransposon',
		#Endogenous Retrovirus
		ERV1 => 'Endogenous Retrovirus',
		ERV2 => 'Endogenous Retrovirus',
		ERV3 => 'Endogenous Retrovirus',
		#Non-LTR Retrotransposon
		SINE => 'SINE Retrotransposon',
		SINE1 => 'SINE Retrotransposon',
		SINE2 => 'SINE Retrotransposon',
		SINE3 => 'SINE Retrotransposon',
		CRE => 'Non-LTR Retrotransposon',
		NeSL => 'Non-LTR Retrotransposon',
		R4 => 'Non-LTR Retrotransposon',
		R2 => 'Non-LTR Retrotransposon',
		L1 => 'Non-LTR Retrotransposon',
		RTE => 'Non-LTR Retrotransposon',
		I => 'Non-LTR Retrotransposon',
		Jockey => 'Non-LTR Retrotransposon',
		CR1 => 'Non-LTR Retrotransposon',
		Rex1 => 'Non-LTR Retrotransposon',
		RandI => 'Non-LTR Retrotransposon',
		Penelope => 'Non-LTR Retrotransposon',
		# Caulimoviridae 
		Caulimoviridae => 'Caulimoviridae', 
		#-----------------------------+
		# SIMPLE REPEAT               |
		#-----------------------------+
		SAT => 'Satellite',
		MSAT => 'Satellite'
		);

    $cl = $RBClass{$qry} ||
	"UNKNOWN REPBASE NAME";
    return $cl;

}


sub GetRepClass
{

    #-----------------------------+
    # LOAD THE REPEAT             |
    # CLASSIFICATION BLAST OUTPUT |
    #-----------------------------+

    my $RetClass;                  # The class name that can be returned
    #my $HitCount = '0';            # Initialize hit count to zero to fetch best 
                                   # hit later

    my $RepBLAST = $_[0];          # The repeat blast file to parse
    #my $TotHits = $_[1];            # The total hits recorded to date 
    #                               # for the query sequence

    my $LoadRec;
    my $QryName;
    my $RepClass;
    
    my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
					  '-file'   => $RepBLAST,
					  '-signif' => $MaxE,
					  '-min_query_len' => $MinQryLen,
					  '-score' => $MinScore ) ||
					      die "Could not open BLAST input file:\n$RepBLAST.\n";
    
    #-----------------------------+
    # SINCE THERE IS ONLY ONE FILE|
    # SENT TO BLAST THERE SHOULD  |
    # BE ONLY ONE RESULT          |
    #-----------------------------+
    while ($BlastResult = $BlastReport->next_result())
    {
	$QryName = $BlastResult->query_name;
	$BlastDB = $BlastResult->database_name;

	my $NumHits = $BlastResult->num_hits;
	my $HitCount = '0';

	#-----------------------------+
	# THERE COULD BE MULTIPLE HITS|
	# FOR EVERY FILE, THIS WILL   |
	# JUST SHOW THE HITS IN THE DB|
	# SENT TO THE SUBFUNCTION     |
	#-----------------------------+
	while ( $BlastHit = $BlastResult->next_hit())
	{

	    #-----------------------------+
	    # TREP REPEAT DATABASE        |
	    #-----------------------------+
	    if ($BlastDB =~ 'TREP_8')
	    {
		$HitCount++;

		my $HitId = $BlastHit->accession()  || "UnkAcc";

		$RepClass = &GetTREPClass($HitId);
		$Name = &GetTREPName($HitId);

		# Commented out 08/10/2006
		# Don't show output on quiet runs
		#if (! $quiet)
		#{
		#    print $QryName."\n";
		#    print "\t".$HitId.":".$Name."\n";
		#    print "\t".$RepClass."\n";
		#}

	    }

	    #-----------------------------+
	    # REPBASE PLANTS              |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'RB_pln')
	    {

		$HitCount++;
		#$Class = "UNK";
		#$Subclass = "UNK";
		#$Superfamily = "UNK";
		$Name = $BlastHit->name();
		
		my @SpName = (split /\#/, $Name);
		my $CatSearch = trim($SpName[1] || "NONE");
		
		$RepClass = &GetRBClass($CatSearch); 

	    }
	    #-----------------------------+
	    # TIGR ORYZA REPEAT DATABASE  |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'os_rep')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name);
	    }

	    #-----------------------------+
	    # WESSLER LAB REPEAT DATABASE |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'Wessler')
	    {
		$HitCount++;
		my $WesName = $BlastHit->name();
		my @SpName = (split /\#/, $WesName);
		$Name = $SpName[0];
		
		my $WesClass = $SpName[1];

		#$Desc = $BlastHit->description();
		$RepClass = &GetWesClass($WesClass);

		# DEBUG PRINT
		# Not printed on 'quiet' runs
		if (! $quiet)
		{

		    if ($RepClass  =~ "UNK")
		    {
			print $Name."\n";
			print "\t".$WesClass."\n";
			print "\t".$RepClass."\n";
		    }
		}

	    }

	    #-----------------------------+
	    # SAN MIGUEL REPEAT DATABASE  |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'SanMiguel')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = "UNK-SanMiguel";
	    }

	    #-----------------------------+
	    # TIGR ZEA MAYS REPEAT        |
	    # DATABASE                    |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'zm_rep')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name); 
	    }

	    #-----------------------------+
	    # TIGR GRAMINEAE REPEAT       |
	    # DATABASE                    |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'gram_rep')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		#$RepClass = "UNK-ZeaRepeat";
		$RepClass = &GetTIGRClass($Name); 
	    }

	    #-----------------------------+
	    # CALC THE LENGTH OF HIT      |
	    #-----------------------------+
	    my $TotHitLen = "0";
	    while ( my $BlastHSP = $BlastHit->next_hsp())
	    {
		$TotHitLen = $TotHitLen + $BlastHSP->length('total');
	    }

	    # IF the all hits flag was selected, show 
	    # all blast hits otherwise just show the first/best hit
	    if ($all)
	    {
		print METAOUT "<TR>".
		    "<TD align=right>".$BlastDB."</TD>".      # Repeat Database
		    "<TD>".$Name."</TD>".         # Repeat class
		    "<TD>".$RepClass."</TD>".     # Repeat name
		    "<TD>".$NumHits."</TD>".      # Number of total hits
		    "<TD>".$BlastHit->raw_score."</TD>".
		    "<TD>".$BlastHit->significance."</TD>".
		    "</TR>\n";
	    }else{
		if ($HitCount == "1")
		{
		    
		    # May just temporarily print the output to the metadata table
		    print METAOUT "<TR>".
			#"<TD align=right>".$BlastDB."</TD>".     # Repeat Database
			"<TD align=right><A HREF=".$RepBLAST.">"  # Link to blastoutput
			.$BlastDB."</A></TD>".
			"<TD>".$Name."</TD>".                    # Repeat class
			"<TD>".$RepClass."</TD>".                # Repeat name
			"<TD align=center>".$NumHits."</TD>".    # Number of total hits
			"<TD align=center>".$BlastHit->raw_score."</TD>". 
			"<TD align=center>".$BlastHit->significance."</TD>".  # E value
			#"<TD align=center>".$BlastHit->length."</TD>".  # length of hit sequence
			"<TD align=center>".$TotHitLen."</TD>";  # length of hit sequence
			#"<TD>LOC:".$BlastHit->locus."</TD>";         # Locus Name
		    if ($desc)
		    {
			print METAOUT "<TD>BIT:".$BlastHit->bits."</TD>".          # Bit score
			    "<TD>LEN:".$BlastHit->length."</TD>".        # Length
			    "<TD>ACC:".$BlastHit->accession."</TD>".     # Accession ID
			    "<TD>HSP:".$BlastHit->hsps."</TD>".          # Num hsps
			    "<TD>DES:".$BlastHit->description."</TD>";
		    }
		    
		    print METAOUT "</TR>\n";
		}
	    }
	    
	} # End of while BlastResult-next_hit

    } # End of while BlastReport next_result

    
} # End of LoadRepClass subfun


sub trim($)
	 {
	     my $string = shift;
	     $string =~ s/^\s+//;
	     $string =~ s/\s+$//;
	     return $string;
	 }
# Left trim function to remove leading whitespace

sub ltrim($)
	  {
	      my $string = shift;
	      $string =~ s/^\s+//;
	      return $string;
	  }
# Right trim function to remove trailing whitespace

sub rtrim($)
	  {
	      my $string = shift;
	      $string =~ s/\s+$//;
	      return $string;
	  }

sub RepHmmerRun
{
    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my $PanName = $_[0];           # The name of the individual query ie PPAN058CON001
    my $QryPath = $_[1];           # The path to to the qry fasta file
    my $WorkDir = $_[2];           # The work dir,
                                   # Will need to make a hmm_class dir here
    my $class = $_[3];             # The class of repeast to search
                                   # vars can include MITE,MULE

    #-----------------------------+
    # VAR SCOPE AND INITIALIZE    |
    #-----------------------------+
    $ProcNum = 0;                  # Process number starts at zero
    my ( $HmmCmd, $DbPath, $OutPath );  # Declare scope for varaiables used later
		  
    #-----------------------------+
    # HMMER MODELS                |
    #-----------------------------+    
    # FIRST DETERMINE THE APPROPRIATE DIR GIVEN
    # THE MODEL SET THAT IS BEING SEARCHED (MITE/MULE)
    if ($class =~ "MITE")
    {
	$ModDir = "/home/jestill/HMMData/db/hmm/mite_models/";
    }elsif ($class =~ "MULE"){
	$ModDir = "/home/jestill/HMMData/db/hmm/mule_models/";
    }elsif ($class =~ "TPASE"){
	$ModDir = "/home/jestill/HMMData/db/hmm/tpase_models/";
    }elsif ($class =~ "PFAM"){
	$ModDir = "/home/jestill/HMMData/db/hmm/pfam/";
    }else{
	$ModDir = "/home/jestill/HMMData/db/hmm/mite_models/"; # Default is to just use MITES
    }
    
    # Open the appropriate dir and load files 
    opendir( MODDIR, $ModDir );
    my @Mod = grep !/^\.\.?$/, readdir MODDIR ;
    closedir( MODDIR );    

    #-----------------------------+
    # CREATE A HMMER OUTPUT DIR   |
    #-----------------------------+
    my $HmmOutDir = $WorkDir.$class;
    mkdir $HmmOutDir, 0777 unless (-e $HmmOutDir); # set permissions
    
    # Determine the total of HMMER queries that will be run
    # Since a single sequence is passed to the subfun this
    # will just be the same as the number of models to test
    # against.
    my $LenMod =  @Mod;
    my $NumProc = $LenMod;
    
    for $IndMod (@Mod)
    {
	
	$ProcNum++; # Increment the process number
	print "HMM ".$class." Process ".$ProcNum." of ".$NumProc."\n";
	
	# HMMER MODEL PATH
	$ModPath = $ModDir.$IndMod;
	$ModFile = $ModPath;

	$OutPath = $HmmOutDir."/".$PanName."_".$IndMod.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if (-e $ModFile ) 
	{
	    #print "DB: $IndDb exists\n";
	}
	else 
	{die "Can not find model:\n$IndMod\n"; }
	
	#-----------------------------+
	# DOES THE HMM QRY SEQ EXIST  |
	#-----------------------------+
	if (-e $QryPath)
	{
	    #print "QRY: $QryPath exists\n";
	}
	else
	{die "Can not find qry file:\n$QryPath\n";}
	
	#------------------------------+
	# PRINT THE HMMER COMMAND      |
	#------------------------------+
	$HmmCmd = "hmmsearch " . 
	    "--domT 2 $ModFile $QryPath >$OutPath";
	#print wrap("\t", "\t", $HmmCmd );
	#print "\n";
	
	#------------------------------+
	# RUN THE BLAST COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
	    system ($HmmCmd);
	}
	
	
    } # End of for each database loop


    # ONCE THE ENTIRE SET HAS RUN IT IS TIME TO PARSE
    print "PARSING THE HMM_".$PanName."_".$class." OUTPUT\n";
    &RepHmmerParse ( $PanName, $WorkDir, $HmmOutDir, $class, $QryPath );

} # END OF RepHmmerRun

sub RepHmmerParse
{
#-----------------------------+
# PARSE THE OUTPUT FROM A     |
# HMMER RUN AGAINST A SET OF  |
# HMM PROFILES                |
#-----------------------------+ 

    #-----------------------------+
    # VARIABLES                   |
    #-----------------------------+
    # THE FOLLOWING SHOULD BE PASSED TO THE SUBFUNCTION
    my $PanName = $_[0];
    my $WorkDir = $_[1];
    my $HmmDir = $_[2];
    my $class = $_[3];
    my $QrySeqFile = $_[4];  # The fasta file path for the query sequence
                             # this can be used to fetch the sequence of the
                             # putative MITE/MULE. It should be possible to
                             # get this from the 

    # File to write the parsed output to
    my $HmmOutFile = $WorkDir.$PanName."_".$class.".txt";

    if ($class =~ "MITE")
    {
	$dataset = "hmm_mite";
    }
    elsif ($class =~ "MULE")
    {
	$dataset = "hmm_mule";
    }    
    elsif ($class =~ "PFAM")
    {
	$dataset = "pfam";
    }   
    elsif ($class =~ "TPASE")
    {
	$dataset = "tpase";
    }else{
	$dataset = "UNK";
    }   

    my $Len = "";                  # Length of the hit
    my $BestName = "";             # Name of the best hit
    my $FilePath;                  # The file path for the inidividual HMM output
    my $TotRes = '0';              # Total number of results for the sequence
    my $BestBit = '0';             # Best BIT Score
    my $BestHit;                   # The name of the Best Hit
    my $BestEval = '';             # Initialize BestEval to NUL
    my $BestLen = '';              # Length of the hit with best bit score

    #-----------------------------+
    # LOAD FILES TO PARSE INTO    |
    # THE @HmmFiles ARRAAY        |
    #-----------------------------+
    opendir( HMMDIR, $HmmDir );
    my @HmmFiles = grep !/^\.\.?$/, readdir HMMDIR ;
    closedir( HMMDIR );
    
    # Open up an Output file
    open ( OUT, ">".$HmmOutFile);
    print OUT "SEQNAME      \tSTART\tEND\tSCORE\tEVAL\tHITNAME\tLENGTH\tSEQ\n";
    print OUT "=============\t=====\t===\t=====\t====\t=======\t======\t===\n";
    
    #-----------------------------------------------------------+
    # FOR EACH FILE IN THE ARRAY OF FILES                       |
    #-----------------------------------------------------------+
    for $IndFile (@HmmFiles)
    {
	
	$FilePath = $HmmDir."/".$IndFile;
	
        # OPEN THE HMM RESULT AS HMMER RESULTS OBJECT
	my $HmmRes = new Bio::Tools::HMMER::Results ( -file => $FilePath ,
						      -type => 'hmmsearch') 
	    || die "Could not open file\n$FilePath\n";
	
	my $NumRes = $HmmRes->number;
	$TotRes = $TotRes + $NumRes;
	
	# ONLY PRINT OUTPUT FOR QUERIES WITH MATCHES
	if ($NumRes >> 0) #08/07/2006
	{	              #08/07/2006
	    foreach $seq ( $HmmRes->each_Set ) 
	    {
		foreach $domain ( $seq->each_Domain ) 
		{		
		    #my $CurBit = $domain->bits;
		    my $CurName = $domain->hmmname;
		    $CurName =~ m/.*\/(.*)\.hmm/;   # Returns what I want
		    $CurName = $1;
		    # RECORD THE NAME AND SCORE OF THE
		    # BEST HIT
		    if ($domain->bits > $BestBit)
		    {
			# ASSUMES BIT SCORE AND E VALUE
			# HAVE THE SAME RANK
			$BestBit = $domain->bits;
			$BestName = $CurName;
			$BestEval = $domain->evalue;
			$BestLen = $domain->end - $domain->start;
		    }
		    
		    print OUT $seq->name."\t";
		    print OUT $domain->start."\t";
		    print OUT $domain->end."\t";
		    print OUT $domain->bits."\t";
		    print OUT $domain->evalue."\t";
		    print OUT $CurName."\t";
		    my $ModLen = $domain->end - $domain->start;
		    print OUT $ModLen."\t";
		    # Get the sequence of the putative MITE
		    # This assumes that the only record in the
		    # fasta file is the record of interest.
		    $HmmQry_seq_in  = Bio::SeqIO->new (
						       '-format' => 'fasta',
						       '-file' => "<".$QrySeqFile ) 
			|| die "Can not open file $QrySeqFile";
		    while( my $Hmm_seqobj = $HmmQry_seq_in->next_seq() )   
		    {
			# Get the repeat object string and print to the 
			# outfile
		    	$HmmRepStr = $Hmm_seqobj->subseq($domain->start,$domain->end);
			print OUT $HmmRepStr;
		    }
		    $HmmQry_seq_in=""; # Reset obj to null
		    print OUT "\n";

		} # End of for each domain in the HMM output file
	    } # End of if greater then zero hits 08/07/2006
	    
	} # End of for each seq in the HMM output file
	
    }
    
    print "\n\n===================================\n";
    print "$class HMMER RESULT\n";
    print "===================================\n";
    print "PAN:       \t".$PanName."\n";
    print "DATASET:   \t".$dataset."\n";
    print "CLASS:     \t".$class."\n";
    print "TOTAL:     \t". $TotRes."\n";
    print "BEST HIT:  \t".$BestName."\n";
    print "BEST BIT:  \t".$BestBit."\n";
    print "BEST EVAL: \t".$BestEval."\n";
    print "===================================\n";

    
    # PRINT APPROPRIATE OUTPUT TO THE METAFILE IF
    # HITS WERE FOUND IN THE DATABASE
    if ($TotRes > 0 )
    {
	print METAOUT "<TR>".
	    "<TD align=right><A HREF=$HmmOutFile>".$dataset.
	    "</A></TD>".      # Repeat Database
	    "<TD>".$BestName."</TD>".                 # Repeat name
	    "<TD>".$class."</TD>".                    # Repeat class
	    "<TD align=center>".$TotRes."</TD>".      # Number of total hits
	    "<TD align=center>".$BestBit."</TD>".
	    "<TD align=center>".$BestEval."</TD>".
	    "<TD align=center>".$BestLen."</TD>".
	    "</TR>\n";
    }
    
    
    close OUT;

    #-----------------------------------------------------------+
    # TO DO
    #-----------------------------------------------------------+
    # - Figure out if there are multiple MITES on a single
    #   contig. This is very likely given the short length
    #   of MITES. Currently just have to manually look at the
    #   output.
 

} # END OF THE REP HMMER PARSE SUBFUNCTION






=head1 HISTORY
The history of program development.
=cut
#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/17/2006
# - Started program
# - Cyto2Fasta subfunction
# - Started AsmAnalysis subfunction
#
# 07/21/2006
# - Additional work to the AsmAnalysis subfucntion
# - Added color text to simplify reading the output
#
# 07/24/2006
# - Additional work to the AsmAnalysis subfunction
# - Added RepBlast subfunction for use in the 
#   AsmAnalysis subfunction
# - Added SrcBlast subfunction
# - Added Blast2Gff subfunction
# - Added GetBacName subfunction
# - Added ApolloConvert subfunction from the
#   RepMaskParse9 program. 
#
# 07/25/2006
# - Copy the asm ace file to ClustName.ace in the
#   main output directory.
# - Added metadata html file with links to files
#   generated from the program.
# - Will parse BLAST to get the best hit 
#   from different databases for the blast output
#   from the repeat blast.
# - Added GetRBClass subfun from jabablast
# - Added GetTIGRClass subfun from jabablast
# - Added GetTREPName subfun from jabablast
# - Added GetTREPCls subfun from jabablast
# - Added GetWesClass subfun from jabablast
# - Added LoadRepClass subfun from jabablast
#   and modified for use in this format
# - Added trim function
# - Added ltrim function
# - Added rtrim function
#
# 07/26/2006
# - Cleaned up some code
# - Added variables to the metadata HTML output
# - Added bit score and e-value to the best
#   hit summary
# - Added capture of variables from the ace format
#   to the AceMetaData array
#
# 07/27/2006
# - Printng @AceMetaData information to the
#   MetaData table
# - Added MyTimeStamp SubFun
# - Changed format of tables showing output
# - Added variable for all hits to be shown
#   in the MetaData HTML file
# - Added variable for more descriptive information
#   to be shown for hits in MetaData HTML file
# - Added link to blast output when clicking on
#   the database name in the MetaData HTML file
#
# 07/28/2006
# - Adding a visualization of the BLAST output
#   using the simple blast2table.pl and the
#   blast-imager.pl scripts from Oreilly
# - This attempt did not work out
#
# 7/21/2005
# - Added ability to include images of the PAN if they
#   exist as a *.png file in the PAN directory.
#
# 08/03/2006
# - Adding BLAST type against repeat set as an option
#   This currently accepts blastn and tblastx
#
# 08/08/2006
# - Been working in HMM approach to MITE/MULE discovery
#   using HMMER as described for rice genome project
# - Changing Program Name to 'RepeatMiner'
#   -  becuase the program uses PANs, BLAST, and HMMER
# - Added RepHmmerRun subfunction
# - Added RepHmmerParse subfunction
#
# 08/09/2006
# - Added MULE models to RepHmmerRun and RepHmmerParse
# - Changed path of models to base of my user dir
# - Added FASTA File output for all derived contigs
# - Added a TPASE prediction from PFAM .. no hits 
#
# 08/10/2006
# - Reformat top to make it easier to do different
#   species. Cenchrus and Pennisetum are together in
#   a single program now. Just need to comment out
#   the vars that are not being used.
#
# 08/11/2006
# - Adding the sequence of the HMM repeat object to the
#   the HMMER parse output. This can not be done directly
#   from the HMMER output file in bioperl so I did this
#   using the start/end on the qry sequence.
# - Added length of the HMM hit and the BLAST hit to the 
#   METADATA output. This is the sum of the HSPs including
#   gaps in the lenght
#
# 08/15/2006
# - Changed the species name to a single variable that will
#   select among a larger group of variables.
#

#-----------------------------------------------------------+
# TO DO LIST
#-----------------------------------------------------------+
# - Working on Firefox mimetypes to open the 
#   Apollo and Ace files in the proper programs.
#   This may require writing shell scripts that
#   launch the proper programs with the needed
#   variables/parameters. May also want to do this
#   with some javascript that launches a shell script
# - Update database with the classification of the contigs?
