#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# RepMiner : jabablast.pl                                   |
# JAMIE'S ALL BY ALL BLAST to IDENTIFY REPEATS              |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 06/14/2006                                       |
# UPDATED: 06/27/2007                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  All by all blast analysis program. Creates matrix and    |
#  text files describing graph ready for analysis in        |
#  the Cytoscape graph visualization program.               |
#                                                           |
# DEPENDENCIES:                                             |
#  -BioPERL                                                 |
#  -NCBI BLAST                                              |
#  -GD                                                      |
#  -MySQL                                                   |
#                                                           |
# USAGE:                                                    |
#  jabablast.pl -h For full usage                           |
#                                                           |
#-----------------------------------------------------------+
# [ ] Add ability to use identity or bit score for the
#     parsing of the BLAST.
#     This will go to the function ParseTabBlast
# [ ] Consider the use of springgraph
#     http://www.chaosreigns.com/code/springgraph/
#


# Set the package name to RepMiner

package RepMiner::jabablast;

print "The RepMiner jabablast program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+

#use Graph::TransitiveClosure::Matrix;
use Graph;
use strict;                    # Gotta behave
use GD;                        # Draw using the GD program
use DBI();                     # Database interface
use Getopt::Std;               # Get options from the command line
use Bio::SearchIO;             # Parse BLAST output

#-----------------------------+
# GENERAL USE PROGRAM VARS    |
#-----------------------------+
my $XCrd;                      # Query coordinate (X coord)
my $YCrd;                      # Hit coordingate (Y coord)
my $QryName;                   # Name of the query

#-----------------------------+
# ALL-BY-ALL BLAST            |
#-----------------------------+
my $AllByAll;                  # Path to AllByAll blast output file
my $A_MinQryLen;               # Minimum query length for all by all blast
my $A_MinScore;                # Minimum Bit Score for all by all blast
my $A_MaxE;                    # Maximum E value for all by all blast
my $ResultCount = 0;           # Result count for BLAST
my $MinHitLen = "20";          # The minimum hit length to consider
my @BlastMatrix;               # Array to hold the blast results
my $BlastResult;               # The BLAST Result 
my $BlastHit;                  # Individual BLAST hit
my $HSP;                       # BLAST HSP Result
my $BlastDB;                   # The blast database that was queried
                               # This package level variable is used 
my $HitId;                     # Unique ID of the BLAST hit
my $SpHInfo;                   # Information related to the BLAST hit
my @HitDesc;                   # Initially split the hit description
my @HitInfo;                   # Array to hold the split of the hit info
my $HitLength;

#-----------------------------+
# CLASSIFICATION BLAST        |
#-----------------------------+
# Currently only a repeat based classification is allowed
my $RepBLAST;                  # Path to Blast against repeat DBs
my $MinQryLen;                 # Minimum query length for blast to repDb
my $MinScore;                  # Minimum Bit score for blast to repDb
my $MaxE;                      # Maximum E value for blast to repDb
my $RepBlastDb;                # Name of repeat database used to classify hits

#-----------------------------+
# GENERAL CLASSIFICATION VARS |
#-----------------------------+
# Variables for storing repeat element classification.
my $Class;
my $Subclass;
my $Superfamily;
my $Family;
my $Name;

#-----------------------------+
# MATRIX                      |
#-----------------------------+
my $xsc;                       # X base scale for drawing matrix (integer >1)
my $ysc;                       # Y base scale for drawing matrix (integer >1)
my $pxs;                       # Number of pixels to draw "dot" in matrix
my $MaxVal = 0;                # Max value to determine the coordinates
                               # of the hit matrix

#-----------------------------+
# GRAPH                       |
#-----------------------------+
my $NetDir;                    # Directory for network files
my $NetName;                   # Name of the *.sif file
my $GraphDir;                  # Directions of edges in graph
my $BlastFormat;               # Foramt of the blast file
                               # ie m = 0,8, or 9

#-----------------------------+
# DATABASE VARIABLES          |
#-----------------------------+
my $DbUserName;                # Database user name
my $DbUserPassword;            # Database user password
my $DbName;                    # Database name
my $tblAllByAll;               # Database table to store All by All Blast
my $tblQryCat;                 # Database table to store blast to repDB

#-----------------------------------------------------------+
# COMMAND LINE VARIABLES                                    |
#-----------------------------------------------------------+
my $Usage = "USAGE:\n".
    "jabablast.pl -i AllByAllBlast.blo -r RepElementBlast.blo\n". 
    "-o OutputDir -u UserName -d DatabaseName -Q\n\n".
    "jabablast.pl -C ConfigFile.jcfg\n\n".
    "+-----------------------------------------------------------+\n".
    "| REQUIRED ARGUMENTS                                        |\n".
    "+-----------------------------------------------------------+\n".
    " -i The AllByAllBlast File to parse [String]\n".
    " -o The path of the output directory [String]\n".
    " -u The user name for connection to the database [String]\n".
    " -b The database name to use for the database connection. [String]\n".
    "    This database MUST already exist in you MySQL database.\n".
    "  OR\n".
    " -Z Config file that provides the above information. [String]\n".
    "\n".
    "jabablast.pl -h TO GET FULL ARGUMENTS LIST.\n";

my %Options;
# leave g for GraphViz output path
# j
# p
# t
# v
# w
# z

getopts('a:b:c:d:e:f:hi:k:l:m:n:o:r:s:u:x:y:CE:GL:N:S:MQZ:', \%Options);

#-----------------------------+
# THE FOLLOWING OPTIONS ARE   |
# REQUIRED                    |
#-----------------------------+
# TODO:
# This should be changed to allow for the option of visualizing the
# AllByAll BLAST results without the repeat data categorization path.
# THe output dir should include the / at the end of the DIR PAth
my $PrintHelp = $Options{h};
my $Config = $Options{Z};
if ($Config){
    print "\n\nThe config file is:\n\t$Config\n\n";}
my $quiet = $Options{Q};
my $CreateMatrix = $Options{M};
my $LaunchCytoscape = $Options{C};
my $RunGraphviz = $Options{G};

#-----------------------------+
# PRINT HELP IF REQUESTED     |
#-----------------------------+
if ($PrintHelp)
{
    &PrintHelp;
    exit;
}

#-----------------------------+
# SET VARIABLES               | 
#-----------------------------+
# If a config file path was given set the user variables
# using the config file, otherwise get all options
# from the command line.
if ($Config)
{

    # DATABASE VARIABLES
    $DbUserName = &ParseConfigFile($Config,"DbUserName");
    $DbName = &ParseConfigFile($Config,"DbName");
    $NetDir = &ParseConfigFile($Config,"NetDir");
    
    # INPUT/OUTPUT FILES
    $AllByAll = &ParseConfigFile($Config,"BLAST_AllByAll");
    $RepBLAST = &ParseConfigFile($Config,"BLAST_RepDB");
    $NetName = &ParseConfigFile($Config,"NetName") ||
	"Network"; # Added a default name 
    
    # BLAST RELATED VARIABLES
    $A_MinQryLen = &ParseConfigFile($Config,"A_MinQryLen")
	|| "50";
    $A_MinScore = &ParseConfigFile($Config,"A_MinScore")
	|| "150";
    $A_MaxE = &ParseConfigFile($Config,"A_MaxE")
	|| "1.0e-05";
    
    # CATEGORIZATION BLAST VARS
    $MinQryLen = &ParseConfigFile($Config,"MinQryLen")
	|| "20";
    $MinScore = &ParseConfigFile($Config,"MinScore")
	|| "50";
    $MaxE = &ParseConfigFile($Config,"MaxE")
	|| "1.0e-03";
    
    # ALL BY ALL MATRIX DRAWING VARIABLES
    $xsc = &ParseConfigFile($Config,"xsc")
	|| "2";  # The X coordinate scaling factor
    $ysc = &ParseConfigFile($Config,"ysc") 
	|| "2";  # The Y coordinate scaling factor
    $pxs = &ParseConfigFile($Config,"pxs") 
	|| "4";  # Pixel size of the matched dots

    # DATABASE TABLES
    $tblAllByAll = &ParseConfigFile($Config,"ABATable") ||
	"tblAllByAll";
    $tblQryCat = &ParseConfigFile($Config,"RepeatTable") ||
	"tblRepeatID";
    
}else{

    $DbUserName = $Options{u} ||
	die "You must provide a database user name\n$Usage\n";
    $DbName = $Options{b} ||
	die "You must provide a database name to connect to";
    $AllByAll = $Options{i} || 
	die "You must provide an AllByAll BLAST output path\n$Usage\n";
    $NetDir = $Options{o} ||
	die "You must provide an output Directory path\n$Usage\n";
    $NetName = $Options{f} || "Network"; # Default network name is Network

    #-----------------------------+
    # GRAPH OPTIONS               |
    #-----------------------------+
    $GraphDir = $Options{d} || "0";

    #-----------------------------+
    # BLAST OPTIONS               |
    #-----------------------------+
    $BlastFormat = $Options{m} || "8";
    # ALL BY ALL BLAST VARS
    $A_MinQryLen = $Options{l} || "50";
    $A_MinScore = $Options{s} || "150";
    $A_MaxE = $Options{e} || "1.0e-05";

    #-----------------------------+
    # CATEGORIZATION BLAST VARS   |
    #-----------------------------+
    # Currently limited to repeat database blast, other
    # node categorizatioin algorithms should be considered
    $RepBLAST = $Options{r};
    $MinQryLen = $Options {L} || "50";
    $MinScore = $Options{S} || "50";
    $MaxE = $Options{E} || "1.0e-03";
    $RepBlastDb = $Options{N};      # Name of repeat DB
                                    # only needed for tab blast
    
    #-----------------------------+
    # DATABASE OPTIONS            |
    #-----------------------------+
    # Table with All by All BLAST results        
    $tblAllByAll = $Options{a} || "tblAllByAll";      
    # Table with repeat ID BLAST results
    $tblQryCat = $Options{c} || "tblRepeatID";

    #-----------------------------+
    # MATRIX OPTIONS              |
    #-----------------------------+
    $xsc = $Options{x} || "2";    # The X coordinate scaling factor
    $ysc = $Options{y} || "2";    # The Y coordinate scaling factor
    $pxs = $Options{p} || "4";    # Pixel size of the matched dots
    
    
} # END OF GETTING INFO FROM COMMAND LINE

#-----------------------------------------------------------+
# CHECK USER VARIABLES BEFORE CONTINUING WITH THE PROGRAM   |
#-----------------------------------------------------------+
#
# TODO: Turn this outupt and question off with the -q flag
#       I can currently use an incorrected password to quit the program.
print "NETDIR:\n\t$NetDir\n";
print "NETNAME:\n\t$NetName\n";

print "ALLBYALL:\n\t$AllByAll\n";
print "A_MinQryLen:\n\t$A_MinQryLen\n";
print "A_MinScore:\n\t$A_MinScore\n";
print "A_MaxE:\n\t$A_MaxE\n";

print "DBUSER NAME:\n\t$DbUserName\n";
print "DBNAME:\n\t$DbName\n";
print "tblAllByAll:\n\t$tblAllByAll\n";

if ($RepBLAST){
    print "MinQryLen:\n\t$MinQryLen\n";
    print "MinScore:\n\t$MinScore\n";
    print "MaxE:\n\t$MaxE\n";
    print "REPBLAST:\n\t$RepBLAST\n";
    print "tblQryCat:\n\t$tblQryCat\n";
    print "DB_NAME:\n\t$RepBlastDb\n"
#    print "dbName:
}

if ($CreateMatrix)
{
    print "XSC:\n\t$xsc\n";
    print "YSC:\n\t$ysc\n";
    print "PXS:\n\t$pxs\n";
}

#-----------------------------+
# GET USER PASSWORD           |
#-----------------------------+
print "\nPassword for $DbUserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
$DbUserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $DbUserPassword;

#-----------------------------------------------------------+
# DB I/O CONNECT                                            | 
#-----------------------------------------------------------+
my $dbh = DBI->connect("DBI:mysql:database=$DbName;host=localhost",
		       $DbUserName, $DbUserPassword,
		       {'RaiseError' => 1});

my $RepDbName = "dbRep";
my $RepDB = DBI->connect("DBI:mysql:database=$RepDbName;host=localhost",
			   $DbUserName, $DbUserPassword,
			   {'RaiseError' => 1});

#-----------------------------------------------------------+
# OUTPUT FILES                                              | 
#-----------------------------------------------------------+

# CREATE NETWORK DIR IF NEEDED
# See if it is possible to make the parent dirs
# if required. Will try mkdir -p otherwise
# do a system call to mkdir -p
mkdir $NetDir, 0777 unless (-e $NetDir);

#-----------------------------+
# OUTPUT FILES                |
#-----------------------------+
my $SumOut = $NetDir."SumInfo.txt";
my $SifOut = $NetDir.$NetName.".sif";
my $NA_BACOut = $NetDir."BAC.NA";
my $NA_RepClass = $NetDir."RepClass.NA";
my $NA_RepName = $NetDir."RepName.NA";
my $NA_SeqData = $NetDir."Seq.NA";
my $EA_BitScore = $NetDir."BitScore.EA";
my $EA_Sig = $NetDir."HitSign.EA";
my $EA_PID = $NetDir."PID.EA";
my $GraphOut = $NetDir."ColorMatrix.png";
my $HSPSif = $NetDir."HSP.sif";
my $HSPFI = $NetDir."HSPFI.EA"; 
my $HSPLen = $NetDir."HSPLen.EA";
my $DrosOut = $NetDir."DrosGPI.NA";

# SUMMARY INFORMATION
open (SUMOUT, ">$SumOut") ||
    die "Can not open $SumOut\n";     

print SUMOUT "+----------------------------------------------------------+\n";
print SUMOUT "| RepMiner Summary Outupt File                             |\n";
print SUMOUT "+----------------------------------------------------------+\n";


# NETWORK SIF FILE
# For BLAST HITS
open (SIFOUT, ">$SifOut") ||   #Network *.SIF file
    die "Can not open $SifOut\n"; 

# NETWORK SIF FILE FOR
# BLAST HSP
open (HSPOUT, ">$HSPSif") ||
    die "Can not open $HSPSif\n";

# NODE ATTRIBUTE : BAC
open (BACOUT, ">$NA_BACOut") ||
    die "Can not open $NA_BACOut\n";  
print BACOUT "BAC\n";

# NODE ATTRIBUTE : REPEAT CLASS
open (REPOUT, ">$NA_RepClass") ||
    die "Can not open $NA_RepClass\n" ;
print REPOUT "RepeatClass\n";

open (REPNAME, ">$NA_RepName") ||
    die "Can not open $NA_RepName\n";
print REPNAME "RepeatName\n";

# NODE ATTRIBUTE : SEQUENCE
open (SEQOUT, ">$NA_SeqData") ||
    die "Can not open $NA_SeqData\n";
print SEQOUT "SeqData\n";

# EDGE ATTRIBUTE : BIT SCORE
open (BITOUT, ">$EA_BitScore") ||
    die "Can not open $EA_BitScore\n";
print BITOUT "BitScore\n";

# EDGE ATTRIBUTE : SIGNIFICANCE
open (SIGOUT, ">$EA_Sig") ||
    die "Can not open $EA_Sig";
print SIGOUT "HitSignificance\n";

# EDGE ATTRIBUTE : PERCENT IDENTITY
open (PIDOUT, ">$EA_PID") ||
    die "Can not open $EA_PID";
print PIDOUT "PercentIdentity\n";

# HSP EDGE ATTRIBUTE : FRACTION IDENTICAL
open (HSPFI, ">$HSPFI") ||
    die "Can not open $HSPFI\n";
print HSPFI "HSPFracIdent\n";

# HSP EDGE ATTRBIUTE : TOTAL LENGTH
open (HSPLEN, ">$HSPLen") ||
    die "Can not open $HSPLen\n";
print HSPLEN "HSPLength\n";

#-----------------------------+
# CHECK DATABASE AND CREATE   |
# TABLES IF NEEDED            | 
#-----------------------------+
print "Database setup.\n";
&DbSetup;

#-----------------------------+
# PARSE THE ALL BY ALL BLAST  |
# AND LOAD TO THE BLAST MATRIX|
#-----------------------------+
#&LoadAllByAll;
if ($BlastFormat == '8') {

    # Parse AxA to SIF Files
    &ParseTabBLAST ($AllByAll);

    # Parse AxA to Graph Object and Classify
    &ParseTabBLAST2Graph ($AllByAll);

}
# AT WORK HERE THIS WILL LOAD TAB DELIMITED ALL BY ALL BLAST

#-----------------------------+
# PARSE THE REPEAT BLAST TO   |
# FETCH THE REPEAT CATEGORIES |
# OF THE DB                   |
#-----------------------------+
# If a repeat blast database is not provided it will
# not run the load rep class subfunction.
# If a RepBlast path was passed at the command line 
if ($RepBLAST) {
    #&LoadRepClass;
    # For tab delimited BLAST of repeat database

    #
    #&LoadTabRepClass($RepBlastDb);

    #                               $RepBlastDb
    &LoadTabRepClassNew( $RepBLAST, $RepBlastDb);

}

#-----------------------------+
# PARSE THE BLASTX AGAINST    |
# THE DROSOPHILA              |
# GAG/POL/INTEGRASE PROTEINS  |
#-----------------------------+
#print "ATTEMPTING TO LOAD DROS FILE\n";
#&LoadDrosGPI ($DrosBLAST,$DrosOut);

#-----------------------------+
# DO QUERY TO MERGE DATA FROM |
# QRY INFORMATION FROM THE    |
# REPEAT DB BLAST TO THE ALL  |
# BY ALL INFORMATION TABLE    |
#-----------------------------+
if ($RepBLAST) {
    &UpdateRepCat;
}

#-----------------------------+
# LAUNCH CYTOSCAPE            |
# NETWORK VIEWEING PROGRAM    |
# WITH THE SIF FILE AND       |
# THE EDGE ATTRIBUTES AND     |
# NODE ATTRIBUTES THAT WERE   |
# CREATED                     |
#-----------------------------+
if ($LaunchCytoscape)
{
#    &LaunchCytoscapeOld ($SifOut);
    &LaunchCytoscapeNew ($SifOut);
#    &LaunchCytoscape_2_4 ($SifOut);
}

#-----------------------------+
# DRAW COLOR VALUED XY PLOT   |
# OF THE SPARSE MATRIX        |
# REPRESENTING THE RESULTS OF |
# THE ALL BY ALL BLAST        |
#-----------------------------+
if ($CreateMatrix)
{
    &DrawXYPlot;
}

#-----------------------------------------------------------+
# CLOSE TEXT OUTPUT FILES                                   |
#-----------------------------------------------------------+
close SUMOUT;                  # Close the summary output file
close SIFOUT;                  # Close the sif network output file
close BACOUT;                  # Close BAC NA file 
close REPOUT;                  # Close Repeat Class NA file
close BITOUT;                  # Close the BitScore EA file
close REPNAME;                 # Close the RepeatName NA file
close SIGOUT;                  # Close the Hit significance EA file
close HSPOUT;                  # Close the HSP SIF output file
close HSPFI;                   # Close the HSP fraction identical EA file
close HSPLEN;                  # Close the HSP Length EA file
close PIDOUT;                  # Close the PercentIdentity EA file

print "The jabablast program has finished.\n";
exit;

#-----------------------------------------------------------+
#                                                           |
# SUBFUNCTIONS                                              |
#                                                           |
#-----------------------------------------------------------+

sub ParseConfigFile
{
#-----------------------------------------------------------+
# This will parase a configuration text file to set         |
# variables that would normally be set at the command line  |
# This is not the fastest way to parse the config file      |
# but this makes the subfuction reusable.                   |
# If the variable name occurs more then once in the text    |
# file, the last occurrence will be used.                   |
#-----------------------------------------------------------+
    my $ConfigFilePath = $_[0];
    my $VarName = $_[1];
    my $VarValue;
    
    open (CONFILE, $ConfigFilePath) ||
	die "Could not open config file:\n\t$ConfigFilePath";
    
    while (<CONFILE>)
    {
	chomp;                 # Remove newline character
	unless (m/\#.*/)       # Ignore comment lines
	{
	    my @SplitLine = split;
	    if ($SplitLine[0] =~ $VarName){
		$VarValue = $SplitLine[1];}
	}
    }
    close CONFILE;
    return $VarValue;

}

sub GetMaxVal
{
#-----------------------------+
# SELECTS THE MAXIMUM VALUE   |
# FROM A COLUMN AND RETURNS   |
# THE RESULT AS A STRING      |
#-----------------------------+

    my ($col, $table) = @_;
    my ($cur,$result);
    my $SelectSQL = "SELECT MAX($col) FROM $table";
    $dbh->do($SelectSQL);
    $cur = $dbh->prepare($SelectSQL);
    $cur->execute();
    my @row=$cur->fetchrow;
    $result=$row[0];
    $cur->finish();
    return $result;
}

sub UpdateRepCat
{
    
    print "Updating categories in ".$tblAllByAll."\n";
    
    my $UpdateCats = "UPDATE ".$tblAllByAll.",".$tblQryCat.
	" SET ".$tblAllByAll.".qry_cat=".$tblQryCat.".qry_cat".
	" WHERE ".$tblAllByAll.".qry_id = ".$tblQryCat.".qry_id".
	" AND ".$tblAllByAll.".qry_cat IS NULL";
    
    $dbh->do($UpdateCats);

}


sub DbSetup
{
#-----------------------------+
# SET UP THE DATABASE BY      |
# MAKING TABLES IF NEEDED     |
# AND CHECKING IF THE USER    |
# WANTS TO DELETE OLDER TABLES|
#-----------------------------+


    #-----------------------------+
    # ALL BY ALL TABLE            |
    #-----------------------------+
    if (&does_table_exist($tblAllByAll))
    {
	print "The table $tblAllByAll already exits.\n";
	my $question = "Do you want to overwrite the existing table.";
	my $answer = &UserFeedback($question);
	if ($answer =~ "n"){exit;}
	# Could add code to use a different name for the table
	$dbh->do("DROP TABLE ".$tblAllByAll);
    }

    my $CreateAllByAllTbl = "CREATE TABLE .".$tblAllByAll.
	" (".
	" rownum INT(10) NOT NULL AUTO_INCREMENT,".
	" qry_num INT(10),".
	" hit_num INT(10),".
	" qry_id VARCHAR(255),".
	" qry_cat VARCHAR(50),".
	" hit_cat VARCHAR(50),". # Category of the hit to id selfhits
	" num_hits INT(10),".    # The number of hits in the all by all blast
	" KEY (rownum))";

    $dbh->do($CreateAllByAllTbl);


    # Add an index to the qry_id field to speed up queries
    my $IndexAllByAllTbl = "ALTER TABLE ".
	$tblAllByAll." ADD INDEX ( `qry_id` )";
    $dbh->do($IndexAllByAllTbl);

    # Only create query category table if RepBlast is true
    if ($RepBLAST)
    {
	#-----------------------------+
	# QUERY CATEGORY TABLE        |
	#-----------------------------+
	if (&does_table_exist($tblQryCat))
	{
	    print "The table $tblQryCat already exits.\n";
	    my $question = "Do you want to overwrite the existing table?";
	    my $answer = &UserFeedback($question);
	    if ($answer =~ "n"){exit;}
	    # Could add code to use a different name for the table
	    $dbh->do("DROP TABLE ".$tblQryCat);
	}
	
	my $CreateQryCatTbl = "CREATE TABLE ".$tblQryCat.
	    "(".
	    " rownum INT(10) NOT NULL AUTO_INCREMENT,".
	    " qry_id VARCHAR(255),".
	    " qry_cat VARCHAR(50),".
	    " KEY (rownum))";
	
	$dbh->do($CreateQryCatTbl);
	
	# Add an index to the qry_id field to speed up queries
	my $IndexQryCatTbl = "ALTER TABLE ".$tblQryCat.
	    " ADD INDEX ( `qry_id` ) ";
	$dbh->do($IndexQryCatTbl);
    }
}



sub LoadDrosGPI
{
#-----------------------------+
# LOAD BLASTX OUTPUT FROM THE |
# DROSOPHILA GAG/POL/INTEGRASE|
# PROTEIN DOMAINS             |
#-----------------------------+

    my $DrosInFile = $_[0];
    my $DrosOutFile = $_[1];
    my $QryName;
    my $BlastDB;
    my $DHitName;              # Name of Dros Gag/Pol Hit
    my $DrosHitCount = 0;

    open (DROSOUT, ">$DrosOutFile") ||
	die "Can not open output file: $DrosOutFile\n";
    print DROSOUT "DrosGPI\n";

    my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
					  '-file'   => $DrosInFile)
	||
        die "Could not open BLAST input file:\n$DrosInFile.\n";

    while (my $BlastResult = $BlastReport->next_result())
    {
	#print "Searching\n";
	while (my $BlastHit = $BlastResult->next_hit())
	{

	    $QryName = $BlastResult->query_name;
	    $BlastDB = $BlastResult->database_name;
	    my @tmpqry = split(/\|/, $QryName );
	    my $OutName = $tmpqry[0];

	    $DrosHitCount++;
	    print $DrosHitCount."\n";
	    my $DHitName= $BlastHit->name();

	    print DROSOUT $OutName."=".$DHitName."\n";

	    while (my $BlastHSP = $BlastHit->next_hsp())
	    {
		print "\t".$BlastHSP->query_string()."\n";
		print "\t".$BlastHSP->homology_string()."\n";
		print "\t".$BlastHSP->hit_string()."\n";
		print "\n";
	    }
	    
	} # End of while BlastHit next_hit
	
    } # End of while BlastResult next_result

    close DROSOUT;
}

sub LoadRepClass
{

    print "Defining repeat classes.\n";

    #-----------------------------+
    # LOAD THE REPEAT             |
    # CLASSIFICATION BLAST OUTPUT |
    #-----------------------------+
    my $LoadRec;
    my $QryName;
    my $RepClass;

    my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
					  '-file'   => $RepBLAST,
					  '-signif' => $MaxE,
					  '-min_query_len' => $MinQryLen,
					  '-score' => $MinScore ) ||
    die "Could not open BLAST input file:\n$RepBLAST.\n";

    while ($BlastResult = $BlastReport->next_result())
    {
	$QryName = $BlastResult->query_name;
	$BlastDB = $BlastResult->database_name;

	#print "Rep Search".$QryName.":".$BlastDB."\n";

	my $NumHits = $BlastResult->num_hits;
	my $HitCount = "0";

	while ( $BlastHit = $BlastResult->next_hit())
	{

	    #-----------------------------+
	    # TREP REPEAT DATABASE        |
	    #-----------------------------+
	    if ($BlastDB =~ 'TREP_8' ||
		$BlastDB =~ 'TREP_9')
	    {
		$HitCount++;

		my $HitId = $BlastHit->accession()  || "UnkAcc";

		$RepClass = &GetTREPClass($HitId);
		$Name = &GetTREPName($HitId);

		# Don't show output on quiet runs
		if (! $quiet)
		{
		    print $QryName."\n";
		    print "\t".$HitId.":".$Name."\n";
		    print "\t".$RepClass."\n";
		}

	    }

	    #-----------------------------+
	    # REPBASE PLANTS              |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'RB_pln')
	    {

		$HitCount++;
		$Class = "UNK";
		$Subclass = "UNK";
		$Superfamily = "UNK";
		$Name = $BlastHit->name();
		
		my @SpName = (split /\#/, $Name);
		my $CatSearch = trim($SpName[1] || "NONE");
		
		$RepClass = &GetRBClass($CatSearch); 

	    }

	    #-----------------------------+
	    # MIPS DATABASE               |
	    # This follows the RepBASE    |
	    # foramt.                     |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'mips_REdat_4_3')
	    {

		$HitCount++;
		$Class = "UNK";
		$Subclass = "UNK";
		$Superfamily = "UNK";
		$Name = $BlastHit->name();
		
		my @SpName = (split /\#/, $Name);
		my $CatSearch = trim($SpName[1] || "NONE");
		
		$RepClass = &GetRBClass($CatSearch); 

		# Show what class is currently being read
		if (! $quiet)
		{
		    print $Name."\n";
		}
		
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
	    # TIGR BRASSICACEAE REPEATS   |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRbras')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name);
	    }

	    #-----------------------------+
	    # TIGR FABACEAE REPEATS       |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRfab')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name);
	    }

	    #-----------------------------+
	    # TIGR SOLANACEAE REPEATS     |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRsol')
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
		
		my $Desc = $BlastHit->description();
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
		
	    } else {
		#-----------------------------+
		# DATABASE NOT RECOGNIZED     |
		#-----------------------------+
		print "ERROR. The repeat database $BlastDB is not".
		    " recognized by RepMiner.\n";
	    }


	    # CURRENTLY WILL ONLY LOAD THE BEST HIT INTO
	    # THE DATABASE FOR A NICE QUICK AND DIRTY
	    # ATTEMPT TO GET THIS TO WORK
	    if ($HitCount == "1")
	    {


		# Can do the HSP parse here if only wanted for the 
		# first hit.

		# Check to see if a classification already exists for
		# for this record. As it stands now, the later qry_cats
		# will overwrite the newer qry cats.

		# Can change the load record only record the classification

		$LoadRec = "INSERT INTO ".$tblQryCat.
		    " (qry_id, qry_cat) ".
		    " VALUES (".
		    " '".$QryName."',".
		    " '".$RepClass."'".
		    ")";
		$dbh->do($LoadRec);

		my @tmpqry = split(/\|/, $QryName );
		my $OutName = $tmpqry[0]; 

		# TRY TO TRACK DOWN AN ERROR
		
		print "SQL\t\t$LoadRec\n";
		print "\t\t".$OutName."\n";
		print "\t\t".$RepClass."\n";
		print "\t\t".$Name."\n";

		# It appears that cycoscape shows the last one that 
		# was added to the list, and ignores earlier 
		# classifications
		print REPOUT $OutName."=".$RepClass."\n";
		print REPNAME $OutName."=".$Name."\n";

	    }

	} # End of while BlastResult-next_hit

    } # End of while BlastReport next_result

}

# Begin to work with tab delimited -m8 or -m9 format output

sub LoadTabRepClassNew
{
    # Only the first hit from the BLAST is used here
    # Since the variables are not the same as the 
    # SeqIO objects, this code will need to be modified
    # for the different databases.

    # Parsing the BLAST file using native PERL functions for
    # working with text files is much faster then using the
    # BioPERL SeqIO
    my $TabBlastFile = $_[0];
    my $BlastDB = $_[1];
    my $PreQry = "";              # Qry ID of previous hit 
    my $CurQry = "";              # Qry ID of current hit
    my $other; # Var to extra information with the name
    my $RepClass;
    my $Name;
    my $HitCount = 0;

    open (BLASTIN, "<".$TabBlastFile) ||
	die "Can not open BLAST file.$TabBlastFile.\n";

    while (<BLASTIN>)
    {
	chomp;                 # Remove newline character
	unless (m/^\#.*/)       # Ignore comment lines, works with -m 9 output
	{
	    
	    my ($QryId, $SubId, $PID, $Len, 
		$MisMatch, $GapOpen, 
		$QStart,$QEnd, $SStart, $SEnd,
		$EVal, $BitScore) = split(/\t/);

	    # Trim leading white space from Bit score
	    $BitScore =~ s/^\s*(.*?)\s*$/$1/;

	    $CurQry = $QryId;

	    # Only use the top hit 
	    unless ($CurQry =~ $PreQry)
	    {
		#-----------------------------+
		# PARSE INFORMATION DEPENDENT |
		# ON THE REPEAT DATABASE USED | 
		#-----------------------------+ 
		if ($BlastDB =~ 'TREP_8' ||
		    $BlastDB =~ 'TREP_9')
		{
		    $HitCount++;
		    
		    my $HitId = $BlastHit->accession()  || "UnkAcc";
		    
		    $RepClass = &GetTREPClass($HitId);
		    $Name = &GetTREPName($HitId);
		    
		    # Don't show output on quiet runs
		    if (! $quiet)
		    {
			print $QryName."\n";
			print "\t".$HitId.":".$Name."\n";
			print "\t".$RepClass."\n";
		    }
		    
		}
		
		#-----------------------------+
		# REPBASE PLANTS              |
		#-----------------------------+
		elsif ($BlastDB =~ 'RB_pln')
		{
		    
		    $HitCount++;
		    $Class = "UNK";
		    $Subclass = "UNK";
		    $Superfamily = "UNK";
		    $Name = $BlastHit->name();
		    
		    my @SpName = (split /\#/, $Name);
		    my $CatSearch = trim($SpName[1] || "NONE");
		    
		    $RepClass = &GetRBClass($CatSearch); 
		    
		}
		
		#-----------------------------+
		# MIPS DATABASE               |
		# This follows the RepBASE    |
		# foramt.                     |
		#-----------------------------+
		elsif ($BlastDB =~ 'mips_REdat_4_3' || 
		       $BlastDB =~ 'mips')
		{
		    
		    $HitCount++;
		    $Class = "UNK";
		    $Subclass = "UNK";
		    $Superfamily = "UNK";
		    $Name = $BlastHit->name();
		    
		    my @SpName = (split /\#/, $Name);
		    my $CatSearch = trim($SpName[1] || "NONE");
		
		    $RepClass = &GetRBClass($CatSearch); 
		    
		    # Show what class is currently being read
		    if (! $quiet)
		    {
			print "MIPS Name:\t".$Name."\n";
		    }
		    
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
		# TIGR BRASSICACEAE REPEATS   |
		#-----------------------------+
		elsif ($BlastDB =~ 'TIGRbras')
		{
		    $HitCount++;
		    $Name = $BlastHit->name();
		    $RepClass = &GetTIGRClass($Name);
		}
		
		#-----------------------------+
		# TIGR FABACEAE REPEATS       |
		#-----------------------------+
		elsif ($BlastDB =~ 'TIGRfab')
		{
		    $HitCount++;
		    $Name = $BlastHit->name();
		    $RepClass = &GetTIGRClass($Name);
		}
		
		#-----------------------------+
		# TIGR SOLANACEAE REPEATS     |
		#-----------------------------+
		elsif ($BlastDB =~ 'TIGRsol')
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
		    
		    my $Desc = $BlastHit->description();
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
		    } # End of not quiet
		}
		
		#-----------------------------+
		# SAN MIGUEL REPEAT DATABASE  |
		#-----------------------------+
		elsif ($BlastDB =~ 'SanMiguel')
		{
		    $Name = $SubId;
		    ($RepClass,$other) = split(/_/, $Name);
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
		    
		} else {
		    #-----------------------------+
		    # DATABASE NOT RECOGNIZED     |
		    #-----------------------------+
		    #print "ERROR. The repeat database $BlastDB is not".
		    #	" recognized by RepMiner.\n";
		
		    $Name = $SubId;
		    $RepClass = "UNK";
		}

		#-----------------------------+
		# San Miguel Database         |
		#-----------------------------+
#		$Name = $SubId;
#		($RepClass,$other) = split(/_/, $Name);
		

		#-----------------------------+
		# WRITE OUTPUT                |
		#-----------------------------+

		my $LoadRec = "INSERT INTO ".$tblQryCat.
		    " (qry_id, qry_cat) ".
		    " VALUES (".
		    " '".$QryId."',".
		    " '".$RepClass."'".
		    ")";
		$dbh->do($LoadRec);
		
		my @tmpqry = split(/\|/, $QryId );
		my $OutName = $tmpqry[0]; 
		
		# TRY TO TRACK DOWN AN ERROR
		#print "SQL\t\t$LoadRec\n";
		print "\t\t".$OutName."\n";
		#print "\t\t".$RepClass."\n";
		#print "\t\t".$Name."\n";
		# It appears that cycoscape shows the last one that 
		# was added to the list, and ignores earlier 
		# classifications
		print REPOUT $OutName."=".$RepClass."\n";
		print REPNAME $OutName."=".$Name."\n";
		
	    }

	    $PreQry = $CurQry;	    
	    
	} # End of remove newline character
	
    } # End of while BLASTN

}

sub LoadTabRepClass
{
    my $BlastDB = $_[0];

    print "Defining repeat classes using tab delim blast report.\n";
    print "The repeat blast database is:".$BlastDB."\n";

    #-----------------------------+
    # LOAD THE REPEAT             |
    # CLASSIFICATION BLAST OUTPUT |
    #-----------------------------+
    my $LoadRec;
    my $QryName;
    my $RepClass;

    my $BlastReport = new Bio::SearchIO ( '-format' => 'blasttable',
					  '-file'   => $RepBLAST,
					  '-signif' => $MaxE,
					  '-min_query_len' => $MinQryLen,
					  '-score' => $MinScore ) ||
    die "Could not open BLAST input file:\n$RepBLAST.\n";

    #$BlastDB = $BlastResult->database_name;

    while ($BlastResult = $BlastReport->next_result())
    {
	$QryName = $BlastResult->query_name;

	#print "Rep Search".$QryName.":".$BlastDB."\n";

	my $NumHits = $BlastResult->num_hits;
	my $HitCount = "0";

	while ( $BlastHit = $BlastResult->next_hit())
	{

	    #-----------------------------+
	    # TREP REPEAT DATABASE        |
	    #-----------------------------+
	    if ($BlastDB =~ 'TREP_8' ||
		$BlastDB =~ 'TREP_9')
	    {
		$HitCount++;

		my $HitId = $BlastHit->accession()  || "UnkAcc";

		$RepClass = &GetTREPClass($HitId);
		$Name = &GetTREPName($HitId);

		# Don't show output on quiet runs
		if (! $quiet)
		{
		    print $QryName."\n";
		    print "\t".$HitId.":".$Name."\n";
		    print "\t".$RepClass."\n";
		}

	    }

	    #-----------------------------+
	    # REPBASE PLANTS              |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'RB_pln')
	    {

		$HitCount++;
		$Class = "UNK";
		$Subclass = "UNK";
		$Superfamily = "UNK";
		$Name = $BlastHit->name();
		
		my @SpName = (split /\#/, $Name);
		my $CatSearch = trim($SpName[1] || "NONE");
		
		$RepClass = &GetRBClass($CatSearch); 

	    }

	    #-----------------------------+
	    # MIPS DATABASE               |
	    # This follows the RepBASE    |
	    # foramt.                     |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'mips_REdat_4_3' || 
		   $BlastDB =~ 'mips')
	    {

		$HitCount++;
		$Class = "UNK";
		$Subclass = "UNK";
		$Superfamily = "UNK";
		$Name = $BlastHit->name();
		
		my @SpName = (split /\#/, $Name);
		my $CatSearch = trim($SpName[1] || "NONE");
		
		$RepClass = &GetRBClass($CatSearch); 

		# Show what class is currently being read
		if (! $quiet)
		{
		    print "MIPS Name:\t".$Name."\n";
		}
		
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
	    # TIGR BRASSICACEAE REPEATS   |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRbras')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name);
	    }

	    #-----------------------------+
	    # TIGR FABACEAE REPEATS       |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRfab')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		$RepClass = &GetTIGRClass($Name);
	    }

	    #-----------------------------+
	    # TIGR SOLANACEAE REPEATS     |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'TIGRsol')
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
		
		my $Desc = $BlastHit->description();
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
		} # End of not quiet
	    }

	    #-----------------------------+
	    # SAN MIGUEL REPEAT DATABASE  |
	    #-----------------------------+
	    elsif ($BlastDB =~ 'SanMiguel')
	    {
		$HitCount++;
		$Name = $BlastHit->name();
		# Add code to simply name here to the string
		# to the left of the underscore _
		#$Name 
		# Using the short name 'family' as
		# the current RepClass although that
		# is not really correct

		# The following did not split
#		my @SplitName = split(/_/, $Name ) ||
#		    die "Could not split Hit Name: $Name";
#		my $SplitLen = @SplitName;
#		print "LEN:".$SplitLen."\n";
#		$RepClass = $SplitName[0];

		my $other;
		($RepClass,$other) = split(/_/, $Name);


#	    my @tmphitid = split(/\|/, $SubId ) ||
#		die "Could not split SbjID";
#	    $YCrd = $tmphitid[0];


		#f (! $quiet)
		#
		#   #print $QryName."\n";
		#   print "SanMig Name:\t".$Name."\n";
		#   print "SanMig Class:\t".$RepClass."\n";
		#

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
		
	    } else {
		#-----------------------------+
		# DATABASE NOT RECOGNIZED     |
		#-----------------------------+
		print "ERROR. The repeat database $BlastDB is not".
		    " recognized by RepMiner.\n";
	    }


	    # CURRENTLY WILL ONLY LOAD THE BEST HIT INTO
	    # THE DATABASE FOR A NICE QUICK AND DIRTY
	    # ATTEMPT TO GET THIS TO WORK
	    if ($HitCount == "1")
	    {

		# Can do the HSP parse here if only wanted for the 
		# first hit.

		# Check to see if a classification already exists for
		# for this record. As it stands now, the later qry_cats
		# will overwrite the newer qry cats.

		# Can change the load record only record the classification

		$LoadRec = "INSERT INTO ".$tblQryCat.
		    " (qry_id, qry_cat) ".
		    " VALUES (".
		    " '".$QryName."',".
		    " '".$RepClass."'".
		    ")";
		$dbh->do($LoadRec);

		my @tmpqry = split(/\|/, $QryName );
		my $OutName = $tmpqry[0]; 

		# TRY TO TRACK DOWN AN ERROR
		print "SQL\t\t$LoadRec\n";
		print "\t\t".$OutName."\n";
		print "\t\t".$RepClass."\n";
		print "\t\t".$Name."\n";
		# It appears that cycoscape shows the last one that 
		# was added to the list, and ignores earlier 
		# classifications
		print REPOUT $OutName."=".$RepClass."\n";
		print REPNAME $OutName."=".$Name."\n";

	    } # End of if HitCount ==1

	} # End of while BlastResult-next_hit

    } # End of while BlastReport next_result

}

sub ParseTabBLAST
{
#-----------------------------+
# LOAD A TAB DELIMITED ALL BY |
# ALL BLAST REPORT            |
#-----------------------------+
# This is the tab delimited text file produced from 
# the m8 or m9 option in BLAST, use of the m9 option
# will required that lines starting with # must
# be ignored, I therfore have the regexp (m/^\#.*/).
# The available information:
#  - $QryId ----> Id of the query sequence 
#  - $SubId ----> Id of the subect seqence (from the database
#  - $PID ------> Percent ID
#  - $Len ------> Alignment length
#  - $MisMatch -> Number of mismatches
#  - $GapOpen --> Number of gap openings
#  - $QStart ---> Start of alignment on query sequence
#  - $QEnd -----> End of alignment on query sequence
#  - $SStart ---> Start of aligment on Subject
#  - $SEnd -----> End of alignment on Subject 
#  - $EVal -----> E value of hit
#  - $BitScore -> Bit Score of hit
    
    my $StartTime = time;
    my @tmpqryid = "";    
    
    # Will try to figure out which one is quicker
    my $TabBlastFile = $_[0];

    # TOO MANY VAR FOR TESTING OUTPUT AND DEBUG
#    my $TooMany = 100;
    my $HitNum = 0;
    my $NumEdges = 0;

    # Concatenated IDs are used to return only a single edge for each
    # Query Subject pair instead of each HSP
    my $CurCatID = ""; # Previous concat ID
    my $PreCatID = ""; # Current concat ID

    # Vars to hold the concatenation of qry and subject


    open (BLASTIN, "<".$TabBlastFile) ||
	die "Can not open BLAST file.$TabBlastFile.\n";
    
    while (<BLASTIN>)
    {
	chomp;                 # Remove newline character
	unless (m/^\#.*/)       # Ignore comment lines, works with -m 9 output
	{
	    $HitNum++;
#	    if ($HitNum > $TooMany) 
#	    {
#		my $EndTime = time;
#		my $TotalTime = $EndTime - $StartTime;
#		print "TOTAL TIME\n\n".$TotalTime."\n";
#		exit;
#	    } 
	    
	    my ($QryId, $SubId, $PID, $Len, 
		$MisMatch, $GapOpen, 
		$QStart,$QEnd, $SStart, $SEnd,
		$EVal, $BitScore) = split(/\t/);
	    # Trim leading white space from Bit score
	    $BitScore =~ s/^\s*(.*?)\s*$/$1/;


# For some strange reason PERL interpets this to be a split
# to an implicit @_ so this does not work
#	    my @tmpqryid = split(/\|/, $QryId ) ||
#		die "Could not split QryID";
#	    $XCrd = $tmpqryid[0];
#	    my @tmphitid = split(/\|/, $SubId ) ||
#		die "Could not split SbjID";
#	    $YCrd = $tmphitid[0];


	    # The following is a temp fix, this will probably
	    # only work as expected for two vars split by a pipe character
	    my ($XCrd, $SomeCrap) = split(/\|/, $QryId);
	    $XCrd = int($XCrd);
	    
	    my ($YCrd, $MoreCrap) = split(/\|/, $SubId);
	    $YCrd = int($YCrd);

            # THE FOLLOWING ONLY PRINTED FOR DEBUG PURPOSES
	    #print "QryID was:".$QryId."\n";
	    #print "HitID was:".$SubId."\n";
	    #print "X:".$XCrd."\n";
	    #print "Y:".$YCrd."\n";
	    #print "$XCrd\t$YCrd\n";

	    #-----------------------------+
	    # WRITE EDGE INFORMATION TO   |
	    # SIF AND *.EA FILES          |
	    #-----------------------------+
	    # The following code is a bit of a kludge that depends
	    # on cytoscape's intepretation of sif files.
	    # It would be better to write this as a subfunction
	    # that would allow the same code to be used for
	    # multiple local alignment programs.
	    # This will be kept here for now


	    # Added 05/16/2007 
	    $CurCatID = $XCrd.":".$YCrd;
	    # Draw a new edge unless it has already been drawn, this assumes
	    # the the order HSPs is ordered such that all HSPs for a subject
	    # sequence will occur in the same order. This only returns the
	    # information for the first hit encountered in the tab blast
	    # output. -- 05/17/2007
	    unless ( $CurCatID =~ $PreCatID)
	    {
		# Increment the number of edges then decide
		# how to draw the graph
		$NumEdges++;         
		
		if ($GraphDir == '0')
		{
		    #-----------------------------+
		    # UNDIRECTED x < y            |
		    #-----------------------------+
		    if ($XCrd < $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		    
		} elsif ($GraphDir == '1'){
		    
		    #-----------------------------+
		    # UNDIRECTED x != y           |
		    #-----------------------------+ 
		    # Using (bl) for both directed makes digraph into unigraph
		    # This is sloppy and depends on cytoscape to interpret
		    if ($XCrd != $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '2'){
		    #-----------------------------+
		    # UNDIRECTED ALL x and y      |
		    #-----------------------------+
		    # All x and y
		    print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
		    print BITOUT $XCrd." (bl) ".$YCrd." = ".$BitScore."\n";
		    print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
		    print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    
		} elsif ($GraphDir == '3') {
		    #-----------------------------+
		    # DIRECTED x < y              |
		    #-----------------------------+
		    # Same output as undirected graph
		    # Cytoscape interpets
		    if ($XCrd < $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '4') {
		    #-----------------------------+
		    # DIRECTED x != y             |
		    #-----------------------------+
		    # The use of bot ensures that the x < y is treated as 
		    # different information then x > y.
		    if ($XCrd < $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		    if ($XCrd > $YCrd)
		    {
			print SIFOUT $XCrd."\tbot\t".$YCrd."\n";
			print BITOUT $XCrd." (bot) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bot) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bot) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '5') {
		    #-----------------------------+
		    # DIRECTED ALL x and y        |
		    #-----------------------------+
		    # All x and y
		    # The use of bot ensures that the x < y is treated as 
		    # different information then x > y.
		    if ($XCrd <= $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		    if ($XCrd > $YCrd)
		    {
			print SIFOUT $XCrd."\tbot\t".$YCrd."\n";
			print BITOUT $XCrd." (bot) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bot) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bot) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} # End of if statements for GraphDir

	    } # End of unless CurCatID =~ PrevCatID    
	    $PreCatID = $CurCatID;

	} # End of unless commented out string
    } # End of while statement for working throught BLAST file

    # Print the number of edges found
    print "EDGES:".$NumEdges."\n";
    
} # End of ParseTabBlast subfunction



sub ParseTabBLAST2Graph
{

# Currently this will only work for undirected graph


# For an undirected graph we want the connected_componenest


# For a directed graph we can choose the
# strongly_connected_components   OR 
# weakly_connected_components


#-----------------------------+
# LOAD A TAB DELIMITED ALL BY |
# ALL BLAST REPORT            |
#-----------------------------+
# This is the tab delimited text file produced from 
# the m8 or m9 option in BLAST, use of the m9 option
# will required that lines starting with # must
# be ignored, I therfore have the regexp (m/^\#.*/).
# The available information:
#  - $QryId ----> Id of the query sequence 
#  - $SubId ----> Id of the subect seqence (from the database
#  - $PID ------> Percent ID
#  - $Len ------> Alignment length
#  - $MisMatch -> Number of mismatches
#  - $GapOpen --> Number of gap openings
#  - $QStart ---> Start of alignment on query sequence
#  - $QEnd -----> End of alignment on query sequence
#  - $SStart ---> Start of aligment on Subject
#  - $SEnd -----> End of alignment on Subject 
#  - $EVal -----> E value of hit
#  - $BitScore -> Bit Score of hit
    
    my $StartTime = time;
    my @tmpqryid = "";    

    # Make dir to hold the fasta files
    my $FastDir = $NetDir."fasta/";
    mkdir $FastDir, 0777 unless (-e $FastDir);
    
    # OPEN FILE TO WRITE NODE ATTRIBUTE FILE
    # FOR THE CLUSTER CLASSIFICATIO
    # This will allow for a check of the connected cluster
    # compared to what I can visualize in  Cytoscape
    my $NA_ConClust = $NetDir."ConnectClust.NA";
    open (CONCLUST, ">$NA_ConClust") ||
	die "Can not open $NA_RepName\n";
    print CONCLUST "ConnectedCluster\n";
    
    # Create and undirected graph object, the default is to
    # create a directed graph
    # Set the scope of the graph object to be local
    my $g;
    
    # Undirected graph
    # The graph must be undirected to return the clusters
    # otherwise we must deal with transitive closure approaches
    $g = Graph->new (directed => 0);

    # Will try to figure out which one is quicker
    my $TabBlastFile = $_[0];

    # TOO MANY VAR FOR TESTING OUTPUT AND DEBUG
#    my $TooMany = 100;
    my $HitNum = 0;  
    my $NumEdges = 0; # The number of edges
    my $NumSing = 0;  # The number of singletons
    my $NumMult = 0;  # The number of clustes with multiple seqs
    #my $Num = 0; 
	
    # Concatenated IDs are used to return only a single edge for each
    # Query Subject pair instead of each HSP
    my $CurCatID = ""; # Previous concat ID
    my $PreCatID = ""; # Current concat ID

    # Vars to hold the concatenation of qry and subject

    #-----------------------------+
    # LOAD ALL POSSIBLE NODES     |
    #-----------------------------+
    open (BLASTIN, "<".$TabBlastFile) ||
	die "Can not open BLAST file.$TabBlastFile.\n";

    print "Adding nodes to Graph object\n";
    $HitNum = 0;
    my $PrevX = "NULL";
    my $NumNodes = 0;
    while (<BLASTIN>)
    {
	chomp;                 # Remove newline character
	unless (m/^\#.*/)       # Ignore comment lines, works with -m 9 output
	{
	    $HitNum++;
	    
	    my ($QryId, $SubId, $PID, $Len, 
		$MisMatch, $GapOpen, 
		$QStart,$QEnd, $SStart, $SEnd,
		$EVal, $BitScore) = split(/\t/);
	    # Trim leading white space from Bit score
	    $BitScore =~ s/^\s*(.*?)\s*$/$1/;

	    # The following is a temp fix, this will probably
	    # only work as expected for two vars split by a pipe character
	    my ($X, $SomeCrap) = split(/\|/, $QryId);
	    $X = int($X);

	    # Only load the records that have not already occured
	    unless ($X =~ $PrevX)
	    {
		print "Adding Node: $X\n";
		$NumNodes++;
		$g->add_vertex( $X );	
	    }
	    $PrevX = $X;
	    
	} # End of unless commented out string

    } # End of while statement for working throught BLAST file


    # Close the BLAST file and then reopen to get the
    # edges.
    close BLASTIN;

    print "$NumNodes Nodes\n";

    #-----------------------------+
    # ADD EDGES                   |
    #-----------------------------+
    # Edge weight can be set by
    # $g->add_weighted_edge ($u, $v, $weight);

    open (BLASTIN, "<".$TabBlastFile) ||
	die "Can not open BLAST file.$TabBlastFile.\n";

    while (<BLASTIN>)
    {
	chomp;                 # Remove newline character
	unless (m/^\#.*/)       # Ignore comment lines, works with -m 9 output
	{
	    $HitNum++;
	    
	    my ($QryId, $SubId, $PID, $Len, 
		$MisMatch, $GapOpen, 
		$QStart,$QEnd, $SStart, $SEnd,
		$EVal, $BitScore) = split(/\t/);
	    # Trim leading white space from Bit score
	    $BitScore =~ s/^\s*(.*?)\s*$/$1/;

	    # The following is a temp fix, this will probably
	    # only work as expected for two vars split by a pipe character
	    my ($XCrd, $SomeCrap) = split(/\|/, $QryId);
	    $XCrd = int($XCrd);
	    
	    my ($YCrd, $MoreCrap) = split(/\|/, $SubId);
	    $YCrd = int($YCrd);

	    # Added 05/16/2007 
	    $CurCatID = $XCrd.":".$YCrd;
	    # Draw a new edge unless it has already been drawn, this assumes
	    # the the order HSPs is ordered such that all HSPs for a subject
	    # sequence will occur in the same order. This only returns the
	    # information for the first hit encountered in the tab blast
	    # output. -- 05/17/2007
	    unless ( $CurCatID =~ $PreCatID)
	    {
		# Increment the number of edges then decide
		# how to draw the graph
		$NumEdges++;         
		
		if ($GraphDir == '0')
		{
		    #-----------------------------+
		    # UNDIRECTED x < y            |
		    #-----------------------------+
		    if ($XCrd < $YCrd)
		    {
			print "Adding Edge $XCrd-->$YCrd\n";
			$g->add_edge($XCrd, $YCrd);
			
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '1'){
		    
		    #-----------------------------+
		    # UNDIRECTED x != y           |
		    #-----------------------------+ 
		    # Using (bl) for both directed makes digraph into unigraph
		    # This is sloppy and depends on cytoscape to interpret
		    if ($XCrd != $YCrd)
		    {
			print "Adding Edge $XCrd-->$YCrd\n";
			$g->add_edge($XCrd, $YCrd);
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '2'){
		    #-----------------------------+
		    # UNDIRECTED ALL x and y      |
		    #-----------------------------+
		    # All x and y
		    print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
		    print BITOUT $XCrd." (bl) ".$YCrd." = ".$BitScore."\n";
		    print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
		    print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    
		} elsif ($GraphDir == '3') {
		    #-----------------------------+
		    # DIRECTED x < y              |
		    #-----------------------------+
		    # Same output as undirected graph
		    # Cytoscape interpets
		    if ($XCrd < $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '4') {
		    #-----------------------------+
		    # DIRECTED x != y             |
		    #-----------------------------+
		    # The use of bot ensures that the x < y is treated as 
		    # different information then x > y.
		    if ($XCrd < $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		    if ($XCrd > $YCrd)
		    {
			print SIFOUT $XCrd."\tbot\t".$YCrd."\n";
			print BITOUT $XCrd." (bot) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bot) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bot) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} elsif ($GraphDir == '5') {
		    #-----------------------------+
		    # DIRECTED ALL x and y        |
		    #-----------------------------+
		    # All x and y
		    # The use of bot ensures that the x < y is treated as 
		    # different information then x > y.
		    if ($XCrd <= $YCrd)
		    {
			print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
			print BITOUT $XCrd." (bl) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bl) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bl) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		    if ($XCrd > $YCrd)
		    {
			print SIFOUT $XCrd."\tbot\t".$YCrd."\n";
			print BITOUT $XCrd." (bot) ".$YCrd." = ".
			    $BitScore."\n";
			print SIGOUT $XCrd." (bot) ".$YCrd." = ".$EVal."\n";
			print PIDOUT $XCrd." (bot) ".$YCrd." = ".$PID."\n";
		    } # End of if XCrd > YCrd
		    
		} # End of if statements for GraphDir

	    } # End of unless CurCatID =~ PrevCatID    
	    $PreCatID = $CurCatID;

	} # End of unless commented out string
    } # End of while statement for working throught BLAST file


    #-----------------------------+
    # PRINT STRINGIFIED GRAPH
    #-----------------------------+
    # This actually prints all nodes connected by edges
    # as separated by -, singletons are printed at the end
    #print "The graph is $g\n";
    
    #-----------------------------+
    # COMPUTE TRANSITIVE CLOSURE  |
    #-----------------------------+
    # Transitive Closure only works for directed graphs
#    print "Computing the TransitiveClosure\n";
#    my $tcm = Graph::TransitiveClosure::Matrix->new($g);
#    
#    # The following returns a concatenated list of vertices
#    #print $tcm->vertices;
#
#    my $TestAns = "IDontKnow";
#    if ($tcm->is_transitive('799','800')){$TestAns = "TRUE";}
#    print "IT:\t799-800\t$TestAns\n";
#    $TestAns = "FALSE";
#    if ($tcm->is_transitive('799','100')){$TestAns = "TRUE";}
#    print "IT:\t799-100\t$TestAns\n"; 
    
    #-----------------------------+
    # PRINT GRAPH SUMMARY INFO    |
    #-----------------------------+
    print "Counting nodes and edges\n";
    my $NumObjNodes = $g->vertices;
    print SUMOUT "Nodes:\t$NumObjNodes\n";
    my $NumObjEdges = $g->edges;
    print SUMOUT "Edges:\t$NumObjEdges\n";
    

    #-----------------------------+
    # ADDTIONAL GRAPH INFO        |
    #-----------------------------+
    #print "Determing average path length";
    #my $apl = $g->average_path_length;
    #print SUMOUT "APL:\t$apl\n";

    #-----------------------------+
    # GRAPH DIAMETER              |
    #-----------------------------+
    # Slow
    # Graph diameter should be determined for
    # each of the clusters
#    print "Determing graph diameter\n";
#    my $gd = $g->diameter;
#    print "GraphDiam:\t$gd\n";
#    print SUMOUT "GraphDiam:\t$gd\n";

    
    
#
#    # The vertices and edges can also be loaded into an array
#    my @V = $g->vertices;
#    @V = sort{$a <=> $b} (@V);
#    foreach my $IndV (@V)
#    {
#	print "$IndV\n";
#    }

    # Average path length is the average shortest path
    
    #-----------------------------+
    # AVERAGE PATH LENGTH         |
    #-----------------------------+
#    print "Determing average path lengths\n";
#    my @V = $g->vertices;
#    @V = sort{$a <=> $b} (@V);
#    foreach my $IndV (@V)
#    {
#	print "Determing length for $IndV\n";
#	my $apl = $g->average_path_length($IndV) || "NULL"; 
#	print "$IndV\t$apl\n";
#	print SUMOUT "$IndV\t$apl\n";
#    }
    

    #//////////////////////////////////////////////////////////////
    # NOTE: Can fetch the connect_components_by_vertix
    #       to create subgraphs that can be analyzed separately
    #$i = $g->connected_component_by_vertex($v)
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    #-----------------------------+
    # WORK WITH THE CONNECTED     |
    # COMPONENETS OF THE GRAPH    |
    #-----------------------------+
    # For a unidrected graph we can do the connected components

    print "Determining the connected components\n";
    my $NumFam = 0;
    my $NumClust = 0;
    my @cc = $g->connected_components();
    
    foreach my $ComList (@cc)
    {
	$NumClust++;
	print "CLUST_$NumClust:\n";
	my $NumComp = 0;
	
			
	# TotComp is the total number of components
	my $TotComp = @$ComList;
	print "\tComp:\t$TotComp\n";

	
	if ($TotComp > 1)
	{
	    
	    $NumMult++;

	    my $ClustFastOut = $FastDir.
		"Seqs_CLUST$NumClust.fasta";
	    
	    open (CLUSTFASTA, ">".$ClustFastOut) ||
		die "Could not open file for writing:\n$ClustFastOut\n";

	    # $IndComp in the individual component of the graph
	    for my $IndComp (@$ComList)
	    {
		$NumComp++;
		
		#-----------------------------+
		# DETERMINE AVERAGE PATH      |
		# LENGTH FOR THE NODE         |
		#-----------------------------+
		# This is very slow
		#print "\ttDeterming average path length\n";
		#my $apl = $g->average_path_length($IndComp) || "NULL"; 
		#print "\t$IndComp\t$apl\n";
		#print SUMOUT "\t$IndComp\t$apl\n";
		
		#-----------------------------+
		# PRINT NODE ATTRIBUTE OUT    | 
		#-----------------------------+
		#print REPNAME $OutName."=".$Name."\n";
		print CONCLUST $IndComp."=CLUST_".$NumClust."\n";

		#-----------------------------+
		# FETCH SEQ ID DATA           |
		#-----------------------------+
		my $FetchIDQry = "SELECT id FROM tblSeqData".
		    " WHERE rownum=$IndComp";
		my $IDsth = $dbh->prepare($FetchIDQry);
		$IDsth->execute() || 
		    print "\nERROR: Can not fetch ID for $IndComp\n";
		my @ID_row_ref = $IDsth->fetchrow_array;
		my $ID = $ID_row_ref[0] || "0";

		#-----------------------------+
		# FETCH SEQUENCE DATA         |
		#-----------------------------+
		my $FetchSeqQry = "SELECT seq FROM tblSeqData".
		    " WHERE rownum=$IndComp";
		my $Seqsth = $dbh->prepare($FetchSeqQry);
		$Seqsth->execute() || 
		    print "\nERROR: Can not fetch seq for $IndComp\n";
		my @Seq_row_ref = $Seqsth->fetchrow_array;
		my $SeqString = $Seq_row_ref[0] || "0";
		
		#-----------------------------+
		# WRITE TO FASTA FILE         |
		#-----------------------------+
		#print CLUSTFASTA ">$ID\n";
		# Changed 06/27/2007 to include the rown
		print CLUSTFASTA ">$IndComp|$ID\n";
		print CLUSTFASTA "$SeqString\n";
		
	    } # End of for each individual component
	    
	    # Close the FASTA output file
	    close CLUSTFASTA;
	    
	} else { 
			
	    $NumSing++;

	    for my $IndComp (@$ComList)
	    {
		#-----------------------------+
		# OPEN SINGLETON FASTA FILE   |
		#-----------------------------+
		my $SingFastOut = $FastDir.
		    "Singletons.fasta";
		
		# This will append data to an already existing file
		open (SINGFASTA, ">>".$SingFastOut) ||
		    die "Could not open file for writing:\n$SingFastOut\n";

		#-----------------------------+
		# FETCH SEQ ID DATA           |
		#-----------------------------+
		my $FetchIDQry = "SELECT id FROM tblSeqData".
		" WHERE rownum=$IndComp";
		my $IDsth = $dbh->prepare($FetchIDQry);
		$IDsth->execute() || 
		    print "\nERROR: Can not fetch ID for $IndComp\n";
		my @ID_row_ref = $IDsth->fetchrow_array;
		my $ID = $ID_row_ref[0] || "0";
		
		#-----------------------------+
		# FETCH SEQUENCE DATA         |
		#-----------------------------+
		my $FetchSeqQry = "SELECT seq FROM tblSeqData".
		    " WHERE rownum=$IndComp";
		my $Seqsth = $dbh->prepare($FetchSeqQry);
		$Seqsth->execute() || 
		    print "\nERROR: Can not fetch seq for $IndComp\n";
		my @Seq_row_ref = $Seqsth->fetchrow_array;
		my $SeqString = $Seq_row_ref[0] || "0";
		
		#-----------------------------+
		# WRITE TO FASTA FILE         |
		#-----------------------------+
		#print SINGFASTA ">CLUST_$NumClust|$ID\n";
		# Changed 06/27/2007
		print SINGFASTA ">CLUST_$NumClust|$IndComp|$ID\n";
		print SINGFASTA "$SeqString\n";
		
	    } # END of for each indcomp
	    
	} # End of TotComp > 1
	
	# print "\n";
	#print "\tComp:\t$NumComp\n";
	
	# Print the number of components in each cluster to the 
	# summary output file
	print SUMOUT "CLUST_$NumClust:\t$NumComp\n";
	
    } # End of for IndComp

    print "SING:\t$NumSing\n";
    print "GRPS:\t$NumMult\n";
    print SUMOUT "SING:\t$NumSing\n";
    print SUMOUT "GRPS:\t$NumMult\n";
    
} # End of ParseTabBlast subfunction


sub LoadAllByAll
{
#-----------------------------+
# LOAD RESULT OF ALL BY ALL   |
# BLAST TO THE DATABASE       |
# AND TEMPORARILY TO THE BLAST|
# MATRIX 2D ARRAY             |
#-----------------------------+
#
# THIS IS THE CORE OF THE GRAPH BASED APPROACH
# TO EXPLORING THE RESULTS OF AN ALL BY ALL BLAST
#
    my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
					  '-file'   => $AllByAll, 
					  '-signif' => $A_MaxE, 
					  '-min_query_len' => $A_MinQryLen,
					  '-score' => $A_MinScore ) 
	||
        die "Could not open BLAST input file:\n$AllByAll.\n";
    
    while ($BlastResult = $BlastReport->next_result())
    {
	my ($LoadRec); # Set scope for this subfunction

	my $NumHits = $BlastResult->num_hits;
	my $QryName = $BlastResult->query_name;

	print SUMOUT $QryName."\t".$NumHits."\n"; 

	my @tmpqry = split(/\|/, $QryName );
	$XCrd = $tmpqry[0];                  # First part of name is number
                                             # representing X position
	my $DbQryID = $tmpqry[1];            # Second part of name is unique
                                             # id to allow lookup in the DB

	#-----------------------------+
	# SELECT NODE ATTRIBUTES FROM |
	# THE DATABASE AND PRINT THE  |
	# RESULTS TO THE NODE         |
	# ATTRIBUTE FILES             |
	#-----------------------------+
	#my $QryBACID = &FetchSeqInfo("tblSeqData", "bac", $DbQryID);
	#print "\tBAC: $QryBACID\n";
	#print BACOUT $XCrd."=".$QryBACID."\n";

	my $QrySeqData =&FetchSeqInfo("tblSeqData", "seq", $DbQryID);
	print SEQOUT $XCrd."=".$QrySeqData."\n";


	if ($XCrd > $MaxVal) {$MaxVal = $XCrd;}
	print $XCrd."\n";
	while ( $BlastHit = $BlastResult->next_hit())
	{
	    my $HitName = $BlastHit->name();
	    my $BitScore = $BlastHit->bits();
	    my $HitSig = $BlastHit->significance();

	    my @tmphit = split(/\|/, $HitName );
	    $YCrd = $tmphit[0];


	    #-----------------------------------------------------------+
	    #                                                           |
	    # DEFINE THE GRAPH                                          |
	    #                                                           |
	    # ALGORITHM WORK NEEDED HERE                                |
	    # 01/08/2007                                                |
	    # This graph defintion creates the nodes as well as the 
	    # edges by defining a list of directed edges connecting 
	    # nodes.
	    #                                                           |
	    #-----------------------------------------------------------+
	    # Print output for each all by all hit to the 
	    # sif output file. Only do this for x > Y
	    # to prevent redundancy in the netwrok file.
	    # This should take the graph building variable as an
	    # argument.
	    
	    # I am using XCrd and YCrd instead of i and j in the matrix
	    # terminology.

	    # The following only does the top triangle of the
	    # adjacency matrix.

	    # The following only does the top part of the homology
	    # matrix
	    if ($XCrd > $YCrd)
	    {
		print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
		print BITOUT $XCrd." (bl) ".$YCrd." = ".$BitScore."\n";
		print SIGOUT $XCrd." (bl) ".$YCrd." = ".$HitSig."\n";

	    }

	    #-----------------------------+
	    # The following drew both lines
	    # but it drew them in the same
	    # direction ... I am possible
	    # messed up with the setting
	    # of the X and Y Crds above
	    # The following does the bottom part of the homology matrix
	    # Using this with the above will create a dataset that 
	    # Cytoscape will interpret as a undirected graph
	    if ($XCrd < $YCrd)
	    {
		# I thought that this should do YCrd Then XCrd 
                # to show path in other direction but this is not printing
		# like I thought it would.
		# Changed bl to bot to show this is the bottom part of
		# the matrix.
		print SIFOUT $XCrd."\tbl\t".$YCrd."\n";
		print BITOUT $XCrd." (bl) ".$YCrd." = ".$BitScore."\n";
		print SIGOUT $XCrd." (bl) ".$YCrd." = ".$HitSig."\n";

	    }

	    $ResultCount++;

	    if ($YCrd == $XCrd)
	    {
		# IDENTIFY SELF HITS IN hit_cat FIELD
		$LoadRec = "INSERT INTO ".$tblAllByAll.
		    " (qry_num, hit_num, qry_id, num_hits, hit_cat) ".
		    " VALUES (".
		    " '".$XCrd."',".
		    " '".$YCrd."',".
		    " '".$QryName."',".
		    " '".$NumHits."',".
		    " 'self'".
		    ")";
	    }else{
		# DO NOT IDENTIFY hit_cat
		$LoadRec = "INSERT INTO ".$tblAllByAll.
		    " (qry_num, hit_num, qry_id, num_hits) ".
		    " VALUES (".
		    " '".$XCrd."',".
		    " '".$YCrd."',".
		    " '".$QryName."',".
		    " '".$NumHits."'".
		    ")";
	    }

	    $dbh->do($LoadRec);

	    # TEMP EXIT FOR DEBUG
	    # if ($ResultCount == 100){exit;}

# 01/10/2007
# COMMENTED OUT THE FOLLOWING HSP WORK TO TRY TO
# FIX DIGRAPH OUTPUT
#	    my $HspNum = "0"; # Initialize HSP Num
#	    while ( $HSP = $BlastHit->next_hsp())
#	    {
#
#		$HspNum++;
#		# May need to add hsp# 
#		# where # is 1 to num_hsps for multiple hsps
#		if (! $quiet)
#		{
#		    print "\t".$HSP->frac_identical."\n";
#		    print "\t".$HSP->length('total')."\n";
#		}
#
#		if ($XCrd > $YCrd)
#		{
#
#		    print HSPOUT $XCrd."\thsp$HspNum\t".$YCrd."\n";
#		    print HSPFI $XCrd." (hsp$HspNum) ".$YCrd." = ".
#			$HSP->frac_identical."\n";
#
#		}
#
#	    } #End of Next HSP







	} # End of Next BLAST hit

    } # End of BLAST result next resul   
} # End of Load All By All Subfunction


sub does_table_exist
{
#-----------------------------+
# CHECK IF THE MYSQL TABLE    |
# ALREADY EXISTS              |
#-----------------------------+
# CODE FROM
# http://lena.franken.de/perl_hier/databases.html
# Makes use of global database handle dbh


    #my ($dbh, $whichtable) = @_;
    my ($whichtable) = @_;
    my ($table,@alltables,$found);
    @alltables = $dbh->tables();
    $found = 0;
    foreach $table (@alltables) {
	$found=1 if ($table eq "`".$whichtable."`");
    }
    # return true if table was found, false if table was not found
    return $found;
}


sub how_many_records
{
#-----------------------------+
# COUNT HOW MANY RECORDS      |
# EXIST IN THE MYSQL TABLE    |
#-----------------------------+
# CODE FROM
# http://lena.franken.de/perl_hier/databases.html
# Makes use of global database handle dbh


    #my ($dbh, $whichtable) = @_;
    my ($whichtable) = @_;
    my ($result,$cur,$sql,@row);

    $sql = "select count(*) from $whichtable";
    $cur = $dbh->prepare($sql);
    $cur->execute();
    @row=$cur->fetchrow;
    $result=$row[0];
    $cur->finish();
    return $result;

}


sub FetchSeqInfo
{
#-----------------------------+
# FETCH THE BAC ID GIVEN THE  |
# TABLE AND ID TO FETCH       |
# 06/21/2006                  |
#-----------------------------+
# Returns UNK for search queries that 
# could not be found in the database
    my $result;
    my ($whichtable, $searchField, $searchID ) = @_;
    #select * from tblSeqData where rownum = '1'

    #my $GetBACSQL = "SELECT * FROM ".$whichtable.
	#" WHERE id='".$searchID."'";

    my $GetBACSQL = "SELECT ".$searchField." FROM ".$whichtable.
	" WHERE id='".$searchID."'";

    my $cur = $dbh->prepare($GetBACSQL);
    $cur->execute();
    my @row=$cur->fetchrow;
    $result=$row[0] || "UNK";
    $cur->finish();
    return $result;
}


sub UserFeedback
{
#-----------------------------+
# USER FEEDBACK SUBFUNCTION   |
#-----------------------------+

  my $Question = $_[0];
  my $Answer;

  print "\n$Question \n";

  while (<>)
    {
      chop;
      if ( ($_ eq 'y') || ($_ eq 'Y') || ($_ eq 'yes') || ($_ eq 'YES') )
	{
	  $Answer = "Y";
	  return $Answer;
	}
      elsif ( ($_ eq 'n') || ($_ eq 'N') || ($_ eq 'NO') || ($_ eq 'no') )
	{
	  $Answer = "N";
	  return $Answer;
	}
      else
	{
	  print "\n$Question \n";
	}
    }

}

sub LaunchCytoscapeOld
{
#-----------------------------+
# LAUNCH CYTOSCAPE WITH THE   |
# EDGE ATTRIBUTES AND NODE    |
# ATTRIBUTES THAT HAVE BEEN   |
# CREATED                     |
# Cytoscape V 2.2             |
#-----------------------------+

    # Edge attributes are -j
    # Node attributes are -n

    print "Launching Cytoscape\n";
    my $CyPath = "/home/jestill/Apps/Cytoscape-v2.2/cytoscape.jar";
    my $VpPath = "/home/jestill/Apps/Cytoscape-v2.2/vizmap.props";
    my $PlugPath = "/home/jestill/Apps/Cytoscape-v2.2/plugins";

    # SIF NETWORK FILE
    my $SifFile = $_[0]; # The -i variable
    # NODE ATTRIBUTE FILES
    my $NAs = $NA_BACOut." ".$NA_RepClass." ".$NA_RepName." ".$NA_SeqData;
    # EDGE ATTRIBUTE FILES
    my $EAs = $EA_BitScore." ".$EA_Sig;
#    my $SysCmd = 'java -Xmx512M -jar '.$CyPath.' -i '.$SifFile.
    my $SysCmd = 'java -Xmx2048M -jar '.$CyPath.' -i '.$SifFile.
	    ' -vp '.$VpPath.' cytoscape.CyMain'.
	    ' --JLD '.$PlugPath.
	    #' -n '.$NA_BACOut.
	    ' -n '.$NAs.
	    ' -j '.$EAs.
	    ' $*';

    print "\n\nCMD IS:\n$SysCmd\n\n";
    system ( $SysCmd );

}


sub LaunchCytoscapeNew
{
#-----------------------------+
# NEW VERSION OF THE LAUNCH   |
# CYTOSCAPE FOR VERSION 2.3   |
#-----------------------------+

    # Relevant Cytoscape command line arguments
    # Edge attributes are    -e
    # Node attributes are    -n
    # Network file           -N
    # Cytoscape plugin path  -p
    # Vizmap properties file -V
    # Cytoscape Properties   -P
    # 

    print "Launching Cytoscape\n";
    my $CyPath = "/home/jestill/Apps/Cytoscape-v2.3/cytoscape.jar";
    #my $VpPath = "/home/jestill/Apps/Cytoscape-v2.2/vizmap.props";
    my $PlugPath = "/home/jestill/Apps/Cytoscape-v2.3/plugins";

    # SIF NETWORK FILE
    my $SifFile = $_[0]; # The -N variable
    # NODE ATTRIBUTE FILES
    my $NAs = $NA_BACOut." ".$NA_RepClass." ".$NA_RepName." ".$NA_SeqData;
    # EDGE ATTRIBUTE FILES
    my $EAs = $EA_BitScore." ".$EA_Sig;
#    my $SysCmd = 'java -Xmx512M -jar '.$CyPath.
    my $SysCmd = 'java -Xmx2048M -jar '.$CyPath.
	' -N '.$SifFile.
	#' -V '.$VpPath.
	' cytoscape.CyMain'.
	' -p '.$PlugPath.
	' -n '.$NAs.
	' -e '.$EAs.
	' $*';

    print "\n\nCMD IS:\n$SysCmd\n\n";
    system ( $SysCmd );

}

sub LaunchCytoscape_2_4
{
# Added
# 04/04/2007
#-----------------------------+
# NEW VERSION OF THE LAUNCH   |
# CYTOSCAPE FOR VERSION 2.4   |
#-----------------------------+
# The following is still not working

    # Relevant Cytoscape command line arguments
    # Edge attributes are    -e
    # Node attributes are    -n
    # Network file           -N
    # Cytoscape plugin path  -p
    # Vizmap properties file -V
    # Cytoscape Properties   -P
    # 

    print "Launching Cytoscape\n";

    #my $VpPath = "/home/jestill/Apps/Cytoscape-v2.2/vizmap.props";
#    my $CyPath = "/home/jestill/Apps/Cytoscape-v2.3/cytoscape.jar";
#    my $PlugPath = "/home/jestill/Apps/Cytoscape-v2.3/plugins";

    my $CyPath = "/home/jestill/Apps/Cytoscape_v2.4.0/cytoscape.jar";
    my $PlugPath = "/home/jestill/Apps/Cytoscape_v2.4.0/plugins";

    # SIF NETWORK FILE
    my $SifFile = $_[0]; # The -N variable
    # NODE ATTRIBUTE FILES
#    my $NAs = $NA_BACOut." ".$NA_RepClass." ".$NA_RepName." ".$NA_SeqData;

    my $NAs = $NA_RepClass." ".$NA_RepName." ".$NA_SeqData;
    # EDGE ATTRIBUTE FILES
#    my $EAs = $EA_BitScore." ".$EA_Sig;
# Adding edge flag here
    my $EAs = $EA_BitScore." -e ".$EA_Sig;
    my $javapath = "/usr/java/jre1.5.0_06/bin/java";

    my $SysCmd = $javapath.
	' -Dswing.aatext=true'.
	' -Xmx512M'.
	' -jar '.$CyPath.
	' cytoscape.CyMain'.
	' -N '.$SifFile.
	#' -V '.$VpPath.
	' -n '.$NAs.
	' -e '.$EAs.
	' -p '.$PlugPath.
	' $*';

    print "\n\nCMD IS:\n$SysCmd\n\n";
    system ( $SysCmd );

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

    my $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    my @row=$cur->fetchrow;
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

    my $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    my @row=$cur->fetchrow;
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
    my $code = substr ($qry,4,7) || 
	"BAD FORMAT"; # If the naming is not as I expect
                           # use the Unexpected Code. This will happen
                           # when the substring is outof bounds.
    

    # IF THE TIGR FORMAT IS MALFORMED THEN SHOW
    # THE QUERY STRING THAT WAS SEND FOR LOOKUP
    if ($code =~ "BAD FORMAT")
    {print "BAD TIGR FORMAT\t".$qry."\n";}

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

sub DrawXYPlot 
{
#-----------------------------+
# DRAW COLOR VALUED X-Y PLOT  |
# OF THE SPARSE MATRIX        |
# DESCRIBING THE ALL BY ALL   |
# BLAST RESULTS               |
#-----------------------------+

    #-----------------------------------------------------------+
    # CREATE GRAPHIC OF THE PARSED AND CATEGORIZED              |
    # ALL BY ALL BLAST                                          |
    #-----------------------------------------------------------+

    #-----------------------------+
    # CREATE GD IMAGE OBJECT      |
    # TO HOLD THE GRAPH           |
    #-----------------------------+
    my $MaxX = &GetMaxVal("qry_num", $tblAllByAll);
    my $MaxY = &GetMaxVal("hit_num", $tblAllByAll);
    $MaxX = $MaxX*$xsc;
    $MaxY = $MaxY*$ysc;

    my $img = new GD::Image($MaxX, $MaxY);

    #-----------------------------+
    # ALLOCATE COLORS FOR THE     |
    # REPEAT CATEGORIES           |
    #-----------------------------+
    my $white = $img->colorAllocate      ( 255,  255,255);
    my $colSep = $img->colorAllocate     (   0,   0,   0);
    my $colUnk = $img->colorAllocate     (   0,   0,   0);
    my $colLTR = $img->colorAllocate     ( 255,   0,   0);
    my $colNonLTR = $img->colorAllocate  (   0,   0, 255);
    my $colMITE = $img->colorAllocate    (   0, 255, 255);
    my $colCACTA = $img->colorAllocate   (   0, 255,   0);
    my $colSelf = $img->colorAllocate    ( 200, 200, 200);  
    my $colOryza = $img->colorAllocate   (   0, 255,   0);
    my $colWess = $img->colorAllocate    ( 255,   0,   0);
    my $colSanMig = $img->colorAllocate  (   0,   0, 255);
    my $colZea = $img->colorAllocate     ( 200, 200,   0);

    #-----------------------------+
    # DRAW BLACK DOTS ON THE      |
    # GRAPH FOR EACH PAIRED HIT   |
    #-----------------------------+
    print "ADDING POINTS TO GRAPHIC OUTPUT\n";

    my $NumRecords = &how_many_records($tblAllByAll);

    #-----------------------------+
    # LOOP THROUGH THE DB OF ALL  |
    # BY ALL HITS AND DRAW HITS   |
    # WITH CLASS INDICATED BY     | 
    # COLOR                       |
    #-----------------------------+
    for (my $i = 1; $i <= $NumRecords; $i++) 
    {
	my $Query = "SELECT * FROM ".$tblAllByAll." WHERE rownum=".$i;
	my $sth = $dbh->prepare($Query);
	$sth->execute();
	while(my @row_ref = $sth->fetchrow_array)
	{
	    $XCrd = $row_ref[1] * $xsc;
	    $YCrd = $row_ref[2] * $ysc;
	    my $QryCat = $row_ref[4] || "unk";
	    my $HitCat = $row_ref[5] || "unk";


	    # May need to add general categorization routine here
	    # to translate from 
	    if ($HitCat =~ "self")
	    {
		$img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colSelf);
	    }else{
		if ($QryCat =~ "non-LTR-retrotransposon")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colNonLTR);
		    $img->fill($XCrd, $YCrd, $colNonLTR);
		}
		elsif($QryCat =~ "LTR-retrotransposon")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colLTR);
		    $img->fill($XCrd, $YCrd, $colLTR);
		}
		elsif ($QryCat =~ "CACTA-transposon")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colCACTA);
		    $img->fill($XCrd, $YCrd, $colCACTA);
		}
		elsif ($QryCat =~ "MITE-foldback element")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colMITE);
		    $img->fill($XCrd, $YCrd, $colMITE);
		}
		elsif ($QryCat =~ "UNK-OryzaRepeat")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colOryza);
		}
		elsif ($QryCat =~ "UNK-Wessler")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colWess);
		}
		elsif ($QryCat =~ "UNK-SanMiguel")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colSanMig);
		}
		elsif ($QryCat =~ "UNK-ZeaRepeat")
		{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colZea);
		}else{
		    $img->arc($XCrd, $YCrd, $pxs, $pxs, 0, 360, $colUnk);
		} # End of for deciding colors for nonSelfHits
	    } # End of deciding between self hits and otherwise
	}
    }

    #-----------------------------+
    # WRITE OUTPUT TO FILE        |
    #-----------------------------+
    open (OUTFILE, ">$GraphOut") ||
	die "Can not open $GraphOut";
    binmode OUTFILE;
    print OUTFILE $img->png;
    close OUTFILE;

    #-----------------------------+
    # WRITE THE OUTPUT TO MONITOR |
    # USING IMAGEMAGICK DISPLAY   |
    #-----------------------------+
    my $png_data = $img->png;
    open (DISPLAY,"| display -") || die;
    binmode DISPLAY;
    print DISPLAY $png_data;
    close DISPLAY;


}

#-----------------------------------------------------------+
# FUNCTIONS                                                 |
#-----------------------------------------------------------+

#-----------------------------+
# FUNCTIONS TO REMOVE         |
# WHITE SPACES FROM A STRING  |
#-----------------------------+

# Removoe trailing and 
# Code from:
# http://www.somacon.com/p114.php

# Perl trim function to remove whitespace from the start and end of the string

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

sub PrintHelp
{
    
    my $FullUsage = "EXAMPLES:\n".
	"jabablast.pl -i AllByAllBlast.blo -r RepElementBlast.blo\n". 
	"-o OutputDir -u UserName -d DatabaseName -Q\n\n".
	"OR \n".
	"jabablast.pl -C ConfigFile.jcfg\n\n".
	"+-----------------------------------------------------------+\n".
	"| REQUIRED ARGUMENTS                                        |\n".
	"+-----------------------------------------------------------+\n".
	" -i The AllByAllBlast File to parse [string]\n".
	" -r The BLAST results against known repeats [string]\n".
	" -o The path of the output directory [string]\n".
	" -u The user name for connection to the database [string]\n".
	" -d The database name to use for the database connection. [string]\n".
	"    This database MUST already exist in you MySQL database.\n".
	"  OR\n".
	" -Z Config file that provides the above information. [string]\n".
	"+-----------------------------------------------------------+\n".
	"| ADDITIONAL ARGUMENTS                                      |\n".
	"+-----------------------------------------------------------+\n".
	" -f The name of the output sif file. [string]\n".
	"    default = Network\n".
	" -G Produce GraphViz output".
	" -Q Run program in quiet mode. [boolean flag]\n".
	"    default = Not quiet\n".
	" -B Use the database for classification information [boolean flag]\n".
	"    Many of the following options require a database\n".
	"    Without a database, only an all by all BLAST can be \n".
	"    visualized without classification into repeat categories\n".
	" -C Open Cytoscape to view the output [boolean flag]\n".
	"    default = Cytoscape not opened.\n".
	"+-----------------------------+\n".
	"| MATRIX OPTIONS              |\n".
	"+-----------------------------+\n".
	" The matrix creation portion of jabablast may be dropped\n".
	" in the future to simplify the command line options.\n".
	" -M Create visualization of all by all BLAST matrix [boolean flag]\n".
	"    default = F\n".
	"    This is useful for visualiztion of small sets of ordered seqs.\n".
	" -x Matrix X axis scaling factor [positive integer]\n".
	"    default = 2\n".
	" -y Matrix Y axis scalign factor [positive integer]\n".
	"    default = 2\n".   
	" -p Matrix Pixel size [positive integer]\n".
	"    default =3\n".
	"+-----------------------------+\n".
	"| DATABASE OPTIONS            |\n".
	"+-----------------------------+\n".
	" -a The table name to use for the all by all blast\n".
	"    in the databse[string]\n".
	"    default = tblAllByAll\n".
	" -c The table name to use for the known repeat classification\n".
	"    table in the database [string]\n".
	"    default = tblRepeatID\n".
	"+-----------------------------+\n".
	"| BLAST OPTIONS               |\n".
	"+-----------------------------+\n".
	" -m blast output aligment view. [Integer 0,8, or 9]\n".
	"    This matches the -m flag from NCBI blastall\n".
	"    -m 0 = pairwise\n".
	"    -m 8 = tabular (default)\n".
	"    -m 9 = tabular with comment lines\n".
	"    default = -m 8\n".  
	" -e Max e value for All by All BLAST results.\n".
	"    default = 1.0e-05\n".
	" -s Minimum bit score for All by All BLAST results.[Integer]\n".
	"    default = 50\n".
	" -l Mimimum length of the query sequence in All by All BLAST\n".
	"    default = 150\n".
	" -r The BLAST results against known repeats [String]\n".
	" -E Max e-value for BLAST against repeat database.\n".
	"    default = 1.0e-03".
	" -S Minimum bit score for BLAST against repeat database.\n".
	"    default = 50\n".
	" -L Minimum length of the query sequence against the repeat".
	" database.\n".
	"    default = 50\n".
	"+-----------------------------+\n".
	"| GRAPH OPTIONS               |\n".
	"+-----------------------------+\n".
	" -d Graph edge direction.\n".
	"    Currently only supported for tab delim blast\n".
	"    0 = undirected (i < j)\n".
	"    1 = undirected (i =! j)\n".
	"    2 = undirected (all i,j)\n".
	"        Reciprocal hits reduced to one edge.\n".
	"    3 = directed (i < j)\n".
	"    4 = directed (i =! j)\n".
	"    5 = directed (all i,j)\n".
	" -n Graph node attributes\n".
	"    The graph node attributes to\n".
	" -g Graph algorithm.\n".
	"    Not currently implemented.\n".
	" -k K best blast hits used\n.".
	"    Not currently implemented.\n".
	" -w Graph edge weight\n".
	"    0 unweighted\n".
	"    Graph edge algorithm\n".             # ie K best blast etc.
	"\n";

print $FullUsage;

}

=head1 NAME

jaba_blast.pl - Parse All by All blast results 

=head1 SYNOPSIS

    jabablast.pl -i AllByAll.blo -r RepBlast.blo -o OutDir -u UserName
                 -b dbName -f Network.sif

=head1 DESCRIPTION

The jabablast (I<jamies all-by-all blast>) program allows for
visualizing all by all blast results for use in the identification and
analysis of the repetitive fraction of genome sequence data.
The program can create a graphical representation of the matrix 
of all-by-all BLAST results as well as text files describing the edges 
and nodes of a graph that can be visualized in the Cytoscape 
(L<http://www.cytoscape.org>) graph visualization program. 
This program is a component of RepMiner.

An important component of jabablast is the ability to parse BLAST restuls 
from various repeat databases into a single classification scheme.
The databases that are used by jabablast are focused on databases
that are relevant to plants. 
The repeat databases that jabablast can use include:

=over 2

=item * TREP (L<http://wheat.pw.usda.gov/ITMI/Repeats/>)

=item * RepBase (L<http://www.girinst.org/repbase/update/index.html>)

=item * SanMiguel (L<http://www.genomics.purdue.edu/~pmiguel/projects/retros/>)

=item * TIGR (L<http://www.tigr.org/tdb/e2k1/plant.repeats/>) 

=back

=head1 ARGUMENTS

Not all of the following arguments are currently implemented. 

=head2 Required Arguments

=over 2

=item -i AllByAllBlast.blo

Path to the AllByAllBlast File to parse.

=item -r RepBlast.blo

The BLAST results against known repeats.

=item -o OutDir

The path of the output directory.

=item -u UserName

The user name for connection to the database.

=item -d dbName

The database name to use for the database connection.
This database MUST already exist in you MySQL database.

=item -f Network.sif

The name of the output network text file.[string]
default = Network.sif

=back

=head2 Boolean Arguments

=over2

=item -Q

Run program in quiet mode. [boolean flag]
default = Not quiet.

=item -B

Use the database for classification information [boolean flag]
Many of the following options require a database
Without a database, only an all by all BLAST can be 
visualized without classification into repeat categories

=item -C

Open Cytoscape to view the output [boolean flag]
default = Cytoscape not opened.

=item -G

Produce graph with graphviz format. This will probably
get moved to a separate program. There are a number 
of variables that can be set with GraphViz that would be
useful to set at the command lined. [Added 05/16/2007]

=item -M 

Create visualization of all by all BLAST matrix.
default = Matrix not created.
This is useful for visualiztion of small sets of ordered seqs.

=back

=head2 Matrix Options

The matrix creation portion of jabablast may be dropped
in the future to simplify the command line options.

=over 2

=item -x 2

Matrix X axis scaling factor [positive integer]
default = 2

=item -y 2

Matrix Y axis scaling factor [positive integer].
default = 2

=item -p 3

Matrix Pixel size [positive integer].
default =3

=back

=head2 Database Options

=over 2

=item -a tblAllByAll

The table name to use for the all by all blast
in the databse[string]
default = tblAllByAll

=item -c tblRepeatID

The table name to use for the known repeat classification
table in the database [string]
default = tblRepeatID

=back

=head2 BLAST Options

=over 2

=item -m 8

blast output aligment view. [Integer 0,8, or 9]
This matches the -m flag from NCBI blastall
default = -m 8

-m 0 = pairwise

-m 8 = tabular (default)

-m 9 = tabular with comment lines

=item -e 1.0e-05

Max e value for All by All BLAST results.
default = 1.0e-05

=item -s 50

Minimum bit score for All by All BLAST results.[Integer]
default = 50

=item -l 150

Mimimum length of the query sequence in All by All BLAST
default = 150

=item -E 1.0e-03 

Max e-value for BLAST against repeat database.
default = 1.0e-03 

=item -S 50

Minimum bit score for BLAST against repeat database.
    default = 50

=item -L 50

Minimum length of the query sequence against the repeat database.
default = 50

=back

=head2 Graph Options

=over

=item -d 

Graph edge direction. This is currently only supported for tab delim blast:

-d 0 = undirected (i < j)

-d 1 = undirected (i != j)

-d 2 = undirected (all i,j)
       Reciprocal hits reduced to a single edge.

-d 3 = directed (i < j)

-d 4 = directed (i != j)

-d 5 = directed (all i,j)
       Reciprocal hits drawn as two edges.

=item -n Graph node attributes

The graph node attributes to use for node classification.
Not currently implemented.

=item -g Graph algorithm

Not currently implemented.

=item -k K best blast hits used

Not currently implemented.

=item -w 0

Graph edge weight option:

-w 0 unweighted 

=back

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 06/14/2006

UPDATED: 06/27/2007

=cut

#-----------------------------------------------------------+
# DETAILED HISTORY                                          |
#-----------------------------------------------------------+
#
# 06/14/2006
# - Program started, wrote overview and beguan testing.
# - Played with GD to get syntaxt correct
# - Working on BLAST parser to load results to an array
# - Took output from basic AbyA blast and visualized as
#   x-y dot graph in PNG format.
# - PNG format output to screen working
# - PNG format output to file working
#
# 06/16/2006
# - DbSetup subfunction
# - Added basic database subfunctions:
#    -how_many_records
#    -does_table_exist
# - Added UserFeedback subfunction
# - Created test BLAST datbase an output file from subset
#   of asgr_c data. 199 query seqs.
#
# 06/19/2006
# - Converted load All by All from main function to 
#   a subfunction
# - Started work on LoadRepClass subfucntion for TREP
#   Database
# - Added string manipulation functions
#    -ltrim
#    -rtrim
#    -trim
# - Added update UpdateRepCat
# - Added GetMaxValue subfunction to determine the maximum
#   x and Y positions for the graph
# - Changed code to pull x,y positions from the database
#   instead of the BlastMatrix array
# 
# 06/20/2006
# - Added lines to delineate BACs for the ASGR data
# - Added index to AllByAll and RepClass tables to significantly
#   speed up the query across tables
# - Added ability to parse multiple repeat databases to add
#   as categories. Currently will just parse to UNK-RepDB Name.
# - Testing the ability to parse blast output files of the
#   concatenated output of several repeat databases
# - Added the number of hits to the information that is parsed
#   from the all by all blast.
# - Output number of hits to summary output file
# - Output number of hits to database
# - Added sif output file for viewing of the all by all
#   BLAST results in a network viewing program
#
# 06/21/2006
# - Fixed the export of SIF format network output
# - Added FetchSeqInfo Subfunction, this is use to fetch
#   the BAC from the tblSeqInfo table
# - Output BAC info to cytoscape compatible *.NA
# - Output RepeatClassification to cytoscape compatbile *.NA
# - Working on adding edge attributes to the output
#
# 06/22/2006
# - Changed format such that all network related outoput
#   is stored in a central dir with *.NA and *.NE related
#   files having same names within the dir
# - Added RepeatName to the set of *.NA files
# - Added Hit Significance *.EA output
# - Added HSP SIF output file
# - Added HSP Frac Identical *.EA output
# 
# 06/26/2006
# - Adding a LaunchCytoscape subfunction
#
# 06/27/2006
# - Installed Cytoscape 2.3, will need to make slight changes
#   to command line to allow for changes in cmd line from
#   v 2.2 to v 2.3. Version 2.3 does a much better job rendering
#   very large networks. V 2.2 was crashing with 2k nodes 18k edges.
# - Rename LaunchCytoscape to LaunchCytoscapeOld
# - Added LaunchCytoscapeNew for version 2.3
# - Added ParseTIGRCode Subfunction
#
# 06/28/2006
# - Convert the drawing of the XY plot to a DrawXY subfunction
# - General code housekeeping (Remove dead code etc.)
# - Added outline-view minor mode to the default perl mode and
#   changed format of my coding to follow rules for this.
# - Added GetRBClass subfunction. This takes the class returned
#   from the BLAST against the RepBase Plants repeat data and
#   returns the name in my Repeat Class ontology.
#
# 06/29/2006
# - Fine tuned GetRBClass Subfunction by looking at output
#   in cytoscape
# - Looked at results of alignment of 'helitron' against the 
#   putative known helitron
#
# 06/30/2006
# - Changed TREP classification to lookup in database table
#   by adding GetTREPClass subfunction
#
# 07/03/2006
# - Working on parsing the Wessler Repeat database information.
#   This will be loaded to a database and then searched by unique
#   name.
#
# 07/05/2006
# - Finished GetTREPClass subfunction
# - Finished GetTREPName subfunction
# - GetWesClass function working for the hits that I have in the
#   ASGR dataset. In the future will want to load these to dbRep for
#   for good results.
# - Adding variable input functions, start with quiet variable
#
# 07/06/2006
# - Added LoadDrosGPI subfunction to parse output from blast against
#   the Drosophila elements and show results as node attribute in 
#   Cytoscape
#
# 12/19/2006
# - Added command line variables to the program.
# - Added the MIPS REDAT database to the parsing.
#
# 01/08/2007
# - Cleaning up program by removing program level variables
#   made redundant by command line options.
#
# 01/16/2007
# - Added the ParseConfigFile subfunction.This allows the user
#   to bypass the length set of variables from the command line
#   and use a configuration file instead.
#
# 01/22/2007
# - Added the TabDeliminted BLAST parser subfunction
# - Currently having problem getting the split to an array
#   function to work from within this subfunction
#  
# 03/03/2007
# - Working on adding a PID out file (PERCENT IDENTITY)
#   and allowing the user to select percent identity to
#   be how the hits are filtered (instead of an e-val cutoff)
#
# 03/04/2007
# -Chaning the name of the SIF output file to be a variable
#  that the use can set. This is the $NetName variable
#  The NetName has been added to the Config file and the cmd
#  line input. The default name is Network
#
# 04/02/2007
# - Updating usage statement.
# - Added create matrix as command line boolean flag -M
# - Cleaning up some options to make more user friendly from
#   the command line.
# 
# 04/03/2007
# - Changed name of jabablast_lite.pl to jabablast
# - Added option to print full usage statement with -h flag
# - Adding option to prodice six edge direction options
#   and coded the ParseTabBLAST subfunction to implement this
#   This should be coded as a separate subfunction to allow
#   me to use the same function for the different local alignment 
#   programs
# - Added command line option to not have a repeat database
# - Added command line option to open or not open Cytoscape 
#
# 04/05/2007
# - Moved help message to a subfunction
# - Put all variable defs under Program Variables header
#   and classifed to make easier to get to
# - Added use strict
# - Set package to RepMiner;
# - Started to move repeat classification subfunctions to RepClass.pl
#   This will require that I send the database handle to to 
#   subfunction for some classifications so I will return to this
#   later. I have more pressing issues.
#
# 05/09/2007
# - Added POD documentation
#
# 05/16/2007
# - Added option to produce GraphViz format. This will need 
#   to use some very standard options since ther command 
#   line is pretty much saturated with the info it can take.
# - Moved GraphViz work to a separate program name for no3
# - Added code to reduce multiple HSPs to a single hit
#   for the ParseTabBlast subfunction
#
# 05/17/2007 
# - Added counter of edges drawn to ParseTabBlast
# - Increasd Java memory for opening Cytoscape from 512
#   to 2048. This will only be reasonable on machines
#   with 2GB of mem, but this will work on jlb10
# - Added LoadTabRepClass as a separate funciton to allow
#   for parsing information from repeat databases that have
#   BLAST reports in the -m8 output format. This will need to
#   be incorporated into a single function where the BLAST
#   alignment format is passed to the subfunction. This option
#   requires that a RepeatDBName be passed at the command line
#   in the form of the -N flag
# - Parsing large graphs to cytoscape
#
# 05/18/2007
# - Changed the ParseRepeatTab blast subfunction so that it
#   no longer uses
# - This works for the SanMiguel database
#
# 05/21/2007
# - Added the ParseTabBLAST2Graph subfunction
# - Currently this creates a undirected graph object
#   and prints clusters relationships to STDOUT
# 
# 05/22/2007
# - Testing ParseTabBLAST2Graph with large datasets 
# - Tested clustering of components by created Cluster
#   out *.NA file and mapping onto Cytoscape network
# - Extending ParseTabBLAST2Graph
#    - Added Fasta file output for each clusters
#    - Added graph stats summary file
# 
# 05/24/2007
# - Added average edge length to summary output
#
# 05/25/2007
# - Working on treating singletons and clusters separately
# 
# 06/27/2007
# - Slight modification to the FASTA output to make sure
#   that it will include the numerical ID that was assigned
#   by Jamie.
#
# 11/19/2007
# - Moved POD doucumentation to the end of the program
# 
# 12/10/2008
# - Rename program to jaba_blast.pl
#
#-----------------------------------------------------------+
# TODO                                                      |
#-----------------------------------------------------------+
# 
# NEEDED
# - Calculate identity of the hit base on HSP info if needed. This
#   should range from 0 to 1.0. 0.95 if 95 Percent identity.
# - Currently this uses a Best Hit to determine the Class and name
#   of the transposable element, I Should also add a Majortiy Rules
#   class for those cases or a best k majority rules
# - Perhaps add a second BLAST Parser for doing the MITEs, This will
#   have a lower stringincy, or could perhaps just do this within
#   the current parser. The mite parser could also be based totally
#   on hmmer evaluation of known models.
# 
# WANTED
# - Classify vector of hit scores etcs. This could be a separate
#   function. (ie. Jenk's optimization to assign colors to a 
#   vector of bit scores obtained from en experiment)
# - Add ability to select the minimum number of clusters an
#   identified cluster must have to be printed to FASTA output
# 
# Additional HSP Level information would be useful
# OTHER REPEAT BLASTS TO PARSE
# RepBase Format
#  -Get Repeat Data from EMBL
#    + Load to DB
#    + Make new FASTA file with class info
#    + Parse EMBL to a Hash for searching by name
# How to work with SanMiguel and Wessler/MAGI data.

