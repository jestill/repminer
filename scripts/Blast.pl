#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# LIRIO BLAST                                               |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 05/10/2006                                       |
# UPDATED: 05/10/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run NCBI blast for the asgr sequences against a set of   |
#  repeast databases of interest.                           |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

# Can run this with a test flag to just run through
# the set to make sure that all of the databases and input
# files exist.

print "The program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Allows to get options from the command line
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.
#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
# DIR FOR ALL ASGR BLAST ANALYSIS QUERY SEQUENCES
#$QDir = "/home/jestill/projects/liriodendron/";   # Base dir for query sequences
$QDir = "/home/jestill/projects/RepMiner/rbs/";   # Base dir for query sequences


#$QDir = "/home/jestill/projects/RepMiner/wgs/";   # Base dir for query sequences

$DbDir = "/home/jestill/blast/repeats/";        # Base dir for db sequences
$BlSuf = "-e 0.00001 -a 2";                # Blast suffix
#$BlProg = "tblastx";                            # Blast program to use 
$BlProg = "blastn";                            # Blast program to use 
$ProcNum = 0;                                   # Process number starts at zero
my ( $BlastCmd, $QryPath, $DbPath, $OutPath );  # Declare scope for varaiables used later

#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
my %Options;                  # Options hash to hold the options from the command line
getopts('t', \%Options);      # Get the options from the command line
my $test;
$test = $Options{t};

#-----------------------------+
# BASE NAMES OF FILES TO USE  | 
# AS QUERY SETS FOR BLAST     |
#-----------------------------+
@Qry = ( #"rbs"
	 "maize_wgs"
	 );
		  
#-----------------------------+
# BLAST QUERY DATABASES       |
#-----------------------------+
# TIGR Repeats from : 
# ftp://ftp.tigr.org/pub/data/TIGR_Plant_Repeats/
# These are the repeat database from different sources
@Db = ( #"gram_rep",            # Gramineae v3.1 repeat database from TIGR
#	"os_rep",              # Oryza sativa repeat database from TIGR
#	#"RB_ath",              # A. thaliana repeat database from RepBase
#	#"RB_ory",              # Oryza sativa repeat database from RepBase
#	"RB_pln",              # All plants repeat database from RepBase
#	"TREP_8",              # TREP v.8 database
#	#"SanMiguel",           # SanMiguel Repeat Database
#	"Wessler",             # MAGI Wessler Repeat Database
#	"PM_TIR",              #
#	"TIGRfab",             # TIGR Fabaceae Repeat Database
#	"TIGRbras",            # TIGR Brassiceae Repeat Database
#	"TIGRsol",              # TIGR Solanaceae Repeat Database
#	"vector",
#	"zm_rep",              # Zea mays repeat database from TIGR
	"SanMiguel"           # SanMiguel Repeat Database

	);

my $LenQry = @Qry;
my $LenDb =  @Db;
my $NumProc = $LenQry * $LenDb;

for $IndQry (@Qry)
{

    for $IndDb (@Db)
    {

	$ProcNum++; # Increment the process number
	print "BLAST Process ".$ProcNum." of ".$NumProc."\n";

	$QryPath = $QDir.$IndQry."/".$IndQry.".fasta";
	$DbPath = $DbDir.$IndDb;
	$TestFile = $DbPath.".nhr";
	$OutPath = $QDir.$IndQry."/".$IndQry."_".$IndDb.".blo";

	#-----------------------------+
	# DOES THE BLAST DB EXIST     |
	#-----------------------------+
	if (-e $TestFile ) 
	{
	    #print "DB: $IndDb exists\n";
	}
	else 
	{die "Can not find database:\n$IndDb\n"; }
	
	#-----------------------------+
	# DOES THE BLAST QRY EXIST    |
	#-----------------------------+
	if (-e $QryPath)
	{
	    #print "QRY: $QryPath exists\n";
	}
	else
	{die "Can not find qry file:\n$QryPath\n";}

	#------------------------------+
	# PRINT THE BLAST COMMAND      |
	#------------------------------+
	# Added the -m 8 for tab delim output
#	$BlastCmd = "blastall -p ".$BlProg." -i $QryPath -d $DbPath -o $OutPath $BlSuf";
	$BlastCmd = "blastall -p ".$BlProg." -i $QryPath -d $DbPath -o $OutPath $BlSuf";
	print wrap("\t", "\t", $BlastCmd );
	print "\n";
	#------------------------------+
	# RUN THE BLAST COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
	    system ($BlastCmd);
	}
	


    } # End of for each database loop

} # End of for each query seq loop


exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 05/10/2006
# - Program was started to create an easy way to blast all
#   of the query datasets of interest by all of the repeat
#   databases of interest. The BLAST reports should then
#   be parsed and uploaded to the dbASGR database.
