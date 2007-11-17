#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# ASGR BLAST                                                |
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
$QDir = "/home/jestill/projects/asgr/blast/";   # Base dir for query sequences
$DbDir = "/home/jestill/projects/asgr/blast/asgr_db/";
$BlSuf = "-e 0.001 -a 2";                       # Blast suffix
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
@Qry = ( 
	 "rep_TREP2092" # repeat TREP2092
	 );
		  
#-----------------------------+
# BLAST QUERY DATABASES       |
#-----------------------------+
@Db = (
       "dm_gag_pol"           # Drosophical melanogaster gag pol
	);

# Determine the total of BLAST queries that will be run
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
	$TestFile = $DbPath.".phr";
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
	$BlastCmd = "blastall -p blastx -i $QryPath -d $DbPath".
	    " -o $OutPath $BlSuf";
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
