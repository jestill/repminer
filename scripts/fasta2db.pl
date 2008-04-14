#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# LIRIO ADD NUMBER                                          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/12/2006                                       |
# UPDATED: 07/12/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Add a number prefix to the id line from a fasta file     |
#  to allow for easy parsing of fasta sequences.            |
#                                                           |
# USAGE:                                                    |
#  ParseASGR.pl                                             |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#  -DBI                                                     |
#  -MySQL                                                   |
#                                                           |
#-----------------------------------------------------------+

print "The program has started\n";

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
# FASTA FILE TO PARSE
my $quiet;
#my $SeqFile = "/home/jestill/projects/liriodendron/orig/Lirio_data.fasta";
my $AllOutPath = "/home/jestill/projects/liriodendron/Lirio_data_num";
my $SeqFormat = "fasta";       # The input sequence format
my $SeqNum = 0;                # Var to keep track of nummber of seqs 
my $SeqUniqueId;               # Unique ID attributed to the seq read record
my $NewId;
#my $DbUserName = "jestill";    # Username for MySQL Database
#my $DbName = "dbSanMiguel";        # DB Name for MySQL Database
my $DbUserPassword;            # Sets scope for the password
#my $SeqTbl = "tblSeqData";     # Name of the data to hold the seq data

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
# Declare the Options hash to hold the options from the command line
# Can provide the option for a variable input format here
# q flag is for running in quiet mode if desired.
my %Options;
getopts('i:d:u:o:p:t:q', \%Options);
my $Usage = "Fasta2DB -i InputFilePath -d dbName -u dbUserName [q]";
my $SeqFile = $Options{i} || 
    die "You must provide an input file path\n$Usage\n";
my $AllOutPath = $Options{o} || 
    die "You must provide an output path\n$AllOutPath\n";
my $DbName  = $Options{d} || 
    die "You must the database name\n$Usage\n";
my $DbUserName = $Options{u} || 
    die "You must provide a user name for the database connection.\n$Usage\n";
# THE FOLLOWING VARIABLES ARE NOT REQUIRED AT THE COMMAND LINE
my $SeqTbl = $Options{t} ||
    "tblSeqData"; # Default SeqTable name is tblSeqData
my $quiet = $Options{q};
$DbUserPassword = $Options{p};

#-----------------------------+
# GET THE USER PASSWORD FROM  |
# THE COMMAND LINE            |
#-----------------------------+
# If the user did not give the password at the command
# line then get it here.
unless ($DbUserPassword)
{
    print "\nPassword for $DbUserName\n";
    system('stty', '-echo') == 0 or die "can't turn off echo: $?";
    $DbUserPassword = <STDIN>;
    system('stty', 'echo') == 0 or die "can't turn on echo: $?";
    chomp $DbUserPassword;
}

#-----------------------------+
# FILE/DB IO                  |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$SeqFile",
			    -format => $SeqFormat );

my $outAll = Bio::SeqIO->new(-file   => ">$AllOutPath",
			     -format => "fasta");

my $dbh = DBI->connect("DBI:mysql:database=$DbName;host=localhost",
		       $DbUserName, $DbUserPassword,
		       {'RaiseError' => 1});

# Create the database table to hold the information
&CreateDBTable ( $SeqTbl );
#exit;

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+

while (my $seq = $inseq->next_seq) 
{
    $SeqNum++;
    $SeqUniqueId = $seq->primary_id;
    $NewID = $SeqNum."|".$SeqUniqueId ;
    $Species = "Unknown";
    # THE FOLLOWING MAY NOT BE NEEDED 
    # WITH THE LIRIODENDRON DATA
    #$seq->desc('');

    $seq->primary_id( $NewID );
    $seq->display_id( $NewID );

    print "\tNEW ID: ".$NewID."\n";
    
    #-----------------------------+
    # WRITE SEQUENCE RECORDS OUT  |
    # TO A NEW FILE               |
    #-----------------------------+
    #$outAll->write_seq($seq);
	
    
    #------------------------------+
    # UPLOAD INFO TO THE DATABASE  |
    #------------------------------+
    $UploadData = "INSERT IGNORE INTO ".$SeqTbl.
	" ( species, id, seq, seq_len)".
	" VALUES (".
	" '".$Species."',".
	" '".$SeqUniqueId."',".
	" '".$seq->seq."',".
	" '".$seq->length."'".
	" )";

    $dbh->do($UploadData);

} # END OF THE FOR EVERY SEQUENCE RECORD 

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


#-----------------------------+
# CREATE THE DATABASE TABLE   |
#-----------------------------+
sub CreateDBTable
{
    my $TableName = $_[0];   
    
    # DROP THE TABLE it already exists
    if (&does_table_exist( $dbh, $TableName ))
    {
	$dbh->do("DROP TABLE ".$TableName);
    }

    $qryCreateTable = "CREATE TABLE ".$TableName." (".
	" rownum INT(10) NOT NULL AUTO_INCREMENT,".
	" species char(20),". # Species the seq is derived from
	" id char(255),".      # Long unique ID for the Seq
	" seq MEDIUMTEXT,".   # Sequence for the read
	" seq_len char(5),".  # Length of the high quality seq
	" KEY (rownum))";
    $dbh->do($qryCreateTable);

}

sub does_table_exist
{
    my ($dbh,$whichtable) = @_;
    my ($table,@alltables,$found);
    @alltables = $dbh->tables();
    $found = 0;
    foreach $table (@alltables) {
	#$found=1 if ($table eq $whichtable)
	# The above quit working sometime around 5/15/2005 and was change to the following
	$found=1 if ($table eq "`".$whichtable."`");
    }
    # return true if table was found, false if table was not found
    return $found;
}

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/12/2006
# - Largely modified from ParseASGR3.pl
#
# 04/05/2007
# - Added the ability to specify the password at the command
#   line using the -p flag
#
# 11/19/2007 
# - Cleaned out old head1 and changed name from
#   Fasta2DB.pl to fasta2db.pl
