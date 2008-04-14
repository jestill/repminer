#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# extract_blast.pl - Simple extract of m8 BLAST to text     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com
# STARTED: 02/24/2008                                       |
# UPDATED: 02/26/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
# Extract blast -m8 output to sql ready format or other     |
# general use. This also parses the RepMiner header format  |
# and returns the unique identifier number from the fasta   |
# header.                                                   |
#                                                           |
# USAGE:                                                    |
# extract_blast.pl -i infile -o outfile.txt                 |
#                                                           |
#-----------------------------------------------------------+
# FIRST:
# extract_blast.pl my_blast_out.blo > ltr_ltr.txt
# THEN CREATE TABLE TO HOLD THIS DATA:
# mysql> create table tbl_ltr_ltr ( 
#                qry INT, 
#                sub INT, 
#                bit INT, 
#                INDEX (qry), 
#                INDEX (sub));
# It would be better to use DECIMAL for the bitscore
# mysql> create table tbl_ltr_ltr (
#                qry INT,
#                sub INT,
#                bit DECIMAL (7,2),
#                INDEX (qry),
#                INDEX (sub));

# THEN DO THE FOLLOWING AS ROOT:
# cp ltr_ltr.txt /var/lib/mysql/db_maize_ltr_2008/
# THEN DO
# mysql> LOAD DATA INFILE 'ltr_ltr.txt' INTO TABLE tbl_ltr_ltr;
# THIS WILL RESULT IN THE DATA BEING AVAILABLE IN THE MYSQL DATABASE

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # 
use Getopt::Long;              # Get options from command line
use Bio::SearchIO;             # Parse BLAST output

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

my $blast_opt;                 # BLAST option, how to parse blast
                               # Use numbers for now
                               # 1 --> Bitscore of Best HSP from -m8 BLAST
                               # 2 --> Bitscore of Tiled HSPs from -m8 BLAST

# BOOLEANS
my $verbose = 0;
my $line_num = 0;
my $quiet = 0;
my $show_usage = 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0; 

# Index Vals
my $pre_cat_id ="0";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
# Allow infile and outfile to be 

#my $infile = shift;
#my $outfile = shift;

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    "b|blast-opt=s" => \$blast_opt,
                    # ADDITIONAL OPTIONS
                    "q|quiet"       => \$quiet,
                    "verbose"       => \$verbose,
                    # ADDITIONAL INFORMATION
                    "usage"         => \$show_usage,
                    "version"       => \$show_version,
                    "man"           => \$show_man,
                    "h|help"        => \$show_help,);

#-----------------------------+
# CHECK REQUIRED VARIABLES    |
#-----------------------------+
if ( (!$infile) || (!$outfile) || (!$blast_opt) ) {
    print STDERR "\a";
    print STDERR "Infile must be specified at command line\n" if (!$infile);
    print STDERR "Outfile must be specified at command line\n" if (!$outfile);
    print STDERR "Blast parsing option must be specified at command line\n" 
	if (!$blast_opt);
    exit;
}

if ($blast_opt =~ "1" ) {

    print STDERR "Parsing Best HSP\n" if $verbose;

    open (BLASTIN, "<$infile") ||
	die "Can not open $infile\n";
    
    while (<BLASTIN>) {
	
	$line_num++;
	print STDERR "LINE: $line_num\n" if $verbose;
	chomp;
	
	unless (m/^\#.*/) {       # Ignore comment lines, works with -m 9 output
	    
	    my ($qry_id, $sub_id, $pid, $alen, 
		$mismatch, $gap_open, 
		$qry_start,$qry_end, $sub_start, $sub_end,
		$e_value, $bit_score) = split(/\t/);
	    
	    # Trim leading white space from bit score
	    $bit_score =~ s/^\s*(.*?)\s*$/$1/;
	    
	    my ($x_val, $cruft) = split(/\|/, $qry_id);
	    my ($y_val, $more_cruft) = split(/\|/, $sub_id);
	    
	    my $x_int = int($x_val);
	    my $y_int = int($y_val);
	    
	    my $cur_cat_id = $x_val.":".$y_val;
	    
	    # The following just grabs the best hit data
	    unless ( $cur_cat_id =~ $pre_cat_id) {
		print STDOUT "$x_val\t$y_val\t$bit_score\n";
		$pre_cat_id = $cur_cat_id;
	    } # End of grab best hit
	}
    } # END OF WHILE BLAST IN
} # End of if BLAST Option = 0
elsif ( $blast_opt =~ "2" ) {
    
    #-----------------------------+
    # TILED BITSCORE              |
    #-----------------------------+
    my $count_result = 0;
    print "Parsing tiled HSPs\n" if $verbose;
    
    # Open the BLAST report object
    # blasttable for -m8 output
    my $blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
                                           '-file'   => $infile)
        || die "Could not open BLAST input file:\n$infile.\n";

# The following for error reporting    
#    if ($blast_report) { print "Blast report is valid\n"; }
    
    while (my $blast_result = $blast_report->next_result()) {
	
	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    # Do something with the BLAST data
	    my ($x_val, $cruft) = split(/\|/, $blast_result->query_name);
            my ($y_val, $more_cruft) = split(/\|/, $blast_hit->name);

            my $x_int = int($x_val);
            my $y_int = int($y_val);

	    # PRINTING RAW NAMES
	    #print STDOUT $blast_result->query_name."\t";  # Name of the query seq
	    #print STDOUT $blast_hit->name."\t";           # Name of the hit seq
	    # PRINTING PARSED NAMES
	    print STDOUT "$x_val\t$y_val\t";
	    print STDOUT $blast_hit->bits()."\n";         # Tiled Bitscore of the hit
	    
	} # End of next BLAST hit
	
    } # End of next BLAST result
}
elsif ( $blast_opt =~ "4"  ) {
    
    #-----------------------------+
    # TILED PERCENT ID            |
    #-----------------------------+
    
    my $count_result = 0;
    print "Parsing tiled HSPs\n" if $verbose;
    
    # Open the BLAST report object
    # blasttable for -m8 output
    my $blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
                                           '-file'   => $infile)
        || die "Could not open BLAST input file:\n$infile.\n";
    
    if ($blast_report) { print "Blast report is valid\n"; }
    
    while (my $blast_result = $blast_report->next_result()) {
	
	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    # Do something with the BLAST data
	    my ($x_val, $cruft) = split(/\|/, $blast_result->query_name);
            my ($y_val, $more_cruft) = split(/\|/, $blast_hit->name);

            my $x_int = int($x_val);
            my $y_int = int($y_val);

	    # PRINTING RAW NAMES
	    #print STDOUT $blast_result->query_name."\t";  # Name of the query seq
	    #print STDOUT $blast_hit->name."\t";           # Name of the hit seq
	    # PRINTING PARSED NAMES
	    print STDOUT "$x_val\t$y_val\t";
	    print STDOUT $blast_hit->bits()."\n";         # Tiled Bitscore of the hit
	    
	} # End of next BLAST hit

    } # End of next BLAST result
} else {
    print "\a";
    print "A valid BLAST parsing option was not passed at the command line\n";
}

# END OF PROGRAM
exit;

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 
#
# 03/09/2008
# - Added Getopt long 
# - Added infile and outfile as reuired command line options
# - Added option to select BLAST parsing options
#    1 --> Bitscore of Best HSP from -m8 BLAST output
#    2 --> Bitscore of Tiled HSP from -m8 BLAST output 
#          This option uses the BioPerl BLAST parsing
#    3 --> PID of Best HSPs
#    4 --> PID of Tiled HSPs
# 
# TO DO:
# Add option for PID instead of BIT SCORE
