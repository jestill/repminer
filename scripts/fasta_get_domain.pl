#!/usr/bin/perl -w                                                                          
#-----------------------------------------------------------+                               
#                                                           |                               
# fasta_orient.pl - Put seqs in same orientation
#                                                           |                               
#-----------------------------------------------------------+                               
#                                                           |                               
#  AUTHOR: James C. Estill                                  |                               
# CONTACT: JamesEstill_@_gmail.com                          |                               
# STARTED: 03/24/2008                                       |
# UPDATED: 03/24/2008                                       |
#                                                           |                               
# DESCRIPTION:                                              |                               
#  Given a multi-fasta file compare the sequence to a       |
#  protein database and use the direction of the hit        |
#  to place all sequences in the input seq file in the      |
#  same orientation with respect to the query database.     | 
#  This is designed to convert full length LTR Retros       |
#  to their reverse complement for annotation them          |
#  all in the plus orientation with respect to their        |
#  protein coding regions.                                  |
#                                                           |
# USAGE:                                                    |
#  fasta_orient.pl -i myseqs.fasta -d database              |
#                  -o oriented_seqs.fasta                   |
#                                                           |
# VERSION: $Rev$                                            |                               
#                                                           |                               
# LICENSE:                                                  |                               
#  GNU General Public License, Version 3                    |                               
#  http://www.gnu.org/licenses/gpl.html                     |                               
#                                                           |                               
#-----------------------------------------------------------+

package REPMINER;


#-----------------------------+                                                             
# INCLUDES                    |                                                             
#-----------------------------+                                                             
use strict;
use Getopt::Long;
use Bio::SeqIO;
use IPC::Open2;                # For interaction with prss34                                

#-----------------------------+                                                             
# PROGRAM VARIABLES           |                                                             
#-----------------------------+                                                             
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+                                                             
# VARIABLE SCOPE              |                                                             
#-----------------------------+                                                             
my $infile;                    # Input fasta file
my $outfile;                   # Output fasta file
my $logfile;                   # Log file
my $database;                  # Database of aa seqs to blastx against
my $do_revcom = 0;             # Do the reverse complement

# Valus with default values
#my $e_val = "0.1";

# Booleans                                                                                  
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+                                                             
# COMMAND LINE OPTIONS        |                                                             
#-----------------------------+                                                             
my $ok = GetOptions(# REQUIRED OPTIONS                                                      
                    "i|infile=s"    => \$infile,
                    "o|outdir=s"    => \$outfile,
		    "d|database=s"  => \$database,
                    # ADDITIONAL OPTIONS
#		    "e|e_val=s"     => \$e_val,
		    "l|log=s"       => \$logfile,
                    "q|quiet"       => \$quiet,
                    "verbose"       => \$verbose,
                    # ADDITIONAL INFORMATION
                    "usage"         => \$show_usage,
                    "version"       => \$show_version,
                    "man"           => \$show_man,
                    "h|help"        => \$show_help,);

#-----------------------------+
# OPEN FILES AND HANDLES
#-----------------------------+

# INPUT DNA SEQUENCE FILE
my $seqs_in = Bio::SeqIO->new( -file => "<$infile",
                              -format => 'fasta') ||
    die "Can not open input sequence file one:\n$infile\n";

# LOG FILE
if ($logfile) {
    open (LOGOUT, ">$logfile") ||
	die "Can not open a log file.\n";
    
}

#-----------------------------+
# REORIENT SEQS AS NEEDED     |
#-----------------------------+

while (my $seq = $seqs_in->next_seq) {

    print STDERR "Processing ".$seq->primary_id()."..";
    print LOGOUT "Processing ".$seq->primary_id().".." if ($logfile);
    
    #my $process = open2 (\*RFH, \*WFH, 'blastall -p blastx -d '.$database.' -m 8 -b 1 -e '.$e_val) ||
    my $process = open2 (\*RFH, \*WFH, 'blastall -p blastx -d '.$database.' -m 8 -b 1') ||
	die "Could not do blast.\n";

    # Write the BLAST query string to the Write file handle and close it
    print WFH ">".$seq->primary_id()."\n".$seq->seq()."\n";
    close (WFH);
   
    my $hit_num = 0;
    
    while (<RFH>) {
	$hit_num++;

	# Print BLAST output to STDERR if running in verbose
	print STDERR "\n".$_ if $verbose;

	my ($qry_id, $sub_id, $pid, $len,
	    $mismatch, $gap_open,
	    $qry_start,$qry_end, $sub_start, $sub_end,
	    $e_val, $bit_score) = split(/\t/);
	
#	if ($sub_start > $sub_end) {
	if ($qry_start > $qry_end) {
	    $do_revcom = 1;
	} 
	else {
	    $do_revcom = 0;
	}

    }

    close (RFH);

    # Write sequences to the outfile
    if ($hit_num == 0 ) {
	# If not blast HITS, print as is and write error message to STDERR

	# Err message
	print STDERR "\tNO BLAST HIT\n";
	print LOGOUT "\tNO BLAST HIT\n" if ($logfile);

	# Fasta output
	print STDOUT ">".$seq->primary_id()."\n";
	print STDOUT $seq->seq."\n";
    }
    else {
	if ($do_revcom) {
	    # Err message
	    print STDERR "\tDoing Revcom\n";
	    print LOGOUT "\tDoing Revcom\n" if ($logfile);

	    # Fasta output
	    print STDOUT ">".$seq->primary_id()."\n";
	    print STDOUT $seq->revcom."\n";
	}
	else {
	    # Err message
	    print STDERR "\tNo Revcom\n";
	    print LOGOUT "\tNo Revcom\n" if ($logfile);

	    # Fasta output
	    print STDOUT ">".$seq->primary_id()."\n";
	    print STDOUT $seq->seq."\n";
	}
    }

}

close (STDOUT);
close (LOGOUT) if ($logfile);

#-----------------------------------------------------------+
# HISTORY 
#-----------------------------------------------------------+
# 03/24/2008
# - Main body of program writtern
# - Assumes that the qry sequence will always be flipped
#   with respect to the protein database
# 03/25/2008
# - Modifying standard error file to be easier to parse
#   This will let me easily see which seqs could not
#   be oriented.
# - Adding option to print STDERR info to a logfile
#   if passed at the command line
