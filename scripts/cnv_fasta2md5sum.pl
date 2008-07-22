#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fasta2md5sum.pl - Convert fasta seq to md5sum string  |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/22/2008                                       |
# UPDATED: 07/22/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a RepMiner formatted FASTA file converts the       |
#  sequence records into a md5sum base on the sequence.     |
#  The outfile is a tab delimted text file giving the       |
#  sequence id and the md5sum.                              |
#                                                           |
# USAGE:                                                    |
#  cnv_fasta2md5sum.pl -i infile.fasta -o outfile.txt       |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#                                                           |
#-----------------------------------------------------------+


#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use strict;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

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
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

#-----------------------------+
# FILE IO                     |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$infile",
			    -format => 'fasta' ) ||
    die "Can not open the intput file:\n$infile\n";

if ($outfile) {
    open (TABOUT, ">$outfile") ||
	die "Can not open output file:\n$outfile\n";
}
else {
    open (TABOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output.\n";
}

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+
my $seq_count = 0;
while (my $seq = $inseq->next_seq) 
{
    $seq_count++;

    my @id_split = split (/\|/ ,$seq->primary_id());


    my $md5 = md5_hex($seq->seq());

    if ($verbose) {
	print STDERR "SEQ: ".$seq->primary_id()."\n" if $verbose;
	print STDERR "\t$md5\n" if $verbose;
	print STDERR "\t".$seq->seq()."\n" if $verbose;
    }

    print TABOUT $id_split[0]."\t".$md5."\n";


} # End of while $seq input loop

    
#-----------------------------+
# CLOSE OUTPUT FILE AND       |
# EXIT PROGRAM                |
#-----------------------------+
close (TABOUT);
exit 0;

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 07/22/2008 
# - Program started 
# - Use the digest md5 program to fetch the md5sum and 
#   send the results to a tab delim text file
# - This program does not currently have any user 
#   accessible help 
