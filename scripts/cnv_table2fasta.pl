#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fasta2md5sum.pl - Convert fasta to text file for db   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/24/2008                                       |
# UPDATED: 07/24/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert a three columng table file into a fasta format   |
#  sequence file. Input format is                           |
#   col 1: seq id unique integer                            |
#   col 2: longer seq id name                               |
#   col 3: sequence string                                  |
#                                                           |
# USAGE:                                                    |
#  cnv_table2fasta.pl -i infile.txt -o outfile.fasta        |
#                                                           |
# REQUIREMENTS:                                             |
#  -Uses basic PERL commands. No special modules needed     |
#                                                           |
# OUTPUT:                                                   |
#  Output is a fasta format sequence file:                  |
#  > seq_id_num | long_seq_id                               |
#  sequence_string                                          |
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
if ($infile) {
    open (TABIN, "<$infile") ||
	die "Can not open the input file:\t$infile\n";
}
else {
    open (TABIN, "<$STDIN") ||
	die "Can not open STDIN for input.\n";
}

if ($outfile) {
    open (FASTOUT, ">$outfile") ||
	die "Can not open output file:\n$outfile\n";
}
else {
    open (FASTOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output.\n";
}

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+
while (<TABIN>) {
    chomp;
    my @string_split = split (/\t/, $_);
    print STDERR "Processing: ".$string_split[0]."\n" if $verbose;
    
    print FASTOUT ">".$string_split[0]."|".$string_split[1]."\n".
	$string_split[2]."\n";
    
}

#-----------------------------+
# CLOSE OUTPUT FILE AND       |
# EXIT PROGRAM                |
#-----------------------------+
close (FASTOUT);
exit 0;

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 07/24/2008
# - Program started
