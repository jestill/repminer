#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fasta2md5sum.pl - Convert fasta to text file for db   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/24/2008                                       |
# UPDATED: 10/28/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a RepMiner formatted FASTA file converts the       |
#  sequence records into a tab delimited text file          |
#  suitable for upload to a database.                       |
#  The outfile is a tab delimted text file giving the       |
#  sequence id, and the md5sum.                             |
#  This is designed for processing output of the program    |
#  cnv_ltrstruc2ann.pl or cnv_ltrstruc2fasta.pl.            |
#  The resulting file is suitable for easy upload to        |
#  a databse.                                               |
#                                                           |
# USAGE:                                                    |
#  cnv_fasta2tabinfo.pl -i infile.fasta -o outfile.txt      |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#  -Digest:MD5                                              |
#                                                           |
# OUTPUT:                                                   |
#  Output is four part tab delimited text file:             |
#    col 1: seq id                                          |
#    col 2: feature name                                    |
#    col 3: sequnce string                                  |
#    col 4: md5sum                                          |
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
while (my $seq = $inseq->next_seq) {
    $seq_count++;

    my @id_split = split (/\|/ ,$seq->primary_id());

    my $md5 = md5_hex($seq->seq());

    if ($verbose) {
	print STDERR "Processing: ".$seq->primary_id()."\n";
    }
    

    # The following will iterate across idsplit if there
    # is more than one feature in the split list
    # otherwise will keep the primary ID as the seq_info
    my $seq_header_info;
    my $num_id_parts = @id_split;
    
    if ($num_id_parts < 2 ) {
	$seq_header_info = $seq->primary_id();
    } else {
	foreach my $id_part(@id_split) {
	    $seq_header_info = $seq_header_info.$id_part."\t";
	}
    }

    print TABOUT $seq_header_info.
	$seq->length()."\t".
	$md5."\t".
	$seq->seq()."\t".
	"\n";


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
# 07/24/2008
# - Program started
# 05/20/2009
# - Added sequence length
