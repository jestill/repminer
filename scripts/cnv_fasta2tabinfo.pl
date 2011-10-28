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
    

    # TEMP MODIFICATION TO GET THE LTR_STRUC ID I AM USING


    # ltr struc output is difficult to parse to to the
    # switch of delimiters, sometimes they are present
    # and sometimes they are not.
#    my @split_id = split(/\_/, $id_split[0]);

#    my $len_split_ltr_id = @split_ltr_id;
#    print STDERR "$len_split_ltr_id\n";
#    print "\t".$id_split[0]."\n";
    
#    my $ltr_struc_id = "NULL";
#    if ($id_split[0] =~ m/(.*)_rprt.txt_(.*)/) {
#	#print STDERR $1."\n";
#	$ltr_struc_id = $1;
#    } else {
#	$ltr_struc_id = $id_split[0];
#    }

    
#    print STDERR "Processing: ". $seq_count." : ".$id_split[3]."\n";

#    print TABOUT "$ltr_struc_id\t".
#	$id_split[3]."\t".
#	$id_split[1]."\t".
#	$id_split[4]."\t".
#	$seq->length()."\t".
#	$seq->seq()."\t".
#	$md5."\n";

    # Modified to the following as a kludge 10/28/2011
    # mut need to iterate across idsplit
    print TABOUT $id_split[0]."\t".
	$id_split[1]."\t".
	$seq->length()."\t".
	$md5.
	$seq->seq()."\t".
	"\n";

#    my @src_info = split(/\./, $id_split[3]);
#    my $src_version = $src_info[1];

#    # PRINT TAB DELIM OUTPUT OF NCBI SEQS
#    print TABOUT "gi\t".
#	$id_split[1]."\t".
#	"gb\t".
#	$src_version."\t".
#	$id_split[3]."\t".       
#	$seq->desc()."\t".
#	$seq->seq()."\n";

#    print TABOUT $id_split[0]."\t".
#	$id_split[1]."\t".
#	$seq->seq()."\t".
#	$md5."\n";


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
