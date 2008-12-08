#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_add_num.pl - Add integer prefix to fasta headers    |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/12/2006                                       |
# UPDATED: 12/18/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Add a number prefix to the id line from a fasta file     |
#  to allow for easy parsing of fasta sequences. The process|
#  also provides each sequence with a unique incremental    |
#  id that can be used for reference.                       |
#                                                           |
# USAGE:                                                    |
#  fata_add_num.pl -i infilepath -o outfilepath             |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#  -DBI                                                     |
#  -MySQL                                                   |
#                                                           |
#-----------------------------------------------------------+
#
# NOTE:
# It may just as easy to do this without any references to
# the bioperl modules. Just read through the input file and
# append the incremental number to the beginning of the 
# sequence. An advantage to using the Bioperl path is that
# a different sequence input format could be used.
#

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use DBI;                       # Allows connection to MySQL database
use Bio::SeqIO;                # Allows for treatment of seqs as objects
use Text::Wrap;                # Allows printing of wrapped text
#use Getopt::Std;               # Allows options flags at command line
use Getopt::Long;              # Get options from command line

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
# USAGE STATEMENT
#my $Usage = "fasta_add_num.pl -i InFilePath -o OutFilePath\n";
my $seq_in_format = "fasta";   # The input sequence format
my $seq_out_format = "fasta";  # The output sequence format
my $seq_num = 0;                # Var to keep track of nummber of seqs 
my $seq_unique_id;               # Unique ID attributed to the seq read record
my $new_id;

# BOOLEANS
my $verbose = 0;
my $line_num = 0;
my $quiet = 0;
my $show_usage = 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0; 
my $blast_opt = 0;

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
                    # ADDITIONAL OPTIONS
		    "in-format"     => \$seq_in_format,
		    "out-format"    => \$seq_out_format,
                    "q|quiet"       => \$quiet,
                    "verbose"       => \$verbose,
                    # ADDITIONAL INFORMATION
                    "usage"         => \$show_usage,
                    "version"       => \$show_version,
                    "man"           => \$show_man,
                    "h|help"        => \$show_help,);

#-----------------------------+
# FASTA FILE IO               |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$infile",
			    -format => $seq_in_format ) ||
    die "ERROR Can not open infile:\n $infile\n";

my $outAll = Bio::SeqIO->new(-file   => ">$outfile",
			     -format => $seq_out_format ) ||
    die "ERROR Can not open outfile:\n $outfile\n";

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+

while (my $seq = $inseq->next_seq)  {
    # This process may not work with all FASTA files, some
    # of the information will be left out of the new FASTA file
    # that is produced. It would be good to just capture the entire
    # FASTA header and append the incremented number.
    $seq_num++;
    $seq_unique_id = $seq->primary_id;
    $new_id = $seq_num."|".$seq_unique_id ;

    # Restet the seqid to the new id
    $seq->primary_id( $new_id );
    $seq->display_id( $new_id );

    #-----------------------------+
    # PRINT THE NEW ID            |
    #-----------------------------+
    print STDERR "\tNEW ID: ".$new_id."\n" if $verbose;

    #-----------------------------+
    # WRITE SEQUENCE RECORDS OUT  |
    # TO A NEW FILE               |
    #-----------------------------+
    $outAll->write_seq($seq);
	

} # END OF THE FOR EVERY SEQUENCE RECORD 

exit;



#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

__END__;


=head1 NAME

cnv_blast2sim.pl - Convert BLAST report to similarity file

=head1 VERSION

This documentation refers to cnv_blast2sim.pl version $Rev: 95 $

=head1 SYNOPSIS

=head2 Usage

    cnv_blast2sim.pl -i BlastOutput.bln -o SimFile.txt

=head2 Required Arguments

    -i, --infile    # Path to the input file to parse
    -o, --outfile   # Path to the output similarity file

=head1 DESCRIPTION

Extract blast -m8 or -m9 output to a simple three column
similarity file. This output file is a tab delimited
text file that can be easily imported into a SQL database
or used as an input file for clustering.
As part of the process this script also parses the
RepMiner header format and returns the unique identifier
number from the fasta header.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the BLAST file to be parsed. If no infile argument is used
the program will attempt to parse from STDIN.

=item -o,--outfile

Path of the output file. If not outfile arugment is used, the program
will send data to standard output.

=back

=head1 OPTIONS

=over 2

=item --in-format

The file format of the sequence to prepend a sequence ID to. Be default the
program will expect sequenced to be in the fasta format, however any
valid sequence format accpeted by bioperl may be used.

=item --out-format

The file format of the output sequence file. Be default, fasta format files
will be produced, but any file format acceptable by bioperl may be used.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands.

=back

=head1 EXAMPLES

Examples here

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item Program does not appear to be doing anything

It is possible that the program is waiting for input from STDIN.
Run with the --verbose mode.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This script does not make use of external configuration files
or any options in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * NCBI blastall

The latest version of the NCBI blastall program can be downloaded from:
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST

=back

=head2 Required Perl Modules

=over

=item * Bio::SearchIO

This module is part of bioperl and is required to parse BLAST
output in a format that is tiled across all HSPs.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the RepMiner
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

=over

=item * Limited to m8 or m9 BLAST format

This script is designed to be a lean and fast parser of the 
similarity information from BLAST. It is therefore limited
to using the simple m8 or m9 BLAST alignment format.

=back

=head1 SEE ALSO

The cnv_blast2sim.pl program is part of the repminer package of 
repeat element annotation programs.
See the DAWG-PAWS web page 
( http://repminer.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/repminer/ )
for additional information about this package.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 07/12/2006

UPDATED: 12/07/2008

VERSION: $Rev$

=cut


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/12/2006
# - Started program
#
# 12/18/2006
# - Changed Name to FastaAddNum.pl
# - Added input variable flags
#
# 12/07/2008
# - Adding POD documentation
# - Adding help subfunction
