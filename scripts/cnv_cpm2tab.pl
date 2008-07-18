#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_cpm2tab.pl - Convert MCL CPM file to tab delim txt    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/17/2008                                       |
# UPDATED: 07/17/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert the output from the MCL program                  |
#   Assumes clm info --node-all-measure                     |
#  for now to get a simple parse of this data.              | 
#  Can use -i and -o flags, or STDIN and STDOUT.            |
#                                                           |
# EXAMPLES:                                                 |
# cnv_cpm2tab.pl -i infile.cpm -o outfile.cpm.txt           |
#                                                           |
#-----------------------------------------------------------+
# TODO:

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # 
use Getopt::Long;              # Get options from command line
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path


#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev:$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# BOOLEANS
my $verbose = 0;
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
$infile = shift;
$outfile = shift;

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
                    # ADDITIONAL OPTIONS
                    "q|quiet"       => \$quiet,
                    "verbose"       => \$verbose,
                    # ADDITIONAL INFORMATION
                    "usage"         => \$show_usage,
                    "version"       => \$show_version,
                    "man"           => \$show_man,
                    "h|help"        => \$show_help,);


#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\n$0\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED VARIABLES    |
#-----------------------------+
# Currrently not variables are required.

#-----------------------------------------------------------+
# MAIN BODY                                                 |
#-----------------------------------------------------------+

# Open infile if one is passed, otherwise use STDIN
if ($infile) {
    open (CLMIN, "<$infile") ||
	die "Can not open infile:\n$infile";
}
else {
    open (CLMIN, "<&STDIN") ||
	die "Can not open STDIN for input\n"
}

# Print to outfile if one specified, otherwise to STDOUT
if ($outfile) {
    open (TABOUT, ">$outfile") ||
	die "Can not open output file\n$outfile\n";
} 
else {
    open (TABOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output\n";
}

#-----------------------------------------------------------+
# PARSE CLM TO TAB DELIMITED OUTPUT                         |
#-----------------------------------------------------------+
my $line_num = 0;
while (<CLMIN>) {
    chomp;
    $line_num++;

    # Print line
    my @clm_parts = split;
    my $nm = substr ($clm_parts[0], 3);   # Name of the graph file
    my $ni = substr ($clm_parts[1], 3);   # Node index id
    my $ci = substr ($clm_parts[2], 3);   # Cluster index id
    my $nn = substr ($clm_parts[3], 3);   # Number of neighbors of the node
    my $nc = substr ($clm_parts[4], 3);   # Cluster size
    my $ef = substr ($clm_parts[5], 3);   # Efficiency
    my $em = substr ($clm_parts[6], 3);   # Efficiency Max
    my $mf = substr ($clm_parts[7], 3);   # Mass fraction
    my $ma;   # --
    my $xn = substr ($clm_parts[9], 3);    # Num of neighbors not in the cluster
    my $xc = substr ($clm_parts[10], 3);   # Num cluster nodes not in neigh list
    my $ns;
    my $ti;
    my $to;
    my $al;

    #print "$_\n";
    #print "\t$ni\n";
    #print "\t$ci\n";
    #print "\t$ef\n"; # 
    #print "\t$mf\n";

    # AS A SIMPLIFIED TAB DELIMITED FILE OF MOST IMPORTANT PARTS
    print TABOUT "$ni\t$ci\t$ef\t$mf\n";

#    if ($line_num == 5) {
#	print STDERR "\nDebug exit\n\n";
#	exit 1
#    };
}

# END OF PROGRAM
close (CLMIN);
close (TABOUT);
exit 0;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

__END__;


=head1 NAME

cnv_blast2sim.pl - Convert BLAST report to similarity file

=head1 VERSION

This documentation refers to cnv_blast2sim.pl version $Rev: 46 $

=head1 SYNOPSIS

=head2 Usage

    cnv_blast2sim.pl -i BlastOutput.bln -o SimFile.txt

=head2 Commonly Used Arguments

    -i, --infile    # Path to the input file to parse
    -o, --outfile   # Path to the output similarity file
    -b              # Blast data to use for similairty metric 

=head1 DESCRIPTION

Extract blast -m8 or -m9 output to a simple three column
similarity file. This output file is a tab delimited
text file that can be easily imported into a SQL database
or used as an import file for clustering.
As part of the process this script also parses the
RepMiner header format and returns the unique identifier
number from the fasta header.

=head1 COMMONLY USED ARGUMENTS

=over 2

=item -i,--infile

Path of the BLAST file to be parsed. If no infile argument is used
the program will attempt to parse from STDIN.

=item -o,--outfile

Path of the output file. If not outfile arugment is used, the program
will send data to standard output.

=item -b, --blast-opt

Number indicating the option used to extract a similarity
score from a BLAST result. Valid arguments are:

=over

=item 1.

Use the BITSCORE of the best HSP.

=item 2.

Use a tiled BITSCORE value. This is the default option.

=item 3.

Use the significance value of the best HSP. 
THIS OPTION IS NOT CURRENTLY IMPLEMENTED.

=item 4.

Use a tiled significance value

=back

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

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

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

=over

=item * Limited to m8 or m9 BLAST format

This script is designed to be a lean a fast parser of the 
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

STARTED: 02/24/2008

UPDATED: 07/15/2008

VERSION: $Rev: 46 $

=cut

#-----------------------------------------------------------+
# CHANGELOG                                                 |
#-----------------------------------------------------------+
#
# 03/09/2008
# - Added Getopt long 
# - Added infile and outfile as reuired command line options
# - Added option to select BLAST parsing options
#    1 --> Bitscore of Best HSP from -m8 BLAST output
#    2 --> Bitscore of Tiled HSP from -m8 BLAST output 
#          This option uses the BioPerl BLAST parsing
#    3 --> Significance of the best scoring HSP
#    4 --> Significance of Tiled HSPs
# 
# 07/15/2008
# - Renamed from extract_blast.pl to cnv_blast2sim.pl
#   this name reflects the fact this this converts native
#   m8 BLAST output to a simple three column similarity
#   file.
# - Added the option to print output to STDOUT when the
#   the infile argument is not given
# - Added the option for input from STDIN

