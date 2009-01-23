#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_sim2norm.pl - Normalizes similarity values of matrix  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/23/2009                                       |
# UPDATED: 01/23/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  progname -i infile -o outfile.fasta                      |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
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
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $include_val = 0;              # include the non-normalize value in output
#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"     => \$infile,
                    "o|outfile=s"    => \$outfile,
		    # ADDITIONAL OPTIONS
		    "include-val"    => \$include_val,
		    "q|quiet"        => \$quiet,
		    "verbose"        => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"          => \$show_usage,
		    "test"           => \$do_test,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,);

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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# OPEN FILE HANDLES           |
#-----------------------------+
# DEFAULT TO STDIN
if ($infile) {
    open (SIMIN, "<$infile") ||
	die "Can not open input file:\n$infile\n";
} 
else {
    print STDERR "Expecting input from STDIN\n";
    open (SIMIN, "<&STDIN") ||
	die "Can not access standard input\n";
}

# DEFAULT TO STDOUT
if ($outfile) {
    open (SIMOUT, ">$outfile") ||
	die "Can not open output file:\n$outfile\n";
}
else {
    open (SIMOUT, ">&STDOUT") ||
	die "Can not access standard output\n";
}

#-----------------------------+
# LOAD SIMILARITY VALS TO ARY |
#-----------------------------+
my (@sim_vals);
while (<SIMIN>) {
    chomp;
    my ($row,$col,$val) = split;
    print STDERR "Processing $_\n" if $verbose;
    $sim_vals[$row][$col] = $val;
}

#-----------------------------+
# NORMALIZE SIM VALS          |
#-----------------------------+
# Read through values, use i j from [0] [1] in input 
# this will make this pseudo sparse matrix
# Restricts program in that it must use i j val method
# But this is a general program that convert any set of similarity values
# to normalized values

# Reset file handle
seek SIMIN,0,0;
while (<SIMIN>) {
    my ($i,$j,$val) = split;

    my $norm_sim = $val/get_min( $sim_vals[$i][$i], $sim_vals[$j][$j] );
    
    if ($include_val) {
	print SIMOUT "$i\t$j\t$norm_sim\t$val\n";
    }
    else {
	print SIMOUT "$i\t$j\t$norm_sim\n";
    }

}

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub get_min { 

    # If the values are the same this returns the first value
    my ($x, $y) = @_;
    my $ans;

    if ($x <= $y) {
	$ans = $x;
    }
    else {
	$ans = $y;
    }

    return $ans;

}

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

1;
__END__

=head1 NAME

cnv_sim2norm.pl - Normalized similarity values of matrix

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_sim2norm.pl -i original.sim -o normalized.sim

=head2 Required Arguments

    --infile        # Path to the input file of sim values
                    # STDIN used if not given
    --outfile       # Path to the output file of sim values
                    # STDOUT used if not given

=head1 DESCRIPTION

This program will take an input matrix of similarity values and convert
them to normalized similarity values. The values will be normalized using:

  norm[Xij] = val[Xij] / min ( val[Xii], val[Xjj] )

This will produce a normalized value ranging from 0 to 1.

The input file must be a three column text file where column 1 is the
row position in the matrix, column 2 is the column position in the matrix
and column 3 is the similarity value. The output file will also return a
three column, tab delimited text file where column 1 is the row position
in the matrix, column 2 is the column position in the matrix and column
3 is the normalized similarity value. It is also possible to inlude
the original similarity value in the output using the --include-val option.
The orignal value will be returned in column 4.

This program should work on any positive similarity value on datasets
that do not use 0 for self similarity. This program is suitable for
normalizing BLAST bitscore values, Smith-Waterman similarity values, and
FASTA Z-score values.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file path is not given, the program
will attempt to read input from standard input.

=item -o,--outfile

Path of the output file. If an output fil path is not given, the program
will write output to standard output.

=back

=head1 OPTIONS

=over 2

=item --include-val

Include the original similarity value in the output file. This will return
four columns of data, where the fourth column is the original value before
being normalized.

=item --usage

Print a short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=back

=head1 EXAMPLES

=head2 Typical Use

The typical use for this program is to convert a file of similarity values
to a file of normalized similarity values.

    cnv_sim2norm.pl -i original_values.sim -o normalized_vals.sim

=head2 Workflow Options

The following option is currently not supported.
Because this program can accept input from STDIN, it is possible to
include this as a step to convert BITSCORES from an all by all blast to
normalized similarity values.

    blastall ... | cnv_blast2sim.pl | cnv_sim2norm.pl

=head1 DIAGNOSTICS

=over

=item Can not open input file

The input file path provided by the -i, --infile switch is not valid.

=item Can not open output file

The output file path provided by the -o, --outfile witch is not valid.

=item Expecting input from STDIN

When an input file path is not specified, the program expects input 
to come from STDIN. This default behaviour allows output from repseek
to be piped directly into the cnv_sim2norm.pl program.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of configuration files or variables defined
in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

This program is not dependent on external programs.

=head2 Required Perl Modules

=over 2

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

This program currnetly not currently support accepting input from
standard input.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 CITATION

A manuscript is in preparation describing this software. Currently you 
should cite the repminer website:

   JC Estill, RS Baucom and JL Bennetzen. 2008. RepMiner.
   http://repminer.sourceforge.net

=head1 HISTORY

STARTED: 01/23/2009

UPDATED: 01/23/2009

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 01/23/2009
# -Main body of program written
# To DO:
#  - THis is currnelty not working with STDIN because it is
#    trying to iterate acrosst the STDIN stream twice.
#    For STDIN to work will need load this to sparse matrix
#    @sim_vals = [row_num][i][j][val]
#    and then get the diagonals for
#    if (i==j) {diag_val = val}
#   Then iterate across sparse matrix and use diag_vals to

#   lookup for getting the min. The current method is the
#   least memory intenstive, and will be left for now.
