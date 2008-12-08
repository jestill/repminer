#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_mcl2na.pl - Convert MCL output to node attributes     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/24/2007                                       |
# UPDATED: 12/07/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
# Convert MCL format files to NA files suitable for use in  |
# visualizing MCL cluster results in Cytoscape as node      |
# attributes.                                               |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# TO DO:
# Load the node attributes to a BiSQL based database.
# This will also need to load the mcl meta information to the db

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

# VARIABLES FOR COMMAND LINE
my $infile;
my $outfile;

# INITIALIZE COUNTERS
my $clust_num = 0;

# BOOLEANS
my $show_help = 0;             # Display help
my $quiet = 0;                 # Run the program in quiet mode
                               # will not prompt for command line options
my $show_node_id = 0;          # Include the database node_id in the output
my $show_man = 0;              # Show the man page via perldoc
my $show_usage = 0;            # Show the basic usage for the program
my $show_version = 0;          # Show the program version
my $verbose = 0;    

my $param_set = "DEF";         # assume a default parameter set for mcl

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
                    "i|infile=s" => \$infile,
                    "o|outfile=s"=> \$outfile,
		    # ADDITIONAL OPTIONS
		    "p|param"    => \$param_set,
		    # BOOLEANS
		    "q|quiet"    => \$quiet,
                    "verbose"    => \$verbose,
		    "version"    => \$show_version,
		    "man"        => \$show_man,
		    "usage"      => \$show_usage,
		    "h|help"     => \$show_help,);

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

#-----------------------------------------------------------+
# MAIN BODY OF PROGRAM                                      |
#-----------------------------------------------------------+

#-----------------------------+
# FILE I/O                    |
#-----------------------------+
if ( $infile ) {
    open (IN, "<".$infile) ||
	die "Could not open input file $infile\n";
}
else {
    open (IN, "<&STDIN");
}

#-----------------------------+
# OPEN THE OUTPUT *.NA FILE   |
# AND WRITE THE HEADER        |
#-----------------------------+
if ( $outfile ) {
    open (OUT, ">".$outfile) ||
	die "Could not open outfile $outfile\n"; 
}
else {
    open (OUT, ">&STDOUT")
}
print OUT "MCL_CLUST_".$param_set."\n";


while (<IN>) {

    $clust_num++;
    chomp;
    print "Processing cluster: $clust_num\n" if $verbose;

    # Currently just splitting on the whitespace
    # splitting on tab would be better

    my @seq_ids = split;

    # For each seq record print its cluster id
    # and the name
    for my $seq_id (@seq_ids) {
        # TO DO: split again by pipe to get the number to
        #        use to map node attributes in cytoscape
        my ($id_num,$rest) = split(/\|/,$seq_id);

        print OUT "$id_num=MCL_Clust_$clust_num\n";

        print STDERR "\tCLUST:$clust_num:\t$seq_id" if $verbose;
    }

}

close (OUT);
close (IN);

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

cnv_mcl2na.pl - Convert MCL output to node attributes

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2  Usage

    cnv_mcl2na.pl -i InFile -o OutFile

=head2 Required Arguments

    -i      # Path to the input file
    -o      # Path to the output file
    -p      # Id of the MCL parameter set

=head1 DESCRIPTION

Converts output from the MCL program to a format that can be 
visualized in Cytoscape.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This is the output file from mcl. If no input file
path is specified, cnv_mcl2na.pl will expect input from STDIN.

=item -o,--outfile

Path of the output file. This is the node atrribute *.NA file
that can be opened in Cytoscape. If no output file path is specified, 
cnv_mcl2na.pl will write output to STDOUT.

=item -p,--param

The id of the MCL parameter set that was used for MCL. This will be
appended to the header name of the NA file that is produced. If no value is
specified, cnv_mcl2na.pl will assume that the default parameter set was used 
and will tag the parameter set as DEF.

=back

=head1 OPTIONS

=over 2

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -v, --verbose

Run the program in verbose mode.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 EXAMPLES

To convert an output file from mcl to the Cytoscape node attribute format

 cnv_mcl2na.pl -i mcl_result.mcl -o mcl_node_attribute.NA

The cnv_mcl2na.pl program will assume that the MCL program was run using
default paramters, and the text 'DEF' will be appended to the header 
of the yielding a file like:

 MCL_CLUSTER_DEF
 seq_1   MCL_CLUST_2
 seq_2   MCL_CLUST_9
 seq_3   MCL_CLUST_2

Other parameter set identifiers can be appened to the header using the -p
flag. For example, if you ran MCL with the inflation value of 3 you
may indicate this in your output file name and the p flag as:

 cnv_mcl2na.pl -i mcl_result.mcl -o mcl_node_attribute_i_3.NA -p I_3

This will generate a file like:

 MCL_CLUSTER_I_3
 seq_1   MCL_CLUST_2
 seq_2   MCL_CLUST_9
 seq_3   MCL_CLUST_2

Using the paramter file and different file names will allow you to compare
mcl results from multiple paramter sets within the same Cytoscape session.

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

=over 2

=item MCL

This program use output from the program mcl.
http://micans.org/mcl/. 
"The MCL algorithm is short for the Markov Cluster Algorithm, a fast and 
scalable unsupervised cluster algorithm for graphs based on simulation 
of (stochastic) flow in graphs."

=item Cytoscape

This program generates a node attribute file for the network visualization
program Cytoscape http://cytoscape.org/.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the RepMiner
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

This script has been verified to work with mcl 1.006, 06-058

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/24/2007

UPDATED: 12/07/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/06/2008
# - Updating code a bit to make documentation easier
#
# 12/06/2008
# - Updateing help file
