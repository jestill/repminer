#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_sim2fil.pl - Convert sim matrix to cytoscape format   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/23/2009                                       |
# UPDATED: 01/23/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts similarity matrix data to the sif format for    |
#  the Cytoscape network visualization program.             |
#                                                           |
# USAGE:                                                    |
#  cnv_sim2sif.pl -i infile.sim -o outfile.sif              |
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

# Graph options
my $graph_dir = "0";
my $rel = "sim";                  # Relationship type
my $ea_file;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # GRAPH OPTIONS
		    "r|relate=s"    => \$rel,
		    "d|direction=s" => \$graph_dir, # Graph direction
		    "e|ea=s"        => \$ea_file,   # Edge attribute 
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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
    open (SIFOUT, ">$outfile") ||
	die "Can not open output file:\n$outfile\n";
}
else {
    open (SIFOUT, ">&STDOUT") ||
	die "Can not access standard output\n";
}

# EDGE ATTRIBUTE FILE
if ($ea_file) {
    open (EAOUT, ">$ea_file");
    print EAOUT uc($rel)."\n";
}

#-----------------------------+
# ITERATE THROUGH AND CONVERT |
#-----------------------------+
# What is reported is based on the graph_dir options
while (<SIMIN>) {
    chomp;
    my ($x,$y,$val) = split;

    # Cytoscape only seems to play will with integers
    $val = int($val);

    #-----------------------------+
    # UNDIRECTED x < y            |
    #-----------------------------+
    if ($graph_dir == '0') {

	if ($x < $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t$rel\t".$y."\n";
	    print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);

	} # End of if XCrd > YCrd
	
    }

    #-----------------------------+
    # UNDIRECTED x != y           |
    #-----------------------------+ 
    # Using ($rel) for both directed makes digraph into unigraph
    # This is sloppy and depends on cytoscape to interpret
    elsif ($graph_dir == '1') {
	
	if ($x != $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t$rel\t".$y."\n";
	    print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);
	    
	} # End of if XCrd > YCrd
	
    } 
    
    
    #-----------------------------+
    # UNDIRECTED ALL x and y      |
    #-----------------------------+
    elsif ($graph_dir == '2') {
	
	print STDERR $x."-->".$y."\n" if $verbose;

	# CYTOSCAPE FILES
	print SIFOUT $x."\t$rel\t".$y."\n";
	print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);
	
    } 
    
    
    #-----------------------------+
    # DIRECTED x < y              |
    #-----------------------------+
    elsif ($graph_dir == '3') {
	if ($x < $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    print STDERR "Adding Edge $x-->$y\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t$rel\t".$y."\n";
	    print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);
	    
	} # End of if XCrd > YCrd
	
    }
    
    #-----------------------------+
    # DIRECTED x != y             |
    #-----------------------------+
    elsif ($graph_dir == '4') {
	# The use of bot ensures that the x < y is treated as 
	# different information then x > y by cytoscape. 
	# bot refers to the bottom part of the all-by-all matrix
	
	
	if ($x < $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t$rel\t".$y."\n";
	    print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);

	} # End of if XCrd > YCrd
	
	if ($x > $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t".$rel."_rev\t".$y."\n";
	    print EAOUT "$x (".$rel."_rev) $y = $val\n" if ($ea_file);
	} # End of if XCrd > YCrd
	
    }
    
    
    #-----------------------------+
    # DIRECTED ALL x and y        |
    #-----------------------------+
    elsif ($graph_dir == '5') {
	
	# All x and y
	# The use of bot ensures that the x < y is treated as 
	# different information then x > y.
	
	if ($x <= $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES
	    print SIFOUT $x."\t$rel\t".$y."\n";
	    print EAOUT "$x ($rel) $y = $val\n" if ($ea_file);
	    
	} # End of if XCrd > YCrd
	
	if ($x > $y) {
	    
	    print STDERR $x."-->".$y."\n" if $verbose;
	    
	    # CYTOSCAPE FILES	    
	    print SIFOUT $x."\t".$rel."_rev\t".$y."\n";
	    print EAOUT "$x (".$rel."_rev) $y = $val\n" if ($ea_file);

	} # End of if XCrd > YCrd
	
    } # End of if statements for GraphDir
    
} # End while simin

close SIFOUT;
close SIMIN;
close EAOUT if ($ea_file);

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

1;
__END__

=head1 NAME

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_sim2sif.pl -i infile.sim -o outfile.sif

=head2 Required Arguments

    --infile        # Path to the input file
                    # Attempts to use STDIN if not specified
    --outfie        # Path to the output file
                    # Prints to STDOUT if not specified

=head1 DESCRIPTION

This program will convert a similarity matrix file to the simple interaction
format (sif) used by the Cytoscape network visualization program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.  If an input file path is not given, the program
will attempt to read input from standard input.

=item -o,--outfile

Path of the output file. If an output file path is not given, the program
will write output to standard output.

=back

=head1 OPTIONS

=over 2

=item -r,--relate

The type of relationship represented by the similarity matrix. By default
this will be set to 'sim'. Relationships below the diagonal will have the
value _rev appened to indicate reverse. This assures that all edges will
be drawn. This will be the header used to name the Edge Attributes file that
is created.

=item -e,--ea

The path to an edge attribute file. If this argument is not given, an
edge attribute file will not be generated.
This option will generate an edge attribute file in addition to the sif file.
This ea file will give the values in the similarity matrix.
More information on the Cytsocape edge attribute file format is 
available at:
http://cytoscape.org/cgi-bin/moin.cgi/Cytoscape_User_Manual/Attributes

=item -d,--direction

Graph edge direction.

-d 0 = undirected (i < j)

-d 1 = undirected (i != j)

-d 2 = undirected (all i,j)
       Reciprocal hits reduced to a single edge.

-d 3 = directed (i < j)

-d 4 = directed (i != j)

-d 5 = directed (all i,j)
       Reciprocal hits drawn as two edges.

=item --usage

Short overview of how to use program from command line.

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

The typical use for this program is to convert a similarity file to a 
sif file for Cytoscape.

    cnv_sim2sif.pl -i infile.sim -o outfile.sif

=head2 Generate Edge Attribute Files

For the values of the edges to be visible in cytoscape, you will need to
be generate an edge attributes file.

    cnv_sim2sif.pl -i infile.som -o outfile.sif --ea sim_value.ea 

This will generate an edge attribute file like the following

    EA
    10 (sim) 10 = 4910
    26 (sim) 26 = 5079
    42 (sim) 42 = 2803
    42 (sim) 377 = 71.9
    58 (sim) 58 = 1.352e+04
    ...

You can control the header used in the edge attribute file with the
--relate variable. This sets the string used to describe the relationship in 
the sif file and the ea file. This uppercase translation of --relate will
be used to generate the ea header. For example 

    cnv_sim2sif.pl -i infile.som -o outfile.sif --ea sim_value.ea 
                   --relate blast

Will generate an edge attribute file like the following:

    BLAST
    10 (sim) 10 = 4910
    26 (sim) 26 = 5079
    42 (sim) 42 = 2803
    42 (sim) 377 = 71.9
    58 (sim) 58 = 1.352e+04
    ...

More information on the cytoscape edge attribute file format is available
on the Cytoscape wiki page at:
http://cytoscape.org/cgi-bin/moin.cgi/Cytoscape_User_Manual/Attributes

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

=head2 Required Software

=head2 Required Perl Modules

=over 2

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

=over 2

=item * Node Identifiers are Limited to Integers

The ability to specify edge direction will only work if the first two
values in the sim file (these are the node identifiers) are integers.

=back

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
# - Bulk of program written. Can convert sim to sif with
#   multiple edge direction options.
