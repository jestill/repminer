#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ap2na.pl - Convert Affinity Prop. to Cytoscape NA     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/04/2008                                       |
# UPDATED: 06/04/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert the affinity propagation output to Cytoscape     |
#  node attributies. Will create the AP Cluster membership  |
#  and can also identify the exemplars in the graph.        |
#                                                           |
# USAGE:                                                    |
#  cnv_ap2na.pl -i apfile.txt -o node_attributes.NA         |
#               -h node_attribute_header                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO: Edge attributes as directed graphs


package REPMINER;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
# String variables set at the command line
my $infile;                    # Path to the input file to parse
my $outfile;                   # Path to the cluster id node attribute file
my $header;                    # header to use is the outfile
my $exemplar_out;              # Path to the outfile of exemplars
my $exemplar_header;           # exemplar header
my $edge_out;                  # Path to the outfile for the edge attributes
my $edge_header;               # Header for the edge attribute file

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
                    "o|outfile=s"   => \$outfile,
		    # ADDITIONAL OPTIONS
		    "x|exemplar=s"  => \$exemplar_out,
		    "ex-header=s"   => \$exemplar_header,
		    "e|edge=s"      => \$edge_out,
		    "edge-header=s" => \$edge_header,
		    "header=s"      => \$header,
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
if ($show_usage) {
    print_help("");
}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\n$0:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

#-----------------------------------------------------------+
# FILE HANDLES                                              |
#-----------------------------------------------------------+

#-----------------------------+
# INPUT FILE                  |
#-----------------------------+
open (INFILE, "<$infile") ||
    die "Can not open the affinity propagation input file\n";

#-----------------------------+
# AFFINITY PROP CLUSTERS      |
#-----------------------------+
open (OUTFILE, ">$outfile") ||
    die "Can not open the output file $outfile\n";
if ($header) {
    print OUTFILE "$header\n";
}
else {
    print OUTFILE "AP_CLUSTERS".$infile."\n";
}

#-----------------------------+
# EXEMPLARS                   |
#-----------------------------+
if ($exemplar_out) {
    open (EXOUT, ">$exemplar_out") ||
	die "Can not open the exemplar out file $exemplar_out\n";

    if ($exemplar_header) {
	print EXOUT "$exemplar_header\n";
    }
    else {
	print EXOUT "AP_EXEMPLARS_".$infile."\n";
    }
}

#-----------------------------+
# EDGE ATTRIBUTES             |
#-----------------------------+
# This will attempt to draw directional arrows going from
# a node to its exemplar node. This will may run into
# problems when attemping to add an edge attribute for
# an edge that does not exist
if ($edge_out) {
    open (EDGEOUT, ">$edge_out") ||
	die "Can not open the edge attribute file\n";

    if ($edge_header) {
	print EDGEOUT "$edge_header\n";
    }
    else {
	print EDGEOUT "AP_EDGE_ATTRIBUTE_".$infile."\n";
    }
    
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $i = 0;

while (<INFILE>) {
    
    $i++;
    chomp;

    # Convert input to interger
    # Sometimes the output from the APClust program may be in
    # scientific notation depending on how it was exported

    # The exemplar node id value as an integer
    my $int_val = int($_);
    
    print STDERR "Processing line $i --> $_\n" if $verbose;
    print STDERR "\tINTEGER $int_val" if $verbose;

    #-----------------------------+
    # WRITE AP CLUSTER ID         |
    #-----------------------------+
    unless ($int_val == 0) {
	print OUTFILE $i."=AP_CLUST_".$int_val."\n";
    }
    # What to do with the zero nodes
    # for now will print to the class zero ... don't know
    # what to make of these
    else {
	print OUTFILE $i."=AP_CLUST_".$int_val."\n";
    }
    
    #-----------------------------+
    # WRITE EXEMPLARS TO OUTFILE  |
    #-----------------------------+
    # Exemplars will find themselves as exemplars
    # create a simple binary outfile. 1 is an exemplar
	# otherwise is zero
    if ($exemplar_out) {

	if ($i == $int_val) {
	    print EXOUT $i."=1\n";
	}
	else {
	    print EXOUT $i."=0\n";
	}
    }


    #-----------------------------+
    # WRITE EDGE ATTRIBUTES       |
    #-----------------------------+
    # Two sets of edges are being left out here
    # 1. Hits where the BLAST hit is in the opposite direction
    # 2. Edges where the node is not directly connected
    # To deal with hits that are in the opposite direction, I will
    # include a backwards class.

    if ($edge_out) {

	# The variable joining the edges needs to be the same as the 
	# data type used in the *.sif file. For example, bl below
	# refers to a blast hit
	# Coding like the following to
	# allow for leaving out self hits 
	if ($i == $int_val) {
	    # This edge is its own exemplar
	    print EDGEOUT $i." (bl) ".$int_val." = self\n";
	}
	else {
	    # This edge points to another node as its exemplar
	    # Doing this the following was assures that both are in
	    # the *.sif file
	    print EDGEOUT $i." (bl) ".$int_val." = external\n";
	    print EDGEOUT $int_val." (bl) ".$i." = external_rev\n";
	}

    }
    

} # End of While INFILE


#-----------------------------+
# CLOSE FILE HANDLES          |
#-----------------------------+
close INFILE;
close OUTFILE;
close EXOUT if ($exemplar_out);
close EDGEOUT if ($edge_out);

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}


=head1 NAME

cnv_ap2na.pl. - Convert Affinity Prop. to Cytoscape NA

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_ap2ap.pl -i apinfile.txt -o node_attribute.txt

    --infile        # Path to the affinity propagation input file
    --outfie        # Path to the node attribute output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

=over 2

=item -h, --header

The header to use for the Node Attribute file. If one is not passed at the
command line, the name of the input file is used.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=back

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

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 06/04/2008

UPDATED: 06/04/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/04/2008
# - Started basic program to convert the basic output 
#   from the affinity propagation program to the format
#   needed to display these as Node Attributes in 
#   Cytoscape. This can create the following output files
#    -AP Cluster ID
#        NA file that identifies the affinity propagation
#        cluster id
#    -AP Exemplars
#        NA files that highlights the exemplar nodes
#    -AP Edges
#        EA file that identifies the edges and directions for
#        that allows for the drawing of a node to an exemplar 
#        node
