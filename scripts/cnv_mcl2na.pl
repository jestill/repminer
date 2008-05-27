#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_mcl2na.pl - Convert MCL output to node attributes     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill@gmail.com                            |
# STARTED: 08/24/2007                                       |
# UPDATED: 02/06/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
# Convert MCL format files to NA files suitable for use in  |
# visualizing results in Cytoscape as node attributes       |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+ 
use strict;
use Getopt::Long;



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

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
                    "i|infile=s" => \$infile,
                    "o|outfile=s"=> \$outfile,
		    # ADDITIONAL OPTIONS
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
if ($show_usage) {
    print_help("");
}

if ($show_help || (!$ok) ) {
    #rint_help("full");
    system("perldoc $0");
    exit($ok ? 0 : 2);
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
# MAIN BODY OF PROGRAM                                      |
#-----------------------------------------------------------+

open (IN, "<".$infile) ||
    die "Could not open input file $infile\n";

# Open the output *.NA file for output and
# write the header.
open (OUT, ">".$outfile) ||
    die "Could not open outfile $outfile";


while (<IN>) {
    $clust_num++;
    chomp;
    print "Processing cluster $clust_num\n";
    
    my @full_seq_ids = split;

    while my $full_seq_id (@full_seq_ids) {

    }

    #print OUT $split_in[0]."\n"; 


	
}


#

=head1 NAME

cnv_mcl2na.pl - Convert MCL output to node attributes

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_mcl2na.pl -i InFile -o OutFile
    
    -i,--infile      # Path to the input file
    -o,--outfie      # Path to the output file

=head1 DESCRIPTION

Converts output from the MCL program to a format that can be 
visualized in Cytoscape.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2


=item -i,--infile

Path of the input file. This is the output file from mcl.

=item -o,--outfile

Path of the output file. This is the node atrribute *.NA file
that can be opened in Cytoscape.

=back

=head1 Additional Options

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

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

This program does not use configuration files or make
 use of variables set in the user
environment 

=head1 DEPENDENCIES

=head2 Software

This program use output from the program mcl.
http://micans.org/mcl/

=head1 BUGS AND LIMITATIONS

Currently no bugs are known for this program.

This script has been verified to work with mcl 1.006, 06-058

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/06/2008
# - Updating code a bit to make documentation easier
