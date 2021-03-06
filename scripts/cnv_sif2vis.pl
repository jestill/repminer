#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_sif2vis.pl - Convert sif files to VisAnt Text Files   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/29/2008                                       |
# UPDATED: 05/04/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert SIF format files for Cytoscape to the VisANT     |
#  simple text file format.                                 |
#                                                           |
# USAGE:                                                    |
#  cnv_sif2vis.pl -i infile.sif -o outfile.txt              |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

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


#-----------------------------+
# OPEN THE I/O FILES          |
#-----------------------------+
# STDIN and STDOUT as defaults
if ($infile) {
    open (SIFIN, "<$infile" ) ||
	die "Can not open the sif input file $infile\n";
}
else {
    print STDERR "Expecting input from STDIN\n";
    open (SIFIN, "<&STDIN") ||
	die "Can not accept input from STDIN.";
}

if ($outfile) {
    open (VISOUT, ">$outfile") ||
	die "Can not open the Visant output file $outfile\n";
}
else {
    open (VISOUT, ">&STDOUT") ||
	die "Can not print output to STDOUT\n";
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
while (<SIFIN>) {

    chomp;
    my @sif_in = split;
 
    # Inititally going to assume node 1 to node 2 direction for the file
    print VISOUT $sif_in[0]."\t".  # Node 1
	$sif_in[2]."\t".           # Node 2
	"1\t".                     # Direction 
	$sif_in[1]."\t".           # Method ID
	"1".                       # Edge weight
	"\n";
	
}

close (SIFIN);
close (VISOUT);

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file. If an input file is not specified, the program
will expect intput from STDIN.

=item -o,--outfile

Path of the output file. If an output file path is not specified, the program
will print output to STDOUT.

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

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

The following messages are associated with this program.

=over 2

=item Expecting input from STDIN

If you did not specify and input file with the -i or --infile option, the
program will expect that input will be coming from STDIN. The use of
STDIN will allow the program to be used as part of an analysis pipeline.

=back

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

=head 2 Required Software

Although this conversion program does not rely on external software, the output
is desiged to convert Cytoscape format files to the VISant format files.

=over 2

=item Cytoscape

Cytoscape is a network visualization program that is available at
http://www.cytoscape.org/.

=item VisANT

The VisANT program is a network drawing program available at
http://visant.bu.edu/. VisANT can be downloaded to run locally using a local
installation of java, or can be run as a java web start program from the
VisANT web page.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 05/29/2008

UPDATED: 05/04/2009

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 05/04/2009
# - Added input from STDIN and output to STDOUT
