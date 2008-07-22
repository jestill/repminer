#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_sif2dbtxt.pl - Convert sif files to Databse text file |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/10/2008                                       |
# UPDATED: 06/10/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert SIF format files for Cytoscape to a text file    |
#  that contains a unique ID for each edge. This makes      |
#  the format suitable for upload to a database using the   |
#  edge ID as the unique identifier. Using this unique      |
#  edge ID can allow for the removal of redundancies by     |
#  using the edge ID as the unique key in the database      |
#  table                                                    |
#                                                           |
# USAGE:                                                    |
#  cnv_sif2dbtxt.pl -i infile.sif -o outfile.txt            |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Path to the infile
my $outfile;                   # Path to the outfile
my $outformat = "unix";        # The output format (DOS|UNIX|MAC)

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
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    # ADDITIONAL OPTIONS
		    "out-format=s" => \$outformat,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

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
open (SIFIN, "<$infile" ) ||
    die "Can not open the sif input file $infile\n";

open (TXTOUT, ">$outfile") ||
    die "Can not open the Visant output file $outfile\n";

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
while (<SIFIN>) {

    chomp;

    my @sif_in = split;
 
    print TXTOUT $sif_in[0].$sif_in[1].$sif_in[2]."\t". # Unique EDGE ID
	$sif_in[0]."\t".                                # Node 1
	$sif_in[2]."\t".                                # Node 2
	$sif_in[1];                                     # Attribute 

    if ($outformat =~ "dos") {
	print TXTOUT "\r\n";                            # Did not parse 
    }
    elsif ($outformat =~ "unix") {
	print TXTOUT "\n";
    }
    elsif ($outformat =~ "mac") {
	print TXTOUT "\r";
    }
    # default 
    else {
	print TXTOUT "\n";
    }

	
}

close (SIFIN);
close (TXTOUT);

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

Path of the input file.

=item -o,--outfile

Path of the output file.

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

STARTED: 06/10/2008

UPDATED: 06/10/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/10/2008
# - Program modified from cnv_sif2vis.pl
