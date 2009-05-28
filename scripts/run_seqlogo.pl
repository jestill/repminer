#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# test_seqlogo.pl                                           |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/08/2009                                       |
# UPDATED: 05/08/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Test running the seq logo program for groups of sequence |
#  data that have been extraced from a sequence database.   |
#                                                           |
# USAGE:                                                    |
#  test_seqlogo.pl -i infile -o outfile.png                 |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
#

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
# MAIN PROGRAM BODY           |
#-----------------------------+
# Goal is to run the seq logo program
# Assume seqlogo in path
my $seqlogo_bin = "/home/jestill/apps/weblogo/seqlogo";

# For full seqs do below
#my $seqlogo_cmd = "$seqlogo_bin".
#    " -f $infile".
#    " -o $outfile".
#    "  -F PNG -c -b -a -w 18 -h 5 -n -s \"-26\" -Y -k 1 -t Ji -S";
##./seqlogo -f ji_test.txt -o test_ji -F PNG -c -b -a -w 18 -h 5 -n -s "-26" -Y -k 1 -t Ji -S
#print "CMD: $seqlogo_cmd\n" if $verbose;

# for significant do below
my $seqlogo_cmd = "$seqlogo_bin".
    " -f $infile".
    " -o $outfile".
    "  -F PNG -c -b -a -w 18 -h 5 -n -s \"-26\" -Y -k 1 -t Ji";
system($seqlogo_cmd);


# ALSO NOTE THAT SEQ LOGO WILL ACCEPT STDIN 
# BY SPECIFYING - AS INFILE WITH -f OPTION
#cat ji_test.txt | ./seqlogo -f - -o test_pipe.png


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

STARTED: 05/08/2009

UPDATED: 05/08/2009

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#



# SEQLOGO NOTES
#usage: seqlogo -f <input filename> [OPTIONs with values]
#Creates a logo for the given input filename.
# 
#Available options:
#  -B <bar bits>              Number of bits in bar (real # > 0)
#  -T <tic bits>              Number of bits between tic marks
#  -C <chars per line>        Number of characters per line of logo
#  -d <box shrink factor>     Shrink factor of characters if option c is toggled
#  -E <error bar fraction>    Fraction of error bar to show (real # > 0)
#  -f <input filename>        Input filename
#  -F <format>                Format of output (EPS, GIF, PDF, PNG), - for STDOUT  -h <logo height>           Height of output logo (real # > 0)
#  -k <kind of data>          0 for amino acid, 1 for nucleic acid
#  -l <sequence lower bound>  Lower bound of sequence (integer)
#  -m <sequence upper bound>  Upper bound of sequence (integer)
#  -o <output file>           Name of output file
#  -s <sequence start>        Sequence start number, defaults to 1 (int)
#  -t <titletext>             Text of title, enclosed in "" if more than one word  -w <logo width>            Width of output logo
#  -x <x-axis label>          Label for x-axis
#  -y <y-axis label>          Label for y-axis
# 
#Available toggles (no values associated) bOenc
#  -a       Toggle antialiasing
#  -b       Toggle bar ends
#  -c       Toggle color
#  -e       Toggle error bar
#  -M       Toggle small sample correction
#  -O       Toggle outlining of characters
#  -p       Toggle fineprint
#  -n       Toggle numbering of x-axis
# -S       Toggle stretching of logos to entire length
# -X       Toggle boxing of characters
# -Y       Toggle y-axis
 
