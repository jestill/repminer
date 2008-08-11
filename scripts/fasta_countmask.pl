#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_countmask.pl - count fasta file masked characters   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/11/2008                                       |
# UPDATED: 08/11/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Counts the masked characters in a fasta file. Reports    |
#  the number of masked characters, the total number        |
#  of characters and the final ratio. Assummes lower        |
#  case masking.                                            |
#                                                           |
# USAGE:                                                    |
#  fasta_countmask.pl -infile.fasta -o outreport.txt        |
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

# INFILE
if ($infile) {
    open (FASTIN, "<$infile");
}
else {
    open (FASTIN, "<&STDIN");
}

# OUTFILE
if ($outfile) {
    open (SUMOUT, ">$outfile");
} 
else {
    open (SUMOUT, ">&STDOUT");
}


my $num_mask_char = 0;       # Global number of masked characters
my $num_char = 0;      # Global number of characters
my $line_num;
while (<FASTIN>) {

    $line_num++;
    
    # Skip FASTA headers
    next if m/^>/; 
    chomp;

    # Use the translate regular expression without 
    # replacement to count the lowercase characters
    my $num_line_mask_char = $_ =~ tr/[a-z]//;
    my $num_line_char = length $_;

    
    # Add to the running global count
    $num_mask_char = $num_mask_char + $num_line_mask_char;
    $num_char = $num_char + $num_line_char;

    # Show what is being parsed for debug purposes
    print STDERR $_."\n" if $verbose;
    print STDERR "LWC CHAR: $num_line_mask_char\n" if $verbose;
    print STDERR "TOT CHAR: $num_line_char\n" if $verbose;
    

}

my $percent_mask_char = $num_mask_char / $num_char;


print SUMOUT "Masked:\t$num_mask_char\n";
print SUMOUT "Total:\t$num_char\n";
print SUMOUT "Ratio:\t$percent_mask_char\n";


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

fasta_countmask.pl -  count fasta file masked characters

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    fasta_countmask.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

Counts the number of masked characters in a fasta file.

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

STARTED: 08/11/2008

UPDATED: 08/11/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
