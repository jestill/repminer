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
# UPDATED: 05/29/2009                                       |
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

package REPMINER;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # Play fair
use Getopt::Long;              # Get options from the command line
use DBI;                       # For connecting to MySQL database

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $in_option ="-";               # This specifies input from STDIN
my $out_option = "-";             # This specified output to STDOUT

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# VARIABLES THAT CAN BE CHANGED WITH
# COMMAND LINE OPTIONS
# Assume seqlogo in path
my $seqlogo_bin = $ENV{SEQLOGO_BIN} || "seqlogo";
my $seqlogo_title = "SeqLogo Result";
my $seqlogo_start = "-26";            # 26 derived from LTR_struc output
my $seqlogo_width = "18";
my $seqlogo_height = "5";
my $do_significant = 0;
my $show_all = 0;
my $x_axis_label = "BP From Insertion Site";

my @logo_files_in;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "i|infile|indir=s"   => \$in_option,
                    "o|outfile|outdir=s" => \$out_option,
		    # SEQLOGO OPTIONS
		    "title=s"            => \$seqlogo_title,
		    "seqlogo-bin=s"      => \$seqlogo_bin,
		    "start=s"            => \$seqlogo_start,
		    "height=s"           => \$seqlogo_height,
		    "width=s"            => \$seqlogo_width,
		    "show-all"           => \$show_all,
		    "x-label"            => \$x_axis_label,
		    # ADDITIONAL OPTIONS
		    "q|quiet"            => \$quiet,
		    "verbose"            => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"              => \$show_usage,
		    "version"            => \$show_version,
		    "man"                => \$show_man,
		    "h|help"             => \$show_help,);

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
# LOAD INPUT FILES TO ARRAY   |
#-----------------------------+
my $indir;
# INPUT IS A FILE
if (-f $in_option) {
    push (@logo_files_in, $in_option);
}
# INPUT IS A FOLDER
elsif (-d $in_option) {
    opendir( DIR, $in_option ) || 
	die "Can't open directory:\n$in_option"; 
    @logo_files_in = grep /\.txt$|\.text$/, readdir DIR ;
    closedir( DIR );
    $indir = $in_option;

    unless ($indir =~ /\/$/ ) {
	$indir = $indir."/";
    }

}
else {
    die "The intput does not appear to bea file or directory";
}

#-----------------------------+
# DEAL WITH OUTPUT OPTIONS    |
#-----------------------------+
my $outdir = $out_option;
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

unless (-e $outdir) {
    print "creating dir: $outdir\n";
    mkdir $outdir ||
	die "Could not creat the output dir:\n$outdir\n";
}

#-----------------------------+
# ITERATE THROUGH EACH FILE   |
# IN THE ARRAY                |
#-----------------------------+
for my $infile (@logo_files_in) {

    print STDERR "Processing $infile\n" if $verbose;
    
    #-----------------------------+
    # Get the root name of the    |
    # file to draw                |
    #-----------------------------+
    my $name_root;
    if ($infile =~ m/(.*)\.text$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($infile =~ m/(.*)\.txt$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    else {
	$name_root = $infile;
    }

    my $outfile = $outdir.$name_root;
    my $infile_path = $indir.$infile;

    my $seqlogo_cmd = "$seqlogo_bin".
	" -f $infile_path".
	" -o $outfile".
	" -t \"$seqlogo_title\"". # quotes allow for spaces in title name
	" -s \"$seqlogo_start\"".
	" -w $seqlogo_width".
	" -h $seqlogo_height".
	" -F PNG -c -b -a -n".
	" -Y -k 1";

   # k is kind of data 1 for DNA
   # a is to draw antialias
   # n is to toggle numbering of x-axis
    
    if ($x_axis_label) {
	$seqlogo_cmd = $seqlogo_cmd." -x "." \"$x_axis_label\"";
    }
    
    
    # This will graphically show all bases, not just significant bases
    if ($show_all) {
	$seqlogo_cmd = $seqlogo_cmd." -S";
    }
    
    system($seqlogo_cmd);
    
    
}




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

UPDATED: 05/29/2009

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
 
