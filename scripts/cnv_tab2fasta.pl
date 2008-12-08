#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cvn_txt2fasta.pl - Converts tab delim text to fasta       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/01/2008                                       |
# UPDATED: 12/07/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#                                                           |
# Converts a tab delimitted text file into the fasta        |
# format. The last column is assumed to be the sequence     |
# data. Preceeding columns are added to the fasta header    |
# file separated by the pipe delimiter.                     |
#                                                           |
# USAGE:                                                    |
#  cnv_tab2fasta.pl -i infile.txt -o outfile.fasta          |
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
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

# INFILE
if ($infile) {
    open (TXTIN, "<$infile");
}
else {
    print STDERR "Waiting for input from STDIN\n";
    open (TXTIN, "<&STDIN");
}

# OUTFILE
if ($outfile) {
    open (FASTOUT, ">$outfile");
} 
else {
    open (FASTOUT, ">&STDOUT");
}


my $line_num = 0;
while (<TXTIN>) {

    chomp;

    $line_num++;

    print STDERR "Processing $line_num" if $verbose;
    my @in_cols = split;
    my $num_in_cols = @in_cols;
    
    for (my $i = 0; $i<$num_in_cols; $i++ ) {

	# Variable number column parse to FASTA
	if ($i == 0) {
	    # First column - Start FASTA header
	    print FASTOUT ">";
	    print FASTOUT $in_cols[$i];
	}
	elsif ($i == $num_in_cols - 1) {
	    # Last column - Print sequence string and newlines
	    print FASTOUT "\n";
	    print FASTOUT $in_cols[$i];
	    print FASTOUT "\n";
	}
	else {
	    # All other columns - prefix with pipe
	    print FASTOUT "|";
	    print FASTOUT $in_cols[$i];
	}

    }

}


close (FASTOUT);
close (TXTIN);

exit;

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

cnv_tab2fasta.pl - Convert tab delim text file to fasta format

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_tab2fasta.pl -i infile.txt -o outfile.fasta

=head2 Required Arguments

    -i        # Path to the input text file
    -o        # Path to the output fasta file

=head1 DESCRIPTION

Converts a tab delimitted text file into the fasta format. The last
column in the tab delimited text file is assumed to be the sequence data. 
Preceeding columns
are added to the fasta header file separated by the pipe
delimiter. It is possible to have variable number of columns for each
sequence record, as long as the last column is always the sequence string.
The typical use for this program will be to convert a text dump for a 
database or Microsoft Excel to a fasta format text file.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the tab delimited text input file. If a file is not specified by
the --infile option, the program will attempt to read input from STDIN.

=item -o,--outfile

Path of the fasta format output file. If a file is not specified by the
--outfile option, the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --usage

Short program description of usage.

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

=head1 EXAMPLES

=head2 Typical Use

In the typical use case, you will want to convert a tab delimitted text
file produced from a database dump to the fasta file format

 cnv_tab2fasta.pl -i infile.txt -o outfile.fasta

This will convert an intput text file in the following format

 seq_1   maize   chrom_1 ATGCAA...
 Seq_2   maize   chrom_1 GCCATC...
 seq_3   sorghum chrom_1 TCAGCC...

To the following format

 >seq_1|maize|chrom_1
 ATGCAA...
 >seq_2|maize|chrom_1
 GCCATC...
 >seq_3|sorghum|chrom_1
 TCAGCC...

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item Program does not appear to be doing anything

It is possible that the program is waiting for input from STDIN.
Run with the --verbose mode. In these cases, you should see the text
'Waiting for input from STDIN'.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This script does not make use of external configuration files
or any options in the user environment.

=head1 DEPENDENCIES

This program is not dependent on external programs or any perl modules outside
of the Getopt::Long which is included in most perl distributions.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the RepMiner
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head1 SEE ALSO

The cnv_tab2fasta.pl program is part of the repminer package of 
repeat element annotation programs.
See the DAWG-PAWS web page 
( http://repminer.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/repminer/ )
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 11/01/2008

UPDATED: 12/07/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/01/2008
# - Don't actually know when this was writtien, probably
#   around november 2008
#
# 12/07/2008
# - Updating POD documentation
