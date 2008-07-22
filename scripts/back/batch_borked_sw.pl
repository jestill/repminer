#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_prss34.pl - Run all by all prss34 on fasta file     |
#                   May need to switch to dir of seqs       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/13/2008                                       |
# UPDATED: 02/16/2008                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Given a multi-fasta file run the prss34 program to       |
#  generate a distance matrix. This will currently create   |
#  output in the format needed for the apclust binary.      |
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
use Bio::SeqIO;
use IPC::Open2;                # For interaction with prss34

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outdir;

# Index vals
my $i=0;
my $j=0;

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
                    "o|outdir=s"  => \$outdir,
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
# CHECK FOR REQUIRED COMMAND  |
# LINE OPTIONS                |
#-----------------------------+


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}


#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
        die "ERROR: Could not create the output directory:\n$outdir";
}

#-----------------------------+
# OUTPUT FILE HANDLES         |
#-----------------------------+
# VECOUT is the vector of nodes
# SIMOUT is the vector of similarity scores
open (VECOUT, ">".$outdir."vec.txt") ||
    die "Can not open vector output file\n";

open (SIMOUT, ">".$outdir."similarity.txt") ||
    die "Can not open similarity output file\n";

open (IDOUT, ">".$outdir."vecname.txt") ||
    die "Can not open the vector output name file\n";

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
# Lazily just duplicate the infile in memory

my $seqs_1 = Bio::SeqIO->new( -file => "<$infile",
			      -format => 'fasta') ||
    die "Can not open input sequence file one:\n$infile\n";

# Need to open the second object of the sequence
# below
#my $seqs_2 = Bio::SeqIO->new( -file => "<$infile",
#			      -format => 'fasta') ||
#    die "Can not open input sequence file two:\n$infile\n";
    

#-----------------------------+
#                             |
#-----------------------------+
while (my $seq1 = $seqs_1->next_seq) {

    $i++;

    # Print output to vector dir
    print VECOUT "$i\n";
    print IDOUT "$i\t".$seq1->primary_id()."\n";
    
    my $seqs_2 = Bio::SeqIO->new( -file => "<$infile",
				  -format => 'fasta') 
	|| die "Can not open input sequence file two:\n$infile\n";

    while (my $seq2 = $seqs_2->next_seq) {
	
	$j++;
	
	my $qry_str = ">".$seq1->primary_id."\n".$seq1->seq()."\n".
	    ">".$seq2->primary_id."\n".$seq2->seq()."\n";
	
	#print STDERR $qry_str if $verbose;

	#-----------------------------+
	# Try to run prss34 with pipes|
	#-----------------------------+
#	my $cmd = $qry_str." | prss34\n"; 
#	# Backtikcs can capture the output of a program
#	$date = `/bin/date`;

	#-----------------------------+
	# TRYING WITH SYSTEM CMD
	#-----------------------------+
#	my $cmd = " prss34";
#	system ($cmd);

	#-----------------------------+
	# TRYING WITH OPEN 2
	#-----------------------------+

	print "QUERY:\n$qry_str\n";

	my $pid = open2(*Reader, *Writer, "prss34 -q - - -m 9" );
	print Writer "$qry_str";
	my $got = <Reader>;

	print "GOT:\n\n$got\n";
	
	exit;

# CMD AS FILE HANDLE
#	open (PRSSRUN, "| prss34 ") ||
#	open (PRSSRUN, "| $cmd") ||
#	    die "Can not run the prss34 program\n";
#	print PRSSRUN $qry_str;
#
#	close (PRSSRUN);
#
    } # End of interate across second seq object

} # End of iterate across first seq object

#-----------------------------+
# CLOSE OPEN FILES            |
#-----------------------------+
close (VECOUT);
close (SIMOUT);
close (IDOUT);

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

STARTED: 02/13/2008

UPDATED: 02/16/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
