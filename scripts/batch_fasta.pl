#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_prss34.pl - Run all by all prss34 on dir of seqs    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/13/2008                                       |
# UPDATED: 02/18/2008                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Given a directory of fasta filesrun the prss34 program to|
#  to generate a distance matrix. This will currently create|
#  output in the format needed for the apclust binary.      |
#  In the future this could be used to store data in the    |
#  BioSQL graph/tree extension.                             |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;
my $dbdir;
my $seq1_path;
my $seq1_id;
my $seq2_path;
my $seq2_id;
my $program = "fasta34";

# Array to hold all Z scores for query sequence
my @z_scores;

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

# Vars needing script level scope
my $seq_id_qry;                # Seq id of the query seq
my $seq_id_subj;               # Seq id of the subject seq

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # ADDITIONAL OPTIONS
		    "d|dbdir=s"    => \$dbdir,
		    "p|program=s"  => \$program, 
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
# CHECK FOR REQUIRED COMMAND  |
# LINE OPTIONS                |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
        " command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
        " command line\n" if (!$outdir);
    exit 1;
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
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
# GET FASTA FILES FROM INDIR  |
#-----------------------------+
opendir( DIR, $indir ) ||
    die "Can't open directory:\n$indir";
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# ERROR IF NOT FASTA FILES    |
#-----------------------------+
if ($count_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
        "$indir\n".
        "Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# OUTPUT FILE HANDLES         |
#-----------------------------+
# The following files are for clustering
# VECOUT is the vector of nodes
# SIMOUT is the matrix of similarity scores
# IDOUT
# PREOUT is the preferences file. This may need to
# be created using a database approach

open (SIMOUT, ">".$outdir."similarity.txt") ||
    die "Can not open similarity output file\n";

open (IDOUT, ">".$outdir."vecname.txt") ||
    die "Can not open the vector output name file\n";

open (PREFOUT, ">".$outdir."preferneces.txt") ||
    die "Can not open the preferences file\n";

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

for my $ind_file_1 (@fasta_files) {

    $i++;
    $seq1_path = $indir.$ind_file_1;
    $j=0;

    # GET SEQ1 ID FROM FASTA FILE ID
    # IT is a good id to use short ids for this reason
    # See the DAWG-PAWS program fasta_shorten.pl for a script to do this
    if ($ind_file_1 =~ m/(.*)\.masked\.fasta$/) {
	$seq1_id = "$1";
    }
    elsif ($ind_file_1 =~ m/(.*)\.fasta$/ ) {	    
	$seq1_id = "$1";
    }
    elsif ($ind_file_1 =~ m/(.*)\.fa$/ ) {
	$seq1_id = "$1";
    }
    else {
	$seq1_id = $ind_file_1;
    }
    
    # Print the data to the VECOUT file
    print IDOUT "$seq1_id\n";

    # Short exit for testing
    #if ($i == "3") { exit };

    #-----------------------------+
    # FOR EVERY FILE IN DB DIR    |
    #-----------------------------+
    for my $ind_file_2(@fasta_files) {

	$j++;
 	$seq2_path = $indir.$ind_file_2;
	
	if ($ind_file_2 =~ m/(.*)\.masked\.fasta$/) {
	    $seq2_id = "$1";
	}
	elsif ($ind_file_2 =~ m/(.*)\.fasta$/ ) {	    
	    $seq2_id = "$1";
	}
	elsif ($ind_file_2 =~ m/(.*)\.fa$/ ) {
	    $seq2_id = "$1";
	}
	else {
	    $seq2_id = $ind_file_2;
	}
	
	print STDERR "Processing $i:$j\n" if $verbose;

	# GENERATE THE COMMAND FOR RUNNING FASTA
	my $cmd = "$program $seq1_path $seq2_path -B -m 9 -H -Q\n";

	# OPEN FILE HANDLE TO GET THE OUTPUT FROM PRSS
	open (ALIGN, "$cmd |") ||
	    die "Can not run the fasta program requested\n";

	# Parse the output
	my $in_align = 0;
	
	# $instr is the input string as read from the prss output
	while (my $instr = <ALIGN>) {
	    
	    # This will only get the first reported alignment
	    if ($in_align==1) {
		chomp ($instr);

		# The alignment string is below
		print "$instr\n" if $verbose;

		# Data has two columns
		my @data = split(/\t+/, $instr) ;
		
		# Test that data has two cols
		
		#-----------------------------+
		# FIRST SET OF VALS           |
		#-----------------------------+
		# The basic info set of data
		my @vals_1 = split(/\s+/, $data[0]) ; 
		my $num_vals_1 = @vals_1;

		my $z_score = $vals_1[$num_vals_1 - 2];
		my $e_val = $vals_1[$num_vals_1 - 1];

		#-----------------------------+
		# SECOND SET OF VALS          |
		#-----------------------------+
		# The extended values set of data
		my @vals_2 = split(/\s+/, $data[1]) ; 
		my $num_vals_2 = @vals_2;
		
		my $p_id = $vals_2[0];
		my $p_sim = $vals_2[1];
		my $sw = $vals_2[2];
		
		my $query_start = $vals_2[4];  # an0
		my $query_end = $vals_2[5];    # ax0

		#-----------------------------+
		# PRINT OUTPUT                |
		#-----------------------------+
                print STDOUT "$seq1_id\t";     # Use seq1 id from fasta file
		print STDOUT "$seq2_id\t";     # Use seq2 id from fasta file
		print STDOUT "$z_score\t";     # Z score
		print STDOUT "$e_val\t";       # E value
		print STDOUT "$p_id\t";        # Percent ID
		print STDOUT "$p_sim\t";       # Percent Similarity
		print STDOUT "$sw\t";          # Smith Waterman Score
		print STDOUT "$query_start\t"; # Start align in query seq
		print STDOUT "$query_end\n";   # End align in query seq

		#-----------------------------+
		# PRINT SIMILARITY OUTPUT     |
		#-----------------------------+
                print SIMOUT "$seq1_id\t";     # Use seq1 id from fasta file
		print SIMOUT "$seq2_id\t";     # Use seq2 id from fasta file
		print SIMOUT "$z_score\n";     # Z score, the similarity score

		#-----------------------------+
		# PUSH VALS TO Z SCORES ARRAY |
		#-----------------------------+
		push (@z_scores, $z_score);

		# No longer in info for best alignment
		$in_align=0;
	    }
	    

	    # Determine if we are in the alignment info
	    my @data = split (' ' , $instr);  # Split by white spaces
	    if ($data[0]) {
		if ($data[0] =~ "The") {
		    $in_align=1;
		    #print "$_";
		}
	    }
	    # End of determine if we are in the alignment info

	} # End of while ALIGN

    } # End of for each file in the database directory

    #-----------------------------+
    # GET THE MEDIAN Z VALUE      |
    #-----------------------------+
    # Using median PERL code from 
    # http://dada.perl.it/shootout/moments.perl.html
    
    @z_scores = sort { $a <=> $b } @z_scores;
    my $n = scalar(@z_scores);
    my $mid = int($n/2);
    my $z_median = ($n % 2) ? $z_scores[$mid] : 
	($z_scores[$mid] + $z_scores[$mid-1])/2;
    
    # Z median is sent out to the preferences file
    print STDERR "z_score median: $z_median\n";
    print PREFOUT "$z_median\n";

    # Empty the Z Value array
    @z_scores = ();

} # End of for each file in the query dir

#-----------------------------+
# CLOSE OPEN FILES            |
#-----------------------------+
close (SIMOUT);
close (IDOUT);
close (PREFOUT);

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


sub median {

}

=head1 NAME

batch_fasta.pl - Run the fasta program in batch mode.

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_fasta.pl -i indir -o OutFile

    --indir        # Path to the input file
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

The batch_fasta.pl program does not require a configuration file or 
depend on variables set in the user's environment.

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

UPDATED: 02/18/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/16/2008
# - Working version of the batch_prss34 program written
# - Uses a directory of single fasta files instead of a 
#   single file that is a multi fasta file.
# - Currently will just return the first reported local
#   alignment between two sequences
#
# 02/17/2008
# - Parse output to tab delim file
# - Modifed from running prss34 by default to running
#   any of the fasta programs. The program to be used can
#   now be specified by the -p,--program command line option
#   The default program is fasta34
#
# 02/18/2008
# - Working on getting the sequence ID from the fasta file name
# - Testing on large database
# - Cleaning up code after first SVN commit
# - Adding code to write vectors and matrix files
#
# TO DO:
# - Write to external file
# - Get multiple local alignments if desired
# - May need to switch to m10 output, m 9 is messy to
#   parse since the query name can have inappropriate charactersv
# - Use dbdir to allow for a database directory that is separate from
#   the query dir. This will make it easier to split the program
#   across nodes?? Would as
# - Should extract length of the query sequence and report this as well
# - May want to add ability to sort the fasta file as integer names
