#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# jaba_fasta.pl - Run all by all fasta on dir of seqs       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/13/2008                                       |
# UPDATED: 12/10/2008                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Given a directory of fasta filesrun the prss34 program to|
#  to generate a distance matrix. This will currently create|
#  output in the format needed for the apclust binary.      |
#  In the future this could be used to store data in the    |
#  BioSQL graph/tree extension.                             |
#                                                           |
# VERSION: $Rev$                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO:
# - Add option to select the edge attribute values that are
#   returned
# - add option to retrun the smith-waterman score as opposed
#   to the z-score. This can be done as a "score" option or
#   as a --z-score options
# - However the use of the --score option will allow
#    -$p_id ----> percent identity
#    -$p_sim ---> percent similarity
#    -$sw  -----> smith waterman
#    -$e_val ---> e value will need to be transformed
#    -$z_score -> default value used
#   --score
#        z -----> The z score, the default value used
#        pid ---> The percent id
#        psim --> The percent similarity
#        sw ----> The smith waterman score
#        nsw ---> The normalized smith waterman score
#        nle----> Negative log e value
#        e -----> The e value
#        all ---> Everything
# 
# - Upload the results as edge attributes to a database
#   however it may make sense to fetch the sequences from the database
#   as part of the process. 
# - Accept single fasta sequence as input, this could use
#   pipes with the seqio object to iterate across the db
# - The database size can be specified with the -Z command, this is
#   used for expectation value in FASTA, and SSEARCH. I think
#   that prss34 generates this using iterations
# - Write a sif file, this may be a problem, and would require threshold
#   values be passed at the command line.
# - --cytoscape option
#     - write the sif, na and EA files for cytoscape
# --db option
#    - upload the results as a network object to a db
# --ktup
#    - the ktuple argument MUST be the third command line argument
# - Allow for -i to pass a file, this will be split into a number of
#   fasta files in the input directory.


package REPMINER;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use Bio::SeqIO;
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
my $indir;
my $outdir;
my $dbdir;                     # Allows to split the job by query seq subsets
my $seq1_path;
my $seq1_id;
my $seq2_path;
my $seq2_id;
my $program = "prss34";
my $base_param = "-B -m 9 -H -Q"; # These are the base parameters
my $full_param;                   # The full paramater set passed to cmd
my $param;                        # Additional parameters passed to prog

# -B  ----> show the normalized score as a z score instead of a bitscore
# -m 9 ---> use the short aligment option
# -H -----> omit the histogram
# -Q -----> quiet, do not prompt for input

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
		    "param=s",     => \$param,
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

unless ($dbdir) {
    $dbdir = $indir;
}


#-----------------------------+
# BUILD THE PARAMETER STRING  |
#-----------------------------+

if ($param) {
    $param =~ s/^\s+//; #remove leading space
    $param =~ s/\s+$//; #remove trailing spaces
    $full_param = $base_param." ".$param;
}
else {
    $full_param = $base_param;
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
    die "ERROR: Can't open directory:\n$indir";
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
# GET FASTA FILES FROM DBDIR  |
#-----------------------------+
opendir( DIR, $dbdir ) ||
    die "Can't open directory:\n$dbdir";
my @db_fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_db_files = @db_fasta_files;

#-----------------------------+
# ERROR IF NOT FASTA FILES    |
#-----------------------------+
if ($count_db_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the db directory\n".
        "$dbdir\n".
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
    for my $ind_file_2(@db_fasta_files) {

	$j++;
 	$seq2_path = $dbdir.$ind_file_2;
	
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
	#my $cmd = "$program $seq1_path $seq2_path -B -m 9 -H -Q\n";
	my $cmd = "$program $seq1_path $seq2_path $full_param\n";
	
	# OPEN FILE HANDLE TO GET THE OUTPUT FROM FASTA
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
                print STDOUT "$seq1_id\t";     # 1 Use seq1 id from fasta file
		print STDOUT "$seq2_id\t";     # 2 Use seq2 id from fasta file
		print STDOUT "$z_score\t";     # 3 Z score
		print STDOUT "$e_val\t";       # 4 E value
		print STDOUT "$p_id\t";        # 5 Percent ID
		print STDOUT "$p_sim\t";       # 6 Percent Similarity
		print STDOUT "$sw\t";          # 7 Smith Waterman Score
		print STDOUT "$query_start\t"; # 8 Start align in query seq
		print STDOUT "$query_end\n";   # 9 End align in query seq

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
    # If a general approach for similarity values was used
    # the z value array would need to hold other types of values as well


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

sub median {

}

__END__;

=head1 NAME

jaba_fasta.pl - Run all by all fasta on dir of seqs

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    jaba_fasta.pl -i indir -o OutFile

=head2 Required Arguments

    -i       # Path to the dir containing the files to process
    -o       # Path to the output dir where the results will be placed

=head1 DESCRIPTION

Given a directory of fasta files, run the prss34 program to
to generate a distance matrix. This will currently create
output in the format needed for the apclust binary.
In the future this could be used to store data in the
BioSQL graph/tree extension. The current version will only run 
prss34, and will output the z score values. Future versions will incorporate
more options with respect to the fasta command line options and the 
scores that are returned. Full description of the fasta suite of 
programs is available in fasta3x.doc in the directory your fasta programs
are stored in. Currently all values will be printed to STOUT in the 
tab delim format. 

=over 2

=item 1. id of qry

This is derived from the file name be removing the fasta extension

=item 2. id of hit sequence

This is derived from the file name by removing the fasta extension.

=item 3. z score

=item 4. e value

=item 5. percent id

=item 6. Percent Similarity

=item 7. Smith Waterman Score

=item 8. Start algin in query sequence

=item 9. End align in query sequence

=back

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the input directory containing the fasta files to process.
These files must end with the 'fasta' or 'fa' file extension to be 
recognized. 

=item -o,--outdir

Path of the output directory where the results will be placed. If this output
directory does not exist, it will be created.  Three files will be placed in
this output directory:

=over 2

=item * similarity.txt

The three colmn data to generate the matrix of similarity values.

=item * vecname.txt

The labels for the vector describing the sequences used to generate
the similarity values.

=item * preferences.txt

The preferences file for use in the affinity propagation program.

=back

=back

=head1 OPTIONS

=over 2

=item -p,--program

The program from the fasta package to use. By default this is set to fasta34,
but prss34 could also be used. The valid options for this include:

=over 2

=item * fasta34

Generates a similarity score between two sequences. 

=item * prss34

This is slower then the fasta algorithm.  This generates a Smith-Waterman local
similarity score for the two sequences and tests for signficance using a
Monte-Carlo analysis. The result is a z-score based on permutations.

=item * ssearch34

Generates a Smith-Waterman score between two sequences. This algorithm is
about 10 times slower then fasta34, but is more sensitive for full length
comparisions.

=back

=item -d,--dbdir

The directory to use as the database directory. By default, the input dir
will serve as both the query sequence directory and the database directory.
Having a separate database directory allows to split the job by query
sequence sets that will all query the same db dir.

=item --param

A parameter string that will be passed to the FASTA program. 
This paramter string is the place to specifiy options such as the database 
size, or to change the score values used. 
By default the following paramters will
always be used '-B -m 9 -H -Q'. Any additional parameter options will
be appended to this string. For the full information on parameters
that can be modified, see the the fasta3x.doc file that was installed 
in your installation of the fasta package. 

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

=head1 EXAMPLES

=head2 Typical Use

To compare all the fasta sequences in a directory of files to themselves
you would run the following command

 jaba_fasta.pl -i indir/ -o outdir/

=head2 Send Scores to an External File

By default the scores that are produced by FASTA are printed to STDOUT.
This can be redirected to a file in Unix using >:

 jaba_fasta.pl -i indir/ -o outdir/ > score_file.txt

This file can be uploaded to a database or parsed to select scores other 
then the z score for use as a metric of similarity.

=head2 Use an Alternative FASTA program

By default the batch fasta program will use the prss34 program. This
option can be modified with the --program option. For example
the following command will use the fasta program to generate the
alignment score between sequences.

 jaba_fasta.pl -i indir/ -o outdir/ --program fasta34

=head2 Modify FASTA Parameters

It may become necessry for you to modify the paramters that you want
to use to run the fasta program. By default the following paramters will
always be used '-B -m 9 -H -Q'. Any additional parameter options will
be appended to this string. For example, to modifiy the match mismatch
scores with '-r' and add owercase masking '-S': 

 jaba_fasta.pl -i indir/ -o outdir/ --param '-r +3/-2 -S'

To specifgy the database size (in number of sequences) you would
use the '-Z' option. For example, if there are 12,000 sequences in
in the all by all fasta search:

 jaba_fasta.pl -i indir/ -o outdir/ -program fasta34 --param '-Z 12000'

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=head2 ERROR: Input/Output directory not specified at the command line

You are required to specifiy a directory for input using the '-i' option
and you must specify a directory for output using the '-o' option.

=head2 ERROR: Can't open directory DIRNAME

The directory that you specified may not exist. Check that the directory
does exist, and use the full directory path, for example if

 jaba_fasta.pl -i in_dir -o out_dir

fails to work try using full path to the loation of that file. For example,
if the directories are in your base home direcotry on Unix/Linux you should
try:

 jaba_fasta.pl -i /home/username/in_dir -o /home/username/out_dir

=head2 Can not find the program file

It is possible that the fasta programs are not in your user path, and you
will need to specify the full path to the fasta program that you want to 
use with the --program option. 
For example, if your prss34 program is in your home direcotry in
a subdirectory called apps ('/home/username/apps/'). You could try
the following command:

 jaba_fasta.pl -i in_dir -o out_dir --program /home/username/apps/prss34

=head1 CONFIGURATION AND ENVIRONMENT

The jaba_fasta.pl program does not require a configuration file or 
depend on variables set in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * FASTA

This program requires that the FASTA suite of programs be installed. More
information on the FASTA package is available from 
(http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml)
The package can be downloaded from:
http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml.

=back

=head2 Required Perl Modules

=over

=item * Bio::SearchIO

This module is part of bioperl and is required to parse BLAST
output in a format that is tiled across all HSPs.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the RepMiner
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

=over 2

=item Fasta Extensions

The fasta file extensions are currently limited to .fasta and .fa.
Files in the input directory that do not end in fasta or fa will be ignored.

=back

=head1 SEE ALSO

The jaba_fasta.pl program is part of the repminer package of 
repeat element annotation programs.
See the RepMiner web page 
( http://repminer.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/repminer/ )
for additional information about this package.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 02/13/2008

UPDATED: 12/10/2008

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
# 12/08/2008
# - Updating POD documentation
# - Added new print_help subfunction
#
# 12/10/2008
# - Added the option to use a sepearte db dir, this will useful
#   for running this on the cluster machine
# - Added the --param option, can specify parameter set to be sent
#   to the fasta program. This would be the place to specify the 
#   database size for the fasta34 and 
#   The paramter option also allows for changing tupule size,
#   scoring, and db size to best fit the data at hand.
# - Changed default program from fasta34 to prss34
#
# TO DO:
# - Get multiple local alignments if desired
# - May need to switch to m10 output, m 9 is messy to
#   parse since the query name can have inappropriate characters
# - Use dbdir to allow for a database directory that is separate from
#   the query dir. This will make it easier to split the program
#   across nodes?? Would as
# - Should extract length of the query sequence and report this as well
# - May want to add ability to sort the fasta file as integer names
