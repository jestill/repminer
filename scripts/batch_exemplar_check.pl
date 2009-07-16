#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_exemplar_check.pl - All-by-all check of exemplars   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# NOTE:                                                     |
# Expects the input header to be in the RepeatMasker format |
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/04/2008                                       |
# UPDATED: 07/16/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of fasta files, do an bl2seq for all   |
#  unique pairs of sequences, report sequences above a      |
#  threshold value as well list the similarities for all    |
#  sequences in the input directory. This reports both      |
#  the tiled bitscore as well as the ratio of the bitscore  |
#  of the pair of sequences divided by smaller bitscore     |
#  of either sequence blasted against itself.               |
#                                                           |
# USAGE:                                                    |
#  batch_bl2seq.pl -i indir  -o summary_file.txt            |
#                                                           |
# VERSION: $Rev: 435 $                                      |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# Notes:
#  - the index for the matrix holding the results of the of 
#    pairwise blast starts at 1. This leaves the [i,0] and 
#    [0,j] to serve as labels.
#  - This currently does the blast for softmasked sequences
#    by default.

package DAWGPAWS;

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
#use Bio::Tools;
use Bio::SearchIO;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 435 $ =~ /(\d+)/;

# VARS WITH GLOBAL SCOPE
my $indir;                     # Dir holding the fasta files for queries
my $dbdir;                     # Dir holding the fasta files for hits (db)
my $outfile;                   # File path to the report file
my @bitscore_matrix;           # Two dim array to hold bitscore values
my @ratio_matrix;              # Two dim array to hold normalize values


# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;               # Run the program in test mode
my $do_clean = 0;              # Remove the bln output files
my $blastdir;                  # The dir the blast reports will be
                               # moved to when finished.

# VARS WITH DEFAULT VALS
my $thresh_pid = 0.8;          # The threshold percent id
my $thresh_len = 80;           # The threshold length

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
		    "d|dbdir=s"   => \$dbdir,
		    # ADDITIONAL OPTIONS
                    "o|outfile=s" => \$outfile,
		    "blastdir=s"  => \$blastdir,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "clean"       => \$do_clean,
		    "thresh-pid=s"    => \$thresh_pid,
		    "thresh-len=s"    => \$thresh_len,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
# Output file is not required, will print to STDERR if not specified
if ( (!$indir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print_help ("usage", $0 );
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

if ($dbdir) {
    unless ($dbdir =~ /\/$/ ) {
	$dbdir = $dbdir."/";
    }
}
else {
    $dbdir = $indir;
}

if ($blastdir) {
    unless ($blastdir =~ /\/$/ ) {
	$blastdir = $blastdir."/";
    }
}

# Print to outfile if one specified, otherwise to STDOUT
# This is a report output file
if ($outfile) {
    open (REPOUT, ">$outfile");
} 
else {
    open (REPOUT, ">&STDOUT");
}

# PRINT THE HEADER FOR THE TAB DELIM TEXT FILE
print REPOUT "CONTAM\tQUERY\tHIT\tOBS_PID\tSTART\tEND\tLENGTH";
print REPOUT "\n";

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @in_fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $count_in_files = @in_fasta_files;

# Error if no files in input directory
if ($count_in_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}
print STDERR "NUMBER OF IN FILES TO PROCESS: $count_in_files\n" 
    if $verbose;


#-----------------------------+
# Get the FASTA file from the
# directory provided by the
# var $dbdir.
#-----------------------------+
opendir( DIR, $dbdir ) || 
    die "Can't open directory:\n$dbdir"; 
my @db_fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $count_db_files = @db_fasta_files;

#-----------------------------+
# ERROR IF NO FILES FOUND IN  |
# THE DATABASE DIRECTORY      |
#-----------------------------+
if ($count_db_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the db directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}






my $num_square = $count_in_files * $count_db_files;
#my $total_proc = ($num_square - $count_files)/2;
#my $total_proc = (($num_square - $count_files)/2) + $count_files;

my $total_proc = $num_square;

print STDERR "TOTAL BL2SEQ: $total_proc\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
# Does all by all bl2seq for j<j
# Therefore this does not do the self blast or the reciprocals

#-----------------------------+
# RUN THE BLAST 2 SEQ         |
#-----------------------------+
# I will include self hits to put i/j in perspective
my $i_num = 0;
my $proc_num;

# Get all i less then j
# fill the matrix of bitscore with the diagonals

for my $file_i (@in_fasta_files) {
    
    $i_num++;
    
    my $j_num = 0;
    for my $file_j (@db_fasta_files) {
	
	$j_num++;
	$proc_num++;
	
	my $bln_outfile = $file_i."_".$file_j.".bln";
	
	# The following command includes lowercase filtering
	my $bl_cmd = "bl2seq -p blastn -U".
	    " -i ".$indir.$file_i.
	    " -j ".$dbdir.$file_j.
	    " -o $bln_outfile";
	
	print STDERR "Proc $proc_num of $total_proc -- $i_num : $j_num\n"
	    if $verbose;
	
#	    # Only print the following for dbug of blast cmd line
#	    #print STDERR "\t$bl_cmd\n" if $verbose;
	
	# Run the BLAST
	# with Bio:SearchIO

	#////////////////////////////////////////////////////////
	# NOTE: It will also be possible here to run cross_match
	#       if desired and then parse/clean with cross_match
	#       it may be best to first identify potential 
	#       problems with blast and then move this to a
        #       dir and process these problems separately. 
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	system ($bl_cmd);
	
	# Parse the results of the bl2seqoutput
	my $has_contam;
	if (-e $bln_outfile) {
	    
	    # Do a check for has contam
	    $has_contam = &blseq2contam ($bln_outfile, 
					   $file_i, 
					   $file_j, 
					   $outfile, 1 );
	    
	    # Print the output to STDERR
#	    print STDERR "$file_i\t$file_j\t$score_val\n" if $verbose;
	    
#	    # Load the result to the score matrix
#	    # Use zero if theree is not a score value
#	    if ($score_val =~ "NULL" ) {
#		$bitscore_matrix[$i_num][$j_num] = 0;
#	    }
#	    else {
#		$has_contam = 0;
#		$bitscore_matrix[$i_num][$j_num] = int($score_val);
#	    }
	    
	    
	    #-----------------------------+
	    # THE do_clean OPTION WILL    |
	    # GENERATE A CLEAN COPY       |
	    # OF THE EXEMPLAR             |
	    # THE CLEAN FUNCTION WILL     |
	    # PROBABLY BE WRITTEN AS      |
	    # A SUBFUNCTION               |
	    #-----------------------------+
	    if ($do_clean) {
		
	    }

#	$out_file_path = $outdir.$seq_id.".fasta";
	    
	    #-----------------------------+
	    # MOVE THE BLAST RESULT TO    |
	    # THE BLAST DIRECTORY         |
	    #-----------------------------+
	    if ($blastdir) {
		my $move_cmd = "mv $bln_outfile $blastdir";
		system ($move_cmd);
	    }
	    
	} # End of process the bl2seq output
	
	
    }
    
}

#-----------------------------+
# CLOSE OUT PROGRAM           |
#-----------------------------+
# Close the file handle
close REPOUT;

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

sub blseq2score {
#-----------------------------+
# The purpose of this is to   |
# return the tiled bitscore   |
# for the bl2seq result       |
#-----------------------------+

    # blastin - path to the blast input file
    # seq_1   - First sequence in the comparision
    # seq_2   - Second sequence in the comparision
    # tabout  - path to the gff output file
    # append  - boolean append data to existing tab text filetab out

    my ($blastin, $seq_1, $seq_2, $tabout, $append ) = @_;    
    my $bitscore;
    
    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
 					   '-file'   => $blastin) 
 	|| die "Could not open BLAST input file:\n$blastin.\n";
    
    while (my $blast_result = $blast_report->next_result()) {
	
    	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    $bitscore = $blast_hit->bits;
	    
	}
    }

    return $bitscore || "0";

}



sub blseq2contam {
#-----------------------------+
# The purpose of this is to   |
# identified potential cross  |
# contamination in exemplars  |
#-----------------------------+

    # blastin - path to the blast input file
    # seq_1   - First sequence in the comparision
    # seq_2   - Second sequence in the comparision
    # tabout  - path to the gff output file
    # append  - boolean append data to existing tab text filetab out

    my ($blastin, $seq_1, $seq_2, $tabout, $append ) = @_;    
    my $bitscore;
    my $found_contam = 0;  # Was contamination found
    
    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
 					   '-file'   => $blastin) 
 	|| die "Could not open BLAST input file:\n$blastin.\n";
    
    while (my $blast_result = $blast_report->next_result()) {

	# GET THE QUERY FAMILY 
	my $query_name = $blast_result->query_name();
	my $query_family = "unknown";
	if ($query_name =~ m/(.*)\/(.*)/ ) {
	    $query_family = $2;
	}
	
    	while (my $blast_hit = $blast_result->next_hit()) {

	    my $hit_name = $blast_hit->name();
	    my $hit_family = "unknown";
	    if ( $hit_name =~ m/(.*)\/(.*)/ ) {
		$hit_family = $2;
	    }
	    # If unkown just show file name
	    if ($hit_family =~ "unknown") {
		$hit_family = $hit_name;
	    }
	    
#	    print STDERR "$query_name : $hit_name \n";
#	    print STDERR "PROC BLAST $query_family - $hit_family\n" 
#		if $verbose;

	    if ( $hit_family ne $query_family) {
	#	print STDERR "--->Processing $query_family - ".
	#	    "$hit_family\n" if $verbose;
		while (my $blast_hsp = $blast_hit->next_hsp()) {

		    my $obs_pid = $blast_hsp->frac_identical;
		    my $obs_len = $blast_hsp->length('total');
		    
		    if ($obs_len > $thresh_len &&
			$obs_pid > $thresh_pid ) {
			# Using a global REPOUT reprot output file handle
			# here, bad practice. Don't do this at home.
			# The contam statement says who is potentially
			# contamination your query sequence

			my $data_out = 
			    "$hit_family-->$query_family\t".
			    "$query_name\t".
			    "$hit_name\t".
			    "$obs_pid\t".
			    $blast_hsp->start('query')."\t".
			    $blast_hsp->end('query')."\t".
			    "$obs_len\t".
			    "\n";
			print REPOUT $data_out;
			
			# Print positive hits to stderr unless
			# they are already coming to STDOUT
			if ($outfile) {
			    print STDERR $data_out if $verbose;
			}
			
		    } # End of if obs_len > threshold_len
		    
		}
	    }
	    
	}
    }
    

#    return $bitscore || "0";
    # Initially will report if contamination was found
    return $found_contam;

}


1;
__END__

#-----------------------------------------------------------+
# JUNKYARD                                                  |
#-----------------------------------------------------------+


#-----------------------------+
# LOAD NORMAILZED BITSCORE    |
# RATIO TO THE RAIO MATRIX    |
#-----------------------------+
# Determine a normailized bitscore ratio
# using the max (self_i self_j) as the denominator
# with i compared to j as the numerator
#
# This will also take the time to print the values
# above the threshold value

print REPOUT "#-----------------------------------------------------------\n";
print REPOUT " VALS ABOVE THRESHOLD\n";
print REPOUT "#-----------------------------------------------------------\n";

my $i = 0;
my $thresh_count = 0;

for my $file_i (@fasta_files) {
    $i++;

    my $j = 0;

    for my $file_j (@fasta_files) {

	$j++;

	if ($i < $j) {    
	    
	    my $i_self = $bitscore_matrix[$i][$i];
	    my $j_self = $bitscore_matrix[$j][$j];
	    
	    if ( $i_self > $j_self ) {
		$ratio_matrix[$i][$j] = $bitscore_matrix[$i][$j]/$i_self;
	    }
	    else {
		$ratio_matrix[$i][$j] = $bitscore_matrix[$i][$j]/$j_self;
	    }    

	    # REPORT VALUES ABOVE THRESHOLD
	    if ($ratio_matrix[$i][$j] > $thresh_val) {
		$thresh_count++;
		print REPOUT "$file_i\t$file_j\t".
		    $ratio_matrix[$i][$j]."\t".
		    $bitscore_matrix[$i][$j]."\n";
	    }

	} # End of if i < j

    } # End of for file_j
} # End of for file_i

# Print spacer after values above threshold
if ($thresh_count == 0) {
    print REPOUT "No values are above the threshold value of $thresh_val\n";
    print REPOUT "\n";
} 
else {
    print REPOUT "$thresh_count values are above the threshold value of ".
	"$thresh_val";
    print REPOUT "\n";
}

#-----------------------------+
# PRINT RATIO MATRIX TO THE   |
# REPORT OUTPUT FILE          |
#-----------------------------+
# This is as a three column output file
print REPOUT "#-----------------------------------------------------------\n";
print REPOUT "# FULL RATIO DATA MATRIX\n";
print REPOUT "#-----------------------------------------------------------\n";
$i = 0;
for my $file_i (@fasta_files) {
    $i++;

    my $j = 0;

    for my $file_j (@fasta_files) {

	$j++;
	
	if ($i < $j) {    
	    
	    my $round_val = sprintf( "%.4f", $ratio_matrix[$i][$j] );
	    print REPOUT "$file_i\t$file_j\t$round_val\t".
		$bitscore_matrix[$i][$j]."\n";
	    
	} # End of if i < j

    } # End of for file_j
} # End of for file_i



=head1 NAME

batch_bl2seq.pl - Batch bl2seq to compare assemblies

=head1 VERSION

This documentation refers to program version $Rev: 435 $

=head1 SYNOPSIS

  USAGE:
    batch_bl2seq.pl -i InDir -o outfile

    --indir         # Dir containing the fasta files
    --outfile       # Output reoport file file

=head1 DESCRIPTION

The purpose of this script is to try to find clones that are
replicates of one another.
Given a directory of fasta files, do an bl2seq for all
unique pairs of sequences, report sequences above a
threshold value as well list the similarities for all
sequences in the input directory. This reports both
the tiled bitscore as well as the ratio of the bitscore
of the pair of sequences divided by smaller bitscore
of either sequence blasted against itself. If an output
file is not specified at the command line, the output will
be printed to STDOUT.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.
It is expected that these will be single fasta files that 
represent ordered contig assemblies.

=back

=head1 OPTIONS

=over 2

=item -o,--outfile

The report file that is generated. The first section of this 
file reports the sequence pairs that are above the threshold value, the 
remainder of this report will be all the remainder sequence pairs.

=item --thresh

This is the threshold ratio value to report at the header of the
output file. This should be a value between zero and one. The default
value is 0.80. Thus with the default value, all sequences which have 
a bitscore ratio above 80% of the shorter of the two squences will be 
reported.

=item --clean

Specifying clean at the command line will remove the bl2seq reports from
the working directory. By default the bl2seq output files will be left
in the directory.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

This program does not make use of a configuration file. All options
are specified at the command line.

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * bl2seq

The bl2seq command line program is part of the NCBI stand alone blast package. 
If you typle bl2seq at the command line, you should see the options for running
the bl2seq program. 
The executables for stand alone blast are available at
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST

=back

=head2 Required Perl Modules

=over BioPerl

This program relies on the BioPerl package of programs. Information on
downloading and installing bioperl is available at 
http://www.bioperl.org/wiki/Getting_BioPerl.

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Not tested on OS other then Linux.

Although this program should run on either Windows or Mac OSX, it has only
been tested on the linux OS.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 10/04/2008

UPDATED: 10/06/2008

VERSION: $Rev: 435 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 10/04/2008
#  - Program started
#  - Wrote blseq2score subfunction
#
# 10/06/2008
#  - Added POD documentation
#
# 07/16/2009
#  - Modified from an existing batch_bl2seq.pl program
#  - POD docs currently reflect the old program this 
#    was copied from
