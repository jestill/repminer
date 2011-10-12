#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_ltr_annotate.pl - Annotate a dir of fasta files     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 04/04/2008                                       |
# UPDATED: 04/07/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a direcotry of fasta sequences representing full   |
#  length copies of LTR retrotransposons, this will         |
#  annotate the biological regions of the LTR retro.        |
#  to annotate LTR retrotransposons. Results are stored     |
#  in a table or database.                                  |
#                                                           |
# USAGE:                                                    |
#  batch_ltr_annotate.pl 
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
use Bio::SeqIO;                # Input/Output for sequence files
use Getopt::Long;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep
		    clock_gettime clock_getres clock
		    stat );

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                    # Input directory where fasta fiels are stored
my $outdir;                   # Base directory where output files will be placed
my $name_root;                # Root name of output files
my $count_files;              # Number of files in the input dir
my $count_db;                 # Number of query databases
my $count_proc;               # Total number of procedures ...
my $file_num;

# Booleans
my $draw_image =0;
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Default values
my $max_e_val = 0.00001;

my $ltr_feats;                 # Hash reference to ltr_features for qry seq 

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s" => \$outdir,
		    # ADDITIONAL OPTIONS
		    "p|draw"      => \$draw_image, # Draws an image, p for png
		    "e|e_val"     => \$max_e_val,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

#-----------------------------+
# CHECK FOR REQUIRED OPTIONS  |
#-----------------------------+
# infile and outfile required

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
#if ( (!$indir) || (!$outdir) || (!$file_config) || (!$dir_blast_db) ) {
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
#    print STDERR "ERROR: A configuration file was not specified at the".
#	" command line\n" if (!$file_config);
#    print STDERR "ERROR: A BLAST DB directory was not specified at the".
#	" command line\n" if (!$dir_blast_db);
#    print_help ("usage", $0 );
}

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

my $time_start = Time::HiRes::time();

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

$count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------+
# BLAST DATABASES             |
#-----------------------------+
# The following will need to be changed to a config file that uses

my @blast_dbs = ("/db/jlblab/pfam/gag_poaceae",
		 "/db/jlblab/pfam/zf_cchc_poaceae",
		 "/db/jlblab/pfam/rvp_all_poa",
#		 "/db/jlblab/pfam/rvp_poaceae",  # May not be all
		 "/db/jlblab/pfam/rve_poaceae",
		 "/db/jlblab/pfam/rvt_poaceae", # both rvt_1 and rvt_2
		 "/db/jlblab/pfam/rh_poaceae",
		 "/db/jlblab/pfam/chromo_viridiplant",
		 "/db/jlblab/pfam/dros_env",
		 );

$count_db = @blast_dbs;
print STDERR "NUMBER OF DATABASES: $count_db\n" if $verbose;
$count_proc = $count_db * $count_files;
print STDERR "NUMBER OF PROCESSES: $count_proc\n" if $verbose;


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
# Open output files for each database ??
# This would require that I decide the output file based on the 
# feat name that is returned
open (SEQOUT, ">".$outdir."feat_seq_out.fasta") ||
    die "Can not open feature sequence outfile\n";

open (FEATSUM, ">".$outdir."feat_summary.txt") ||
    die "Can not open the feature summary file\n";

# Write feature summary header
print FEATSUM 
    "seq_id\t".                   # Sequence ID
    "seq_len\t".
    "gag_s\tgag_e\t".             # GAG
    "zf_s\tzf_e\t".               # Zinc finger, gag component
    "rvp_s\trvp_e\t".             # Protease/AP
    "rve_s\trve_e\t".             # INTEGRASE
    "rvt_s\trvt_e\t".             # Reverse Transcriptase
    "rh_s\trh_e\t".               # RnaseH
    "chrom_s\tchrom_e\t".         # Chromodomain
    "env_s\tenv_e\n";             # Envelope

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
for my $ind_file (@fasta_files)
{
    $file_num++;

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }
    
    my $infile = $indir.$ind_file;
    my $seqin = Bio::SeqIO->new(-file   => "<$infile",
				-format => 'fasta');
    
    # The following format will also work for a directory 
    # of multiple fasta files

    # For each sequence record in the fasta file
    while (my $inseq = $seqin->next_seq) {


	# Print seq_id to the summary file
	print FEATSUM $inseq->primary_id()."\t";
	print FEATSUM $inseq->length()."\t";

        #-----------------------------+
        # PROCESS THE FEATURES FOR    |
        # EACH BLAST DATABASE         |
        #-----------------------------+     
	foreach my $blastdb (@blast_dbs) {
	    
	    print "Processing:\t".$blastdb."\n" if $verbose;
	    
	    my $seq_id = $inseq->primary_id();
	    my $seq_string = $inseq->seq();

	    my $annotated_rpn = seq_annotate ($seq_id, $seq_string, 
					      $blastdb, $max_e_val );

	    if ($annotated_rpn) {
		# If the feature has a value then print to STDOUT
		print STDOUT "FEAT NAME:\t".
		    $annotated_rpn->{ feat_name }."\n";
		print STDOUT "    START:\t".
		    $annotated_rpn->{ feat_start }."\n";
		print STDOUT "      END:\t".
		    $annotated_rpn->{ feat_end }."\n";
		print STDOUT "    SEQID:\t".
		    $annotated_rpn->{ qry_id }."\n";
		print STDOUT " FEAT_SEQ:\t".
		    $annotated_rpn->{ feat_seq }."\n";
		print STDOUT "\n";
		
		# PRINT SEQUENCE OUTPUT TO FEATURE FILE
		print SEQOUT ">".$annotated_rpn->{ qry_id }."_".
		    $annotated_rpn->{ feat_name }."\n";
		print SEQOUT $annotated_rpn->{ feat_seq }."\n";
		
		# PRINT 
		print FEATSUM $annotated_rpn->{ feat_start }."\t";
		print FEATSUM $annotated_rpn->{ feat_end }."\t";
	    }
	    else {
		# PRINT NULL VALUES TO FEATURE SUMMARY FILE
		# Currently using N to indicate Null. Not present
		print FEATSUM "N\tN\t";
	    }
	    
	} # End of for each blast database

	#-----------------------------+
	# new line to the sum file    |
	#-----------------------------+
	print FEATSUM "\n";	

	#-----------------------------+
	# CAN PRINT IMAGE OUTPUT HERE |
	#-----------------------------+
	# Need temp array to store feature locations
	
	# Reset array vals to null


    } # End of for each sequence record


} # End of for each sequence file

my $time_end = Time::HiRes::time();
my $time_to_run = $time_end - $time_start;

# Print time to run if in verbose mode
print STDERR "TIME: $time_to_run\n" if $verbose;

# CLOSE FILE HANDLES
close (SEQOUT);
close (FEATSUM);

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


sub seq_annotate {

    my $max_num_hits = 1;
    my @ans;              # Answer
    my $vals;             # Hash Reference to values
    my $feat_name;

    # annotate a seq from a datbase
    my ($seq_id, $seq_string, $blastdb, $max_e_val ) = @_;
    # dbh is the database handle where the data are to be store
    # seq_id is the id of the query sequence
    # seq_string is the
    # b should set the number of alignments returned
    my @bl_params = ('b'       => 1,
		     'e-value' => $max_e_val,
		     'program' => 'blastx', 
		     'database' => $blastdb);
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@bl_params);
    
    my $qry_input = Bio::Seq->new(-id=> $seq_id,
				  -seq=>$seq_string );
    
    my $blast_report = $factory->blastall($qry_input);

    # This currently assumes a single query sequence was used
    my $hit_count = 0;
    
    while (my $blast_result = $blast_report->next_result()) {

	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		if ($hit_count < $max_num_hits) {
		    my ($feat_start, $feat_end) = $blast_hsp->range('query');

		    # Print to STDERR IF VERBOSE
#		    print STDERR $feat_name."\n" if $verbose;    
#		    print $feat_start."\n" if $verbose;
#		    print $feat_end."\n" if $verbose;
#		    print STDERR
#			$blast_result->query_name."\n" if $verbose;
#		    print STDERR 
#			$blast_hsp->query_string."\n" if $verbose;
		    $hit_count++;
		    #my @range = $blast_hsp->range('query');

		    # Load values to the hash reference
		    #$vals->{'feat_name'} = $feat_name;
		    $vals->{'feat_name'} = $blast_result->database_name;
		    $vals->{'feat_start'} = $feat_start;
		    $vals->{'feat_end'} = $feat_end;
		    $vals->{'qry_id'} = $blast_result->query_name;
		    $vals->{'feat_seq'} = $blast_hsp->query_string;

		} # End of hit_count
	    }
	}
    } # End of while blast_result

    # Return the feature annotation
    if ($hit_count > 0) {
	return $vals;
    }
    else {
	return 0;
    }

}


=head1 NAME

ltr_annotate.pl - Use blast to Annotate LTR retrotransposon

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

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 04/11/2008
# - Added summary output fil3
# - Adding ability to draw cartoon of ltr internal structures
