#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2sim.pl - Convert BLAST report to similarity file|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 02/24/2008                                       |
# UPDATED: 12/06/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
# Extract blast -m8 output to a simple three column         |
# similarity file. This output file is a tab delimited      |
# text file that can be easily imported into a SQL database |
# or used as an import file for clustering.                 |
# As part of the process this script also parses the        |
# RepMiner header format and returns the unique identifier  |
# number from the fasta header.                             |
#                                                           |
# EXAMPLES:                                                 |
# cnv_blast2sim.pl -i infile -o outfile.txt                 |
#                                                           |
#-----------------------------------------------------------+
# TODO:
#  Making parsing the fasta header to integer an option.
#  The advantage of using a simple integer is that this 
#  reduces the size of the file needed to store the results.
# TO DO:
#  Transform the significance value to a value that can be used for
#  similarity.
#
# FOR USE IN SQL DATABASES
# FIRST:
# cnv_blast2sim.pl my_blast_out.blo > ltr_ltr.txt
# THEN CREATE TABLE TO HOLD THIS DATA:
# mysql> create table tbl_ltr_ltr ( 
#                qry INT, 
#                sub INT, 
#                bit INT, 
#                INDEX (qry), 
#                INDEX (sub));
# It would be better to use DECIMAL for the bitscore
# mysql> create table tbl_ltr_ltr (
#                qry INT,
#                sub INT,
#                bit DECIMAL (7,2),
#                INDEX (qry),
#                INDEX (sub));

# THEN DO THE FOLLOWING AS ROOT:
# cp ltr_ltr.txt /var/lib/mysql/db_maize_ltr_2008/
# THEN DO
# mysql> LOAD DATA INFILE 'ltr_ltr.txt' INTO TABLE tbl_ltr_ltr;
# THIS WILL RESULT IN THE DATA BEING AVAILABLE IN THE MYSQL DATABASE

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # 
use Getopt::Long;              # Get options from command line
use Bio::SearchIO;             # Parse BLAST output
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
#my $blast_opt = 2;             # BLAST option, how to parse blast
                               # Use numbers for now
                               # 1 --> Bitscore of Best HSP from -m8 BLAST
                               # 2 --> Bitscore of Tiled HSPs from -m8 BLAST

# BOOLEANS
my $verbose = 0;
my $line_num = 0;
my $quiet = 0;
my $show_usage = 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0; 
my $blast_opt = "0";

# Index Vals
my $pre_cat_id ="0";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
#my $infile = shift;
#my $outfile = shift;

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    "b|blast-opt=s" => \$blast_opt,  # Has default 0
                    # ADDITIONAL OPTIONS
                    "q|quiet"       => \$quiet,
                    "verbose"       => \$verbose,
                    # ADDITIONAL INFORMATION
                    "usage"         => \$show_usage,
                    "version"       => \$show_version,
                    "man"           => \$show_man,
                    "h|help"        => \$show_help,);


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

#-----------------------------+
# CHECK REQUIRED VARIABLES    |
#-----------------------------+
# Currrently not variables are required.

#-----------------------------------------------------------+
# MAIN BODY                                                 |
#-----------------------------------------------------------+

# Print to outfile if one specified, otherwise to STDOUT
if ($outfile) {
    open (SIMOUT, ">$outfile");
} 
else {
    open (SIMOUT, ">&STDOUT");
}

#-----------------------------------------------------------+
# BITSCORE OF BEST HSP : OPTION 1                           |
#-----------------------------------------------------------+
# This option is not dependent on bioperl
if ($blast_opt =~ "1" ) {

    print STDERR "Parsing Best HSP\n" if $verbose;

    open (BLASTIN, "<$infile") ||
	die "Can not open $infile\n";
    
    while (<BLASTIN>) {
	
	$line_num++;
	print STDERR "LINE: $line_num\n" if $verbose;
	chomp;
	
	# Ignore comment lines
	unless (m/^\#.*/) { 
	    
	    my ($qry_id, $sub_id, $pid, $alen, 
		$mismatch, $gap_open, 
		$qry_start,$qry_end, $sub_start, $sub_end,
		$e_value, $bit_score) = split(/\t/);
	    
	    # Trim leading white space from bit score
	    $bit_score =~ s/^\s*(.*?)\s*$/$1/;
	    
	    my ($x_val, $cruft) = split(/\|/, $qry_id);
	    my ($y_val, $more_cruft) = split(/\|/, $sub_id);
	    
	    my $x_int = int($x_val);
	    my $y_int = int($y_val);
	    
	    my $cur_cat_id = $x_val.":".$y_val;
	    
	    # The following just grabs the best hit data
	    unless ( $cur_cat_id =~ $pre_cat_id) {
		print STDOUT "$x_val\t$y_val\t$bit_score\n";
		$pre_cat_id = $cur_cat_id;
	    } # End of grab best hit
	}
    } # END OF WHILE BLAST IN
} # End of if BLAST Option = 0

#-----------------------------------------------------------+
# TILED BITSCORE : OPTION 2                                 |
#-----------------------------------------------------------+
elsif ( $blast_opt =~ "2" ) {
    
    my $blast_report;
    my $count_result = 0;
    print "Parsing tiled HSPs\n" if $verbose;
    

    # OPEN STDIN OR FILE PATH
    if (!$infile) {
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-fh'   => \*STDIN )
	    || die "Could not open BLAST input from STDIN.\n";
    }
    else {
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-file'   => $infile)
	    || die "Could not open BLAST input file:\n$infile.\n";
    }


    while (my $blast_result = $blast_report->next_result()) {
	
	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    # Do something with the BLAST data
	    my ($x_val, $cruft) = split(/\|/, $blast_result->query_name);
            my ($y_val, $more_cruft) = split(/\|/, $blast_hit->name);

            my $x_int = int($x_val);
            my $y_int = int($y_val);

	    # PRINTING RAW NAMES
	    #print STDOUT $blast_result->query_name."\t";
	    #print STDOUT $blast_hit->name."\t";

#	    # PRINTING PARSED NAMES
#	    print STDOUT "$x_val\t$y_val\t";
#	    print STDOUT $blast_hit->bits()."\n";

	    # 07/15/2008
	    # Changing to the following
	    # PRINTING PARSED NAMES
	    print SIMOUT "$x_val\t$y_val\t";
	    print SIMOUT $blast_hit->bits()."\n";
	    
	} # End of next BLAST hit
	
    } # End of next BLAST result
}

#-----------------------------+
# TILED SIGNIFICANCE : OPT 4  |
#-----------------------------+
elsif ( $blast_opt =~ "4"  ) {
    
    my $blast_report;
    my $count_result = 0;
    print "Parsing tiled HSPs\n" if $verbose;
    
    # Open the BLAST report object
    # blasttable for -m8 output
   # OPEN STDIN OR FILE PATH
    if (!$infile) {
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-fh'   => \*STDIN )
	    || die "Could not open BLAST input from STDIN.\n";
    }
    else {
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-file'   => $infile)
	    || die "Could not open BLAST input file:\n$infile.\n";
    }


    while (my $blast_result = $blast_report->next_result()) {
	
	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    # Do something with the BLAST data
	    my ($x_val, $cruft) = split(/\|/, $blast_result->query_name);
            my ($y_val, $more_cruft) = split(/\|/, $blast_hit->name);

            my $x_int = int($x_val);
            my $y_int = int($y_val);

	    # PRINTING RAW NAMES
	    #print STDOUT $blast_result->query_name."\t";
	    #print STDOUT $blast_hit->name."\t";

	    # PRINTING PARSED NAMES
	    print SIMOUT "$x_val\t$y_val\t";
	    print SIMOUT $blast_hit->significance()."\n";
	    
	} # End of next BLAST hit

    } # End of next BLAST result
} else {
    print "\a";
    print STDERR "A valid BLAST parsing option was not passed at "
	."the command line\n";
}

# END OF PROGRAM
close (SIMOUT);
exit 0;

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

cnv_blast2sim.pl - Convert BLAST report to similarity file

=head1 VERSION

This documentation refers to cnv_blast2sim.pl version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_blast2sim.pl -i BlastOutput.bln -o SimFile.txt

=head2 Required Arguments

    -i, --infile    # Path to the input file to parse
    -o, --outfile   # Path to the output similarity file
    -b, --blast-opt # Blast data to use for similairty metric 

=head1 DESCRIPTION

Extract blast -m8 or -m9 output to a simple three column
similarity file. This output file is a tab delimited
text file that can be easily imported into a SQL database
or used as an input file for clustering.
As part of the process this script also parses the
RepMiner header format and returns the unique identifier
number from the fasta header.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the BLAST file to be parsed. If no infile argument is used
the program will attempt to parse from STDIN.

=item -o,--outfile

Path of the output file. If not outfile arugment is used, the program
will send data to standard output.

=item -b, --blast-opt

Number indicating the option used to extract a similarity
score from a BLAST result. Valid arguments are:

=over

=item 1.

Use the BITSCORE of the best HSP.

=item 2.

Use a tiled BITSCORE value. This is the default option.

=item 3.

Use the significance value of the best HSP. 
THIS OPTION IS NOT CURRENTLY IMPLEMENTED.

=item 4.

Use a tiled significance value

=back

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

Run the program without doing the system commands.

=back

=head1 EXAMPLES

If you have a BLAST file in the -m 8 aligment format, this can be used as input
for the cnv_blast2sim.pl program. The following example shows converting
a BLAST report (blast_result.bln) to a similarity matrix (sim_matrix.txt):

 cnv_blast2sim.pl -i blast_result.bln -o sim_matrix.txt

By default, this will use a tiled bitscore as the metric of similarity.
You may also select to report just the bitscore for the best HSP from
an aligment between two sequence records using the -b 1 options. 
Using the same blast_result.bln report this command would be:

 cnv_blast2sim.pl -i blast_result.bln -o sim_matrix.txt -b 1

The output from BLAST can be passed directly to this parser using pipes. 
The following example shows using the cnv_blast2sim.pl with default parameters.
The output will be printed to STDOUT. Using a fasta file (seqs.fasta) and the 
BLAST formatted database of this file (myDB), the command to run BLAST using
blastn follows:

 blastall -m 8 -p blastn -i seqs.fasta -d myDB | cnv_blast2sim.pl

The output file path (sim.txt) can be specified with the -o option as shown
below:

 blastall -m 8 -p blastn -i seqs.fasta -d myDB | cnv_blast2sim.pl -o sim.txt

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item Program does not appear to be doing anything

It is possible that the program is waiting for input from STDIN.
Run with the --verbose mode.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This script does not make use of external configuration files
or any options in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * NCBI blastall

The latest version of the NCBI blastall program can be downloaded from:
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST

=back

=head2 Required Perl Modules

=over

=item * Bio::SearchIO

This module is part of bioperl and is required to parse BLAST
output in a format that is tiled across all HSPs.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the RepMiner
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

=over

=item * Limited to m8 or m9 BLAST format

This script is designed to be a lean and fast parser of the 
similarity information from BLAST. It is therefore limited
to using the simple m8 or m9 BLAST alignment format.

=back

=head1 SEE ALSO

The cnv_blast2sim.pl program is part of the repminer package of 
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

STARTED: 02/24/2008

UPDATED: 12/06/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# CHANGELOG                                                 |
#-----------------------------------------------------------+
#
# 03/09/2008
# - Added Getopt long 
# - Added infile and outfile as reuired command line options
# - Added option to select BLAST parsing options
#    1 --> Bitscore of Best HSP from -m8 BLAST output
#    2 --> Bitscore of Tiled HSP from -m8 BLAST output 
#          This option uses the BioPerl BLAST parsing
#    3 --> Significance of the best scoring HSP
#    4 --> Significance of Tiled HSPs
# 
# 07/15/2008
# - Renamed from extract_blast.pl to cnv_blast2sim.pl
#   this name reflects the fact this this converts native
#   m8 BLAST output to a simple three column similarity
#   file.
# - Added the option to print output to STDOUT when the
#   the infile argument is not given
# - Added the option for input from STDIN
#
# 12/06/2008
# - Updated help file information
