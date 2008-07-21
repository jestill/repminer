#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_mclerr2tab.pl - Convert MCL STDERR to tab delim text  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/20/2008                                       |
# UPDATED: 07/20/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Parse the mcl stderr files in a directory and return     |
#  the results returned by mcl in a tab delimited           |
#  text file.                                               |
#                                                           |
# EXAMPLES:                                                 |
# cnv_cpm2tab.pl -i infile.cpm -o outfile.cpm.txt           |
#                                                           |
#-----------------------------------------------------------+
# TODO:

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # 
use Getopt::Long;              # Get options from command line
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path


#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev:$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outfile;

# BOOLEANS
my $verbose = 0;
my $quiet = 0;
my $show_usage = 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0; 
my $debug = 0;

# Index Vals
my $pre_cat_id ="0";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
#$indir = shift;
#$outfile = shift;

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|indir=s"    => \$indir,
                    "o|outfile=s"   => \$outfile,
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
if (!$indir) {
    print "An input directory is required.";
    exit;
}


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

#-----------------------------------------------------------+
# MAIN BODY                                                 |
#-----------------------------------------------------------+
# Print to outfile if one specified, otherwise to STDOUT
if ($outfile) {
    open (TABOUT, ">$outfile") ||
	die "Can not open output file\n$outfile\n";
} 
else {
    open (TABOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output\n";
}


#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @err_files = grep /\.mcl\.err$/, readdir DIR ;
closedir( DIR );
my $count_files = @err_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print "\a";
    print "\nERROR: No mcl STDERR files were found in the input directory\n".
	"$indir\n".
	"These files must have the extension *mcl.err .\n\n";
    exit;
}


#-----------------------------------------------------------+
# PARSE CLM TO TAB DELIMITED OUTPUT                         |
#-----------------------------------------------------------+
for my $ind_file (@err_files)
{

    my $file_path = $indir.$ind_file;
    #print "\a";


    # VARAIBLES TO EXTRACT
    my $num_nodes = "";          # The number of nodes in the matrix
    my $mcl_i_val = "";          # The inflation value used in mcl
    my $mcl_s_val = "";          # The scheme value used
    my $mcl_total_time = 0;     # Total time to finish iterations
    my $mcl_num_clusters = "";   # Number of mcl clusters found
    my $mcl_prune_synopsis = ""; # Synopsis score
    my $mcl_prune_mark_1 = "";   #
    my $mcl_prune_mark_2 = "";   #
    my $mcl_prune_mark_3 = "";   #
    my $mcl_num_iter = "";       # Total number of iterations
    my $mcl_out_file = "";
    #my $iter_time = "0";

    # OPEN THE MCL STDERR FILE
    open (MCLERRIN, "<$file_path") ||
	die "Can not open the infile\n$file_path\n";
    
    # PARSE THE FILE
    while (<MCLERRIN>) {
	chomp;
	
	#print "$_\n";

	#-----------------------------+
	# GET ITERAND INFORMATION     |
	#-----------------------------+
	if (m/( *)(\d*)( *)\.{19}( *)(.*)/) {

	    # The following for debugging
	    if ($debug) {
		print STDERR "$_\n";
		print STDERR "One---".$1."-\n" if ($1);	    
		print STDERR "Two---".$2."-\n" if ($2);	    
		print STDERR "Three-".$3."-\n" if ($3);
		print STDERR "Four--".$4."-\n" if ($4);
		print STDERR "Five--".$5."-\n" if ($5);
		print STDERR "Six---".$6."-\n" if ($6);
		print STDERR "Seven-".$7."-\n" if ($7);
	    }

	    if ($5) {
		my @iter_parts = split(/ /, $5);

		my $num_iter_parts = @iter_parts;

		if ($debug) {
		    print STDERR "Parts $num_iter_parts \n";
		    print STDERR "\t nip0: ".$iter_parts[0]."\n" 
			if ($iter_parts[0]);
		    print STDERR "\t nip1: ".$iter_parts[1]."\n" 
			if ($iter_parts[1]);
		    print STDERR "\t nip2: ".$iter_parts[2]."\n" 
			if ($iter_parts[2]);
		    print STDERR "\t nip3: ".$iter_parts[3]."\n" 
			if ($iter_parts[3]);
		}
		
		if ($num_iter_parts == 3) {
		    my $iter_time = $iter_parts[1];
		    $mcl_num_iter = $2;
		    $mcl_total_time = $iter_time + $mcl_total_time;
		    #my @slash_parts = split(/\//, $iter_parts[2]);
		}
		elsif ($num_iter_parts == 4) {
		    my $iter_time = $iter_parts[2];
		    #print "$_\n";
		    $mcl_num_iter = $2;
		    $mcl_total_time = $iter_time + $mcl_total_time;
		    #my @slash_parts = split(/\//, $iter_parts[3]);
		}


	    }
	    
	}

	#-----------------------------+
	# GET NUMBER OF CLUSTERS      |
	#-----------------------------+
	if (m/\[mcl\] (.*) clusters found/) {
	    $mcl_num_clusters=$1;
	}


	#-----------------------------+
	# GET OUTPUT FILE NAME        |
	#-----------------------------+
	if (m/output is in (.*)/) {
	    #print "$_\n";
	    #print $1."\n";
	    $mcl_out_file = $1;
	    # Get s and I value from file name if possible
	    if ($mcl_out_file =~ m/.*\/(.*)i\-(.*)_s\-(.*).mcl/) {
		$mcl_i_val = $2;
		$mcl_s_val = $3;
	    }
	}
	
	#-----------------------------+
	# JURY PRUNING SYNOPSIS       |
	#-----------------------------+
	if (m/jury pruning synopsis: <(.*)>/) {
	    $mcl_prune_synopsis = $1;
	    #print "";
	}

	#-----------------------------+
	# JURY PRUNING MARKS          |
	#-----------------------------+
	if (m/jury pruning marks(.*)<(.*)>/) {
	    my @prune_marks = split(/,/, $2);
	    $mcl_prune_mark_1 = $prune_marks[0];
	    $mcl_prune_mark_2 = $prune_marks[1];
	    $mcl_prune_mark_3 = $prune_marks[2];
	}

	#-----------------------------+
	# GET NUMBER OF NODES         |
	#-----------------------------+
	if (m/matrix with (.*) entries/) {
	    $num_nodes = $1;
	}

    }
    close MCLERRIN;
    
    #-----------------------------+
    # PRINT EXTRACTED DATA        |
    #-----------------------------+
    if ($debug) {
	print STDERR "-----------------------------\n";
	print STDERR " File: ".$mcl_out_file."\n";
	print STDERR "I-Val: ".$mcl_i_val."\n";
	print STDERR "S-Val: ".$mcl_s_val."\n";
	print STDERR "Nodes: ".$num_nodes."\n";
	print STDERR " Iter: ".$mcl_num_iter."\n";
	print STDERR "Clust: ".$mcl_num_clusters."\n";
	print STDERR " Time: ".$mcl_total_time."\n";
	print STDERR "Synop: ".$mcl_prune_synopsis."\n";
	print STDERR "Mark1: ".$mcl_prune_mark_1."\n";
	print STDERR "Mark2: ".$mcl_prune_mark_2."\n";
	print STDERR "Mark3: ".$mcl_prune_mark_3."\n";
	print STDERR "-----------------------------\n";
	print STDERR "\n";
    }
    
    print TABOUT $mcl_out_file."\t";
    print TABOUT $mcl_i_val."\t";
    print TABOUT $mcl_s_val."\t";
    print TABOUT $num_nodes."\t";
    print TABOUT $mcl_num_iter."\t";
    print TABOUT $mcl_num_clusters."\t";
    print TABOUT $mcl_total_time."\t";
    print TABOUT $mcl_prune_synopsis."\t";
    print TABOUT $mcl_prune_mark_1."\t";
    print TABOUT $mcl_prune_mark_2."\t";
    print TABOUT $mcl_prune_mark_3."\t";
    print TABOUT "\n";


}

close (TABOUT);

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

cnv_mclerr2tab.pl - Convert mcl STDERR output to tab delim text

=head1 VERSION

This documentation refers to cnv_mclerr2tab.pl version $Rev:$

=head1 SYNOPSIS

=head2 Usage

    cnv_mclerr2tab.pl -i indir/ -o outfile.txt

=head2 Commonly Used Arguments

    -i, --indir     # Path to the dir with the files to parse
    -o, --outfile   # Path to the output text file

=head1 DESCRIPTION

Parse the mcl stderr files in a directory and return
the results returned by mcl in a tab delimited
text file.

=head1 COMMONLY USED ARGUMENTS

=over 2

=item -i,--indir

Directory containing the mcl stderr output files to parse.

=item -o,--outfile

Path of the output file. If not outfile arugment is used, the program
will send data to standard output.

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

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: Error message

Error messages will go here.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This script does not make use of external configuration files
or any options in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * MCL

This program will parse the output of MCL Markov Clustering program
available at 
( http://micans.org/mcl/ ).

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=192812

=back

=head2 Limitations

=over

=item * Known to work with MCL 1.007-grumpy-gryphon, 08-157.
Earlier or later versions may not be correctly parsed. 

=back

=head1 SEE ALSO

This program is part of the RepMiner suite of programs.
For more information see
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

STARTED: 07/20/2008

UPDATED: 07/20/2008

VERSION: $Rev:$

=cut

#-----------------------------------------------------------+
# CHANGELOG                                                 |
#-----------------------------------------------------------+
#
# 07/20/2008
# - Program started with the design to to parse the
#   mcl.err files in a directory.

