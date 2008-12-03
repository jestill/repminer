#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fasta2tenestdb.pl - Convert fasta to TENest database  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/21/2008                                       |
# UPDATED: 11/21/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a multiple record fasta file, split the file into  |
#  the directory and file structure expected by the         |
#  TENest program.                                          |
#                                                           |
# USAGE:                                                    |
#  cnv_fasta2tenestdb.pl -i infile.fasta -o outdir          |
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
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path
use Bio::SeqIO;                # Use Bio::SeqIO for the reading and 
                               # writing of seq files in different formats

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                      # Input file of seqs
my $outdir;                      # Base output directory
my $ltrfile;                     # The fasta file containing the LTRs

# BOOLEANS
my $test = 0;
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $do_revcom = 0;                # Do the reverse complement
                                  # This assumes that input files are in +
my $do_list = 0;                  # Genrate the ltr_list file

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
		    "l|ltrfile=s" => \$ltrfile,
		    "r|do-revcom" => \$do_revcom,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "do-list"     => \$do_list,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$test,
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
    print STDERR "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
        die "ERROR: Could not create the output directory:\n$outdir\n";
}


#-----------------------------+
# CREATE THE TENEST DIRS      |
#-----------------------------+

my $condir = $outdir."CON-1/";
unless (-e $condir) {
    print STDERR "Creating CON-1 dir ...\n" if $verbose;
    mkdir $condir ||
	die "ERROR: Count not create the CON-1 dir at:\n$condir\n";
}

my $tedir = $outdir."TEs/";
unless (-e $tedir) {
    print STDERR "Creating the TEs dir ...\n" if $verbose;
    mkdir $tedir ||
	die "ERROR: Can not create the TEs dir at:\n$tedir\n";
}

my $ltrdir = $outdir."LTRs/";
unless (-e $ltrdir) {
    print STDERR "Creating the LTRs dir ...\n" if $verbose;
    mkdir $ltrdir ||
	die "ERROR: Can not create the LTRs dir at:\n$ltrdir\n";
}


#-----------------------------+
# OPEN THE LTR LIST FILE      |
#-----------------------------+
if ($do_list) {
    my $ltr_list_out = $outdir."ltr_list";
    open (LTRLIST, ">$ltr_list_out") ||
	die "Can not open LTR list file at:\n$ltr_list_out\n";
}

#-----------------------------------------------------------+
# GENERATE OUTPUT FILES FOR THE INTPUT FILE                 |
#-----------------------------------------------------------+

# create one SeqIO object to read in,and another to write out
my $infileformat = "fasta";

my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);


while (my $inseq = $seq_in->next_seq) {


    print STDERR "Processing: ".$inseq->primary_id()."\n" if $verbose;


    #-----------------------------+
    # PRINT FILE NAME TO LTR LIST |
    #-----------------------------+
    if ($do_list) {
	print LTRLIST $inseq->primary_id()." 0\n";
    }

    #-----------------------------+
    # PRINT FILE TO CON-1 DIR     |
    #-----------------------------+
    my $con_path = $condir.$inseq->primary_id();
    
    # The consenus sequence object
    # The TENest program appears to refer to sequences in the 
    # CON-1 directory as consensus sequences
    my $con = Bio::SeqIO->new('-file' => ">$con_path",
			      '-format' => $infileformat);

    $con->write_seq($inseq);

    #-----------------------------+
    # DO REVCOM IF REQUESTED      |
    #-----------------------------+
    if ($do_revcom) {

	#-----------------------------+
	# POSITIVE ORIENTATION OBJECT |
	#-----------------------------+
	my $rev_p_path = $tedir.$inseq->primary_id();

	my $rev_p = Bio::SeqIO->new('-file' => ">$rev_p_path",
				    '-format' => $infileformat);
	
	$rev_p->write_seq($inseq);

	#-----------------------------+
	# REVERSE COMPLEMENT OBJECT   |
	#-----------------------------+
	my $rev_r_path = $tedir.$inseq->primary_id()."_R";

	my $rev_r = Bio::SeqIO->new('-file' => ">$rev_r_path",
				    '-format' => $infileformat);

	my $rev_com = $inseq->revcom();

	my $rev_com_id = $inseq->primary_id()."_R";
	$rev_com->display_id( $rev_com_id );

	my $rev_com_id = $rev_com->display_id();
	
	print STDERR "\t$rev_com_id\n" if $verbose;

	$rev_r->write_seq($rev_com);
	

    } # End of do the revcom

}

# CLOSE THE LTR LIST FILE HANDLE
if ($do_list) {
    close LTRLIST;
}

#-----------------------------------------------------------+
# GENERATE OUTPUT FILES FOR THE LTRS                        |
#-----------------------------------------------------------+
# This requires that the LTR file was passed as --ltrfile

if ($ltrfile) {
    
    my $ltr_in = Bio::SeqIO->new('-file' => "<$ltrfile",
				 '-format' => $infileformat);
    
    #-----------------------------+
    # GENERATE AN OUTPUT FILE FOR |
    # EACH FILE IN THE LTR FILE   |
    #-----------------------------+
    while (my $inseq = $ltr_in->next_seq) {
	
	my $ltr_path = $ltrdir.$inseq->primary_id()."-0.fasta";
	
	my $ltr_out = Bio::SeqIO->new('-file' => ">$ltr_path",
				      '-format' => $infileformat);
	
	my $ltr_seq_id = $inseq->primary_id()."-0.fasta";
	$inseq->display_id( $ltr_seq_id );
	
	$ltr_out->write_seq($inseq);

    }

}

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

1;
__END__

=head1 NAME

cnv_fasta2tenestdb.pl - Convert fasta to TENest database

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    cnv_fasta2tenestdb.pl -i infile.fasta -o outdir

    --infile        # Path to the input file
    --outdir        # Base directory to hold output files 

=head1 DESCRIPTION

Given a multiple record fasta file, split the file into
the directory and file structure expected by the
TENest program.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outdir

Base output directory where the files will be stored.

=back

=head1 Additional Options

=over 2

=item -l,--ltrfile

The multiple record fasta file that contains the LTRs of the LTR
retrotransposons.

=item -r,--do-revcom

This option will create a reverse complement of the sequences used 
as input. This data is used by TENest to generate the direction
of the insertion relative to the TENest query sequence. The output
of the forward and reverse strands will be stored as separate files
in the TEs directory. This assumes that all of the input files are
in the positive orientation.

=item --do-list

Create the ltr_list file in the output directory.

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

This program does not make use of configuration files or any 
variables in the user environment.

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

STARTED: 11/21/2008

UPDATED: 11/21/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/21/2008
# - Program started. The purpose is take fasta files from
#   a database and genreate a database of files for
#   the TENest program.

