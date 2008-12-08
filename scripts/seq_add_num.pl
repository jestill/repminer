#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq_add_num.pl - Add integer prefix to sequence headers   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/12/2006                                       |
# UPDATED: 12/07/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Add a number prefix to the id line from a fasta file     |
#  to allow for easy parsing of fasta sequences. The process|
#  also provides each sequence with a unique incremental    |
#  id that can be used for reference. Although designed     |
#  to work with fasta format files, any bioperl valid       |
#  sequence file format may be used for input or export.    |
#                                                           |
# USAGE:                                                    |
#  seq_add_num.pl -i infilepath -o outfilepath              |
#                                                           |
# REQUIREMENTS:                                             |
#  -bioperl                                                 |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Bio::SeqIO;                # Allows for treatment of seqs as objects
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
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $infile;
my $outfile;

my $seq_in_format = "fasta";   # The input sequence format
my $seq_out_format = "fasta";  # The output sequence format
my $seq_num = 1;               # Var to keep track of nummber of seqs 
my $seq_unique_id;             # Unique ID attributed to the seq read record
my $new_id;

# BOOLEANS
my $verbose = 0;
my $line_num = 0;
my $quiet = 0;
my $show_usage = 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0; 
my $blast_opt = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+

my $ok = GetOptions(# REQUIRED OPTIONS
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
                    # ADDITIONAL OPTIONS
		    "s|start-num"   => \$seq_num,
		    "in-format"     => \$seq_in_format,
		    "out-format"    => \$seq_out_format,
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
# FASTA FILE IO               |
#-----------------------------+
my $inseq = Bio::SeqIO->new(-file   => "<$infile",
			    -format => $seq_in_format ) ||
    die "ERROR Can not open infile:\n $infile\n";

my $outAll = Bio::SeqIO->new(-file   => ">$outfile",
			     -format => $seq_out_format ) ||
    die "ERROR Can not open outfile:\n $outfile\n";

#-----------------------------+
# PROCESS SEQUENCE FILE       |
#-----------------------------+

while (my $seq = $inseq->next_seq)  {
    # This process may not work with all FASTA files, some
    # of the information will be left out of the new FASTA file
    # that is produced. It would be good to just capture the entire
    # FASTA header and append the incremented number.
    $seq_unique_id = $seq->primary_id;
    $new_id = $seq_num."|".$seq_unique_id ;

    # Restet the seqid to the new id
    $seq->primary_id( $new_id );
    $seq->display_id( $new_id );

    #-----------------------------+
    # PRINT THE NEW ID            |
    #-----------------------------+
    print STDERR "\tNEW ID: ".$new_id."\n" if $verbose;

    #-----------------------------+
    # WRITE SEQUENCE RECORDS OUT  |
    # TO A NEW FILE               |
    #-----------------------------+
    $outAll->write_seq($seq);
	
    # Increment the sequence number
    $seq_num++;


} # END OF THE FOR EVERY SEQUENCE RECORD 

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

seq_add_num.pl - Add integer prefix to sequence headers

=head1 VERSION

This documentation refers to seq_add_num.pl version $Rev$

=head1 SYNOPSIS

=head2 Usage

    seq_add_num.pl -i infile.fasta -o outfile_num.fasta

=head2 Required Arguments

    -i, --infile    # Path to the input file to add prefix integer to
    -o, --outfile   # Path to the output file

=head1 DESCRIPTION

Add a number prefix to the id line from a fasta file
to allow for easy parsing of fasta sequences. The process
also provides each sequence with a unique incremental
id that can be used for reference. Although designed
to work with fasta format files, any bioperl valid
sequence file format may be used for input or export.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the sequence file to add the prefix to the id.
If no input file argument is used, the program will attempt to read from STDIN.

=item -o,--outfile

Path of the output file. If not outfile arugment is used, the program
will send data to standard output.

=back

=head1 OPTIONS

=over 2

=item --in-format

The file format of the sequence to prepend a sequence ID to. Be default the
program will expect sequenced to be in the fasta format, however any
valid sequence format accpeted by bioperl may be used.

=item --out-format

The file format of the output sequence file. Be default, fasta format files
will be produced, but any file format acceptable by bioperl may be used.

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

=back

=head1 EXAMPLES

=head2 Typical Use

To use the program to parse a fasta format file to a fasta format file
where the headers have unique integers as a prefix:

 seq_add_num.pl -i infile.fasta -o out_num.fasta

This will add an integer to the header files such that the following headers:

 >seq_one
 ATGT...
 >seq_two
 ATGT..

Would be modified to:

 >1|seq_one
 ATGT...
 >2|seq_two
 ATGT...

This will assure that in future analysis, all sequence records will have
a single unique integer assigned to that sequence.

=head2 Specify Seq Number Starting Point

It is also possible to specify where the integer should begin incrementing
from. This is useful when you want to combine sequence from two 
separate files, but want the combined sequence records to
have unique identifiers. For example the command:

 seq_add_num.pl -i infile.fasta -o out_num.fasta --start-num 76

would modify the seuqence headers to:

 >76|seq_one
 ATGT...
 >77|seq_two
 ATGT...

=head2 Alternatives to FASTA File Format

It is also possible to use sequence records that are not in the
FASTA file format. The following command will convert a genbank format
file to a fasta format file with integer prefixes starting at one.

 seq_add_num.pl -i infile.gb --in-format genbank -o out_num.fasta

The program assumes that you want the output in FASTA file format, so
it is not necessary to specify FASTA format for output.

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

=head2 Required Perl Modules

=over

=item * Bio::SeqIO

This module is part of bioperl and is required to read and write
sequence records in multiple formats.

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

=item * Long sequence header information may be lost

This program will only use the portion of the fasta header that is
considered to be the primary id. The unique integer will be prepened
to the primary id, and this will be used as the new id.

=back

=head1 SEE ALSO

The seq_add_num.pl program is part of the repminer package of 
repeat element annotation programs.
See the DAWG-PAWS web page 
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

STARTED: 07/12/2006

UPDATED: 12/07/2008

VERSION: $Rev$

=cut


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/12/2006
# - Started program
#
# 12/18/2006
# - Changed Name to FastaAddNum.pl
# - Added input variable flags
#
# 12/07/2008
# - Adding POD documentation
# - Adding help subfunction
