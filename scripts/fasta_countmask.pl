#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_countmask.pl - count fasta file masked characters   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/11/2008                                       |
# UPDATED: 08/14/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Counts the masked characters in a fasta file. Reports    |
#  the number of masked characters, the total number        |
#  of characters and the final ratio. Assummes lower        |
#  case masking.                                            |
#                                                           |
# USAGE:                                                    |
#  fasta_countmask.pl -infile.fasta -o outreport.txt        |
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
my $do_full = 0;               # Do analysis for each sequence

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "full"        => \$do_full,
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
    print "\nfasta_countmask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

# INFILE
if ($infile) {
    open (FASTIN, "<$infile");
}
else {
    open (FASTIN, "<&STDIN");
}

# OUTFILE
if ($outfile) {
    open (SUMOUT, ">$outfile");
} 
else {
    open (SUMOUT, ">&STDOUT");
}


my $num_mask_char = 0;       # Global number of masked characters
my $num_char = 0;      # Global number of characters
my $line_num;

# Individual sequence values
my $seq_id;
my $seq_mask_char=0;
my $seq_char=0;
my $seq_ratio=0;

while (<FASTIN>) {
    
    chomp;
	
    $line_num++;
    
    if  (m/^>/) {

	if ($do_full) {
	    if ($seq_char > 0) {
		$seq_ratio = $seq_mask_char / $seq_char;
		print SUMOUT $seq_id."\t";
		print SUMOUT $seq_mask_char."\t";
		print SUMOUT $seq_char."\t";
		print SUMOUT $seq_ratio."\n";
	    }
	    
	    # Reset counters
	    $seq_mask_char=0;
	    $seq_char=0;
	    $seq_ratio=0;
	}	
	
	# Reset the seq id
	$seq_id = substr ($_, 1);
    }


    # Skip the rest for FASTA headers
    next if m/^>/; 

    # Use the translate regular expression without 
    # replacement to count the lowercase characters
    my $num_line_mask_char = $_ =~ tr/[a-z]//;
    my $num_line_char = length $_;

    # Add to the running count for the sequence
    if ($do_full) {
	$seq_mask_char = $seq_mask_char + $num_line_mask_char;
	$seq_char = $seq_char + $num_line_char;
    }

    # Add to the running global count
    $num_mask_char = $num_mask_char + $num_line_mask_char;
    $num_char = $num_char + $num_line_char;

    # Show what is being parsed for debug purposes
    print STDERR $_."\n" if $verbose;
    print STDERR "LWC CHAR: $num_line_mask_char\n" if $verbose;
    print STDERR "TOT CHAR: $num_line_char\n" if $verbose;
    

}

my $percent_mask_char = $num_mask_char / $num_char;

if ($do_full) {
    # Print the last seq record
    $seq_ratio = $seq_mask_char / $seq_char;
    print SUMOUT $seq_id."\t";
    print SUMOUT $seq_mask_char."\t";
    print SUMOUT $seq_char."\t";
    print SUMOUT $seq_ratio."\n";
}

# The overall summary file
# This goes to STDOUT
print STDOUT "Masked:\t$num_mask_char\n";
print STDOUT "Total:\t$num_char\n";
print STDOUT "Ratio:\t$percent_mask_char\n";

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

__END__

# OLD PRINT HELP SUBFUNCTION

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

fasta_countmask.pl -  count fasta file masked characters

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    fasta_countmask.pl -i InFile -o OutFile

    -i,--infile        # Path to the input file
    -o,--outfie        # Path to the output file

=head1 DESCRIPTION

Counts the number of masked characters in a fasta file. Currently
this is limited to lowercase masked characters. The default output
is three lines of text printed to STDOUT that summarizes the number
of bases that were masked, the total number of bases in the input 
file and the ratio of bases that were masked:

  Masked: 7432395
  Total:  7583908
  Ratio:  0.980021777690341

It is also possible to generate summaries for all of the individual
sequence records in the fasta file using the --full option.

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

=item --full

Provide output for the analysis for each individual sequence.
This will generate a four column text file

=over

=item col 1.

Sequence Id.

=item col 2.

The number of bases that were masked for this individual seqeunce record.

=item col 3.

The total number of bases in the individual sequence record.

=item col 4.

The ratio of bases in the individual sequence record that were masked.

=back

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

There are multiple ways to use this program, and since fasta_countmask.pl
accepts pipes, it can easily be made part of a work flow

To generate a summary of the masked characters in the sequence file
my_seq.fasta.masked that has been masked by RepeatMasker:

    fasta_countmask.pl -i my_seq.fasta.masked

To generate a summary file for all masked files in a directory, you
can pipe the text in using the cat program in unix:

    cat *.masked | fasta_countmask.pl 

You can also get information for each seqeunce in your file using the
--full option.

    cat *.masked | fasta_countmask.pl --full

The results of the above examplee will be send the summary data to STDOUT,
in other words the results are printed to your terminal. 
To send the results to an output file use the -o or --outfile option.

    cat *.masked | fasta_countmask.pl --full --outfile summary_file.txt

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

STARTED: 08/11/2008

UPDATED: 08/14/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 08/14/2008
# - Base program started.
#
# 08/14/2008
# - Added the ability to get data for individual sequences
#   by using the --full switch.
