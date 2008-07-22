#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# blast_getseq.pl - Get best hit protein seq from blast out |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/30/2008                                       |
# UPDATED: 03/30/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Extract the sequence of the best hit from a blast report |
#                                                           |
# USAGE:                                                    |
#  blast_getseq.pl -i blast_report.blo -o bh_seqs.fasta     |
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
use Bio::SearchIO;             # Parse BLAST output

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

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
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


# Check required optiosn
if ( (!$infile) || (!$outfile) ) {
    print "Error: Not all options present at command line\n";
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

blast_getseq ($infile, $outfile, "1");

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


sub blast_getseq {
    

    my ($blastin, $seqout, $max_num_hits ) = @_;
    # $blastin --> path to the blast output file
    # $seqout ---> path to the fasta output file
    # $num_hits -> number of hits to return seqs for  
    
    #-----------------------------+
    # FILE IO                     |
    #-----------------------------+
    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
					   '-file'   => $blastin )
	|| die "Could not open BLAST input file:\n$blastin.\n";

    open (SEQOUT, ">$seqout") ||
	die "Can not open sequence outfile $seqout\n"; 

    #-----------------------------+
    # PROCESS BLAST REPPORT       |
    #-----------------------------+
    while (my $blast_result = $blast_report->next_result())
    {
	# Initialize the hit counter 
	my $hit_count = 0;
	
	while (my $blast_hit = $blast_result->next_hit())
	{
	    
	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {
		if ($hit_count < $max_num_hits) {
		    print STDERR ">".$blast_result->query_name."\n" if $verbose;
		    print STDERR $blast_hsp->query_string."\n" if $verbose;
		    
		    print SEQOUT ">".$blast_result->query_name."\n";
		    print SEQOUT $blast_hsp->query_string."\n"; 

		    # Should not be here but prevents double vals for
		    # split hsps
		    $hit_count++;

		}
		
	    }

	    $hit_count++;

	} # End of while next hit

    } # End of while blast_report next result
    
    close (SEQOUT);

    # Return sub function worked
    return 1;

}

=head1 NAME

blast_getseq.pl - Get best hit protein seq from blast out

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    blast_getseq.pl -i InFile -o OutFile
    
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
