#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2na.pl - Convert blast output to cytoscape NA    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/06/2008                                       |
# UPDATED: 06/06/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert m8 blast output to the cytoscpae na file format  |
#                                                           |
# USAGE:                                                    |
#  cnv_blast2na.pl -i blastfile.bln -o node_attribute.na    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# Only the best blast hit is used to indenitfy the class

package REPMINER;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $na_header;                 # header to use in the NA file

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

my $do_best = 0;               # just return the best blast hit

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "best"        => \$do_best,
		    "header=s"    => \$na_header,
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

#-----------------------------------------------------------+
# FILE I/O HANDLES                                          |
#-----------------------------------------------------------+
open (INFILE, "<$infile") ||
    die "Can not open the blast file for parsing\n";

open (NAOUT, ">$outfile") ||
    die "Can not open the Node attribute file for writing\n";

if ($na_header) {
    print NAOUT "$na_header\n";
}
else {
    print NAOUT "BLAST_".$infile."\n";
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $prev_qry = "";

while (<INFILE>) {
    
    # Ignore comment lines
    next if m/^\#.*/;
    chomp;
    
    my ($qry_id, $sub_id, $pid, $len, 
	$mis_match, $gap_open, 
	$qry_start, $qry_end, $sub_start, $sub_end,
	$e_val, $bit_score) = split(/\t/);
    
    # First check if this is a duplicate hit
    #print STDERR $prev_qry." - ".$qry_id."\n";
    next if $prev_qry =~ $qry_id;
    
    #print STDERR "\tprinting\n";
    # Set the prev_qry id for the next line
    $prev_qry = $qry_id;

    # The following is optimized for the way that I identify 
    # sequences in repminer using pipe delimited numbers in the
    # header of the fasta files
    my @split_qry = split (/\|/, $qry_id);
    
    my $qry_num_id = $split_qry[0];
    
    print STDERR "Processing $qry_num_id\t$sub_id\n" if $verbose;
    print NAOUT $qry_num_id."=".$sub_id."\n";

#    # If just doing the best (ie first blast hit)
#    # chieck to see if we have already printed
#    if ($do_best) {
#	unless ($prev_num_id == $qry_num_id) {
#	    print STDERR "Processing $qry_num_id\t$sub_id\n" if $verbose;
#	    print NAOUT $qry_num_id."=".$sub_id."\n";
#	}
#    }
#    else {
#	print STDERR "Processing $qry_num_id\t$sub_id\n" if $verbose;
#	print NAOUT $qry_num_id."=".$sub_id."\n";
#    }
#
#    $prev_num_id = $qry_num_id;

}

close (INFILE);
close (NAOUT);

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


=head1 NAME

Name.pl - Short program description. 

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

STARTED: 06/06/2008

UPDATED: 06/06/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/06/2008
# - Main body of program written, a simple parser to convert
#   the m8 blastn output to a node attribute file. This 
#   selects the first hit to use to identify the node.
# - Useful to look at blast against existing database of 
#   transposable elements
