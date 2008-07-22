#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2ea.pl - Convert blast output to cytoscape EA    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/18/2008                                       |
# UPDATED: 06/18/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert m8 blast output to the cytoscpae ea file format. |
#  Returns the tiled bit scores as an edge attribute.       |
#                                                           |
# USAGE:                                                    |
#  cnv_blast2na.pl -i blastfile.bln -o edge_attribute.ea    |
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
use strict;                    # 
use Getopt::Long;              # Get options from command line
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
my $ea_header;                 # header to use in the EA file

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
		    "header=s"    => \$ea_header,
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
my $blast_parse = Bio::SearchIO->new( -file   => $infile,
				      -format => 'blasttable') ||
    die "Can not open blast file for parsing\n$infile\n";


open (EAOUT, ">$outfile") ||
    die "Can not open the Node attribute file for writing\n";

if ($ea_header) {
    print EAOUT "$ea_header\n";
}
else {
    print EAOUT "BLAST_".$infile."\n";
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
while (my $result = $blast_parse->next_result) {
    
    my $qry_id;
    my $qry_name = $result->query_name();
    my @qry_parts = split(/\|/,$qry_name);
    my $num_qry_parts = @qry_parts;
    if ($num_qry_parts > 1) {
	$qry_id = $qry_parts[0];
    }
    else {
	$qry_id = $qry_name;
    }
    
    while (my $hit = $result->next_hit) {
	
	my $hit_id;
	my $hit_name = $hit->name();
	my @hit_parts = split(/\|/,$hit_name);
	my $num_hit_parts = @hit_parts;
	
	if ($num_hit_parts > 1) {
	    $hit_id = $hit_parts[0];
	}
	else {
	    $hit_id = $hit_name;
	}

	print STDERR "Processing $qry_id --> $hit_id\n" if $verbose;

	# PRINT TO EDGE ATTRIBUTE FILE
	# edge attribute as the tiled bitscore
	print EAOUT $qry_id." (bl) ".$hit_id." = ".int($hit->bits())."\n";
	#print STDOUT $qry_id." (bl) ".$hit_id." = ".$hit->bits()."\n";
	
    } # End of while blast next hit
    

} # End of while blast next result

close (EAOUT);

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

STARTED: 06/18/2008

UPDATED: 06/18/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/18/2008
# - Program started from teh cnv_blast2na.pl program
# - This expects the blast file to be in m8 form
