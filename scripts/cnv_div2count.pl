#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_div2count - Convert diversity file to count           |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 06/09/2009                                       |
# UPDATED: 06/09/2009                                       |
#                                                           |
# DESCRIPTION:                                              | 
#   Convert a basic diversity file to the count of the      |
#   richness of the values in that bin.                     |
#                                                           |
# LICENSE:                                                  |
#  GNU Lesser Public License                                |
#  http://www.gnu.org/licenses/lgpl.html                    |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use DBI;
use Getopt::Long;

my $infile;                      # Full path to the input file to parse
my $outfile;                     # Full path to the diversity output file

# BOOLEANS
my $verbose = 0;
my $show_help = 0;             # Display help
my $show_man = 0;              # Show the man page via perldoc
my $show_usage = 0;            # Show the basic usage for the program
my $show_version = 0;          # Show the program version
my $quiet = 0;                 # Run the program in quiet mode
                               # will not prompt for command line options
#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
                    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # BOOLEANS
		    "q|quiet"    => \$quiet,
		    "verbose"    => \$verbose,
		    "version"    => \$show_version,
		    "man"        => \$show_man,
		    "usage"      => \$show_usage,
		    "h|help"     => \$show_help,);


#-----------------------------+
# OPEN THE INPUT FILE         |
#-----------------------------+
# Default to STDIN otherwise
if ($infile) {
    open (DIVIN , "<$infile") ||
	die "Can not open gff file for intput\n";
}
else {
    open (DIVIN , "<&STDIN") ||
	die "Can not accept input from STDIN\n";
    print STDERR "Expecting input from STDIN\n";
}

#-----------------------------+
# OPEN THE OUTPUT FILE        |
#-----------------------------+
if ($outfile) {
    open (DIVOUT, ">$outfile") ||
	die "Can not open output file:\n";
}
else {
    open (DIVOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output\n";
}


#-----------------------------+
# LOAD GFF FILE TO DATABASE   |
#-----------------------------+
my $line_num = 0;
while (<DIVIN>) {

    # Increment Line Number
    $line_num++;

    # Skip head information
    next if m/^\#/;

    chomp;
    print STDERR "Processing line $line_num\n" if $verbose;

    my @div_vals = split;
    
    my $div_rich = 0;             # Richness

    for my $ind_val (@div_vals) {
	# Only process numbers
	if ($ind_val =~ m/^[\d]*$/ ) {
	    if ( int($ind_val) > 0 ) {
		$div_rich++;
	    }
	}
    }

    print DIVOUT $div_vals[0]."\t$div_rich\n";

}


#-----------------------------+
# CLOSE FILE HANDLES          |
#-----------------------------+
close (DIVIN);
close (DIVOUT);

exit;




#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
__END__


#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 06/16/2009
# - Program started
