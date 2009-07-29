#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_div2gff - Convert diversity file to gff file          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 06/09/2009                                       |
# UPDATED: 07/29/2009                                       |
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
my $out_val = "rich";            # The output value to report, coverage or richness ..
                                 # Richness (rich), or coverage (cover)
my $feature_name = "rm_result";   
my $feature_source = "RepeatMasker";
my $feature_strand = "+";
my $feature_frame = ".";

# BOOLEANS
my $do_chrom_int = 0;          # Convert chromosome id to an integer
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
		    "v|value=s"   => \$out_val,
		    "n|name=s"    => \$feature_name,
		    "s|source=s"  => \$feature_source,
		    "strand=s"    => \$feature_strand,
		    "frame=s"     => \$feature_frame,
		    # BOOLEANS
		    "chrom-int"  => \$do_chrom_int,
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
    open (GFFOUT, ">$outfile") ||
	die "Can not open output file:\n";
}
else {
    open (GFFOUT, ">&STDOUT") ||
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
    my $div_score = 0;             # Richness


    # Add choice for coverage or richness here
    if ($out_val = "rich") {
	for my $ind_val (@div_vals) {
	    # Only process numbers
	    if ($ind_val =~ m/^[\d]*$/ ) {
		if ( int($ind_val) > 0 ) {
		    $div_score++;
		}
	    }
	}
    }
    elsif ($out_val = "cover") {

    }

    my ($chrom,$span) = split(/\:/, $div_vals[0]);
    my ($start,$end) = split(/\-/,$span);

    # Convert chromosome id to integer if not already an integer to facilitate drawing    
    if ($do_chrom_int) {
	unless ($chrom =~ m/^[\d]*$/ ) {
	    # default output from the agp_analysis program is like the following
	    if ($chrom =~ m/^chr\_(\d*)$/) {
		$chrom = $1;
	    }
	}
    }

    my $feature_attribute = $feature_name."_".$start."_".$end;

    my $gff_out = $chrom."\t".    # chromosome
	$feature_source."\t".     # source of the annotation
	$feature_name."\t".       # feature type
	$start."\t".              # start
	$end."\t".                # end
	$div_score."\t".          # diversity score
	$feature_strand."\t".     # feature strand (+,-, or .)
	$feature_frame."\t".      # feature frame
	$feature_attribute."\n";  # label for the feature

    print GFFOUT "$gff_out";

}


#-----------------------------+
# CLOSE FILE HANDLES          |
#-----------------------------+
close (DIVIN);
close (GFFOUT);

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
# 07/29/2009
# - Modified cnv_div2gff.pl from the program cnv_div2count.pl

