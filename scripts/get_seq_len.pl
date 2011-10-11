#!/usr/bin/perl -w
# Get length of sequences
#-----------------------------+
# 06/11/2009
# James C. Estill
# James.Estill_at_gmail.com

# INCLUDES
use strict;
use Getopt::Long;
use Bio::SeqIO;

# Vars set at command line
my $infile;
my $outfile;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_version = 0;
my $show_man = 0;
my $show_usage = 0;
my $show_help = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
                    "i|infile=s" => \$infile,
		    "o|outfile=s" => \$outfile,
		    # ALTERNATIVE TO --dsn 
		    # ADDITIONAL OPTIONS
		    # BOOLEANS
		    "q|quiet"    => \$quiet,
		    "verbose"    => \$verbose,
		    "version"    => \$show_version,
		    "man"        => \$show_man,
		    "usage"      => \$show_usage,
		    "h|help"     => \$show_help,);

#-----------------------------+
# OPEN THE OUTPUT FILE        |
#-----------------------------+
if ($outfile) {
    open ( RMOUT ,">$outfile") ||
	die "Can not open file for output:\n$outfile\n";
} 
else {
    open ( RMOUT ,">&STDOUT") ||
	die "Can not STDOUT for output:\n";
}


#-----------------------------+
# OPEN INPUT FILE             |
#-----------------------------+
# If an input file is specified use the infile otherwise
# expect input from STDIN
my $seq_in;
if ($infile) {
    $seq_in = Bio::SeqIO->new ( '-format' => 'fasta',
				'-file' => "<$infile" )
	|| die "Can not open infile:\n$infile\n";
}
else {
    $seq_in = Bio::SeqIO->new ( '-format' => 'fasta',
				'-fh' => \*STDIN)
	|| die "Can not open input from STIN\n";
}

#-----------------------------+
# GET SEQUENCE LENGTHS        |
#-----------------------------+
while ( my $seq = $seq_in->next_seq()) {

    print RMOUT $seq->primary_id()." - ".$seq->length()."\n";
    
}

close RMOUT;

exit;
