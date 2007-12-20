# Test perl script
# It will be used to extract seq features from the annotate Arabidopsis genome available from TIGR
# THE URL HERE
#
# Author: James C. Estill
# Contact: jestill@uga.edu
# Started: 7/2/04
# Last Visited: 7/8/04

#
use Bio::SeqIO;
# get command-line arguments, or die with a usage statement
my $usage = "x2y.pl infile infileformat outfile outfileformat\n";
my $infile = shift or die $usage;
my $infileformat = shift or die $usage;
my $outfile = shift or die $usage;
my $outfileformat = shift or die $usage;
	  
# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
                              '-format' => $outfileformat);

 # write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq)
{
      $seq_out->write_seq($inseq);
}

exit;
