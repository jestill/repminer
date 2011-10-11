#!/usr/bin/perl -w
# -----------------------------------------------------------------------+
#                                                                        |
# CHECK SEQUENCE PARAMETERS                                              |
#                                                                        |
# -----------------------------------------------------------------------+
# seqcheck.pl                                                            |
#                                                                        |
# AUTHOR: James C. Estill                                                |
# CONTACT: jestill@uga.edu                                               |
# STARTED: 3/9/06                                                        |
# UPDATED: 3/9/06                                                        |
# DESCRIPTION:                                                           |
# Check some of the properties of the sequences files in a fasta database|
# and report sequences that are outside of expected parameters.          |
#                                                                        |
# USAGE:                                                                 |
# seq_clean.pl /home/jestill/BigFastaFile.fasta fasta                     |
#------------------------------------------------------------------------+

# -------------------------
# INCDLUDES
# -------------------------
use Bio::SeqIO;            # Use Bio::SeqIO for the reading and writing of seq files in different formats
use Cwd;                   # Use the Cwd module to get the current working directory

# --------------------------
# VARIABLES
# --------------------------
my $MinLength = "100";       # Minimum seqence length to report


# --------------------------
# CONSTANTS 
# --------------------------
# Get command-line arguments, or die with a usage statement
my $usage = "Fasta2Dir.pl infile infile_format outfile outfile_format \n";
my $infile = shift or die $usage;          # The full path of the infile that will be transformed
my $infileformat = shift or die $usage;    # Above is the format that the input file is in. 
                                           # Any valid bioperl format:
                                           # (ie. abi, ace, gcg, genbank, fasta, swiss, tigr etc.) 
my $outfile = shift or die $usage;         # The output file format
my $outfileformat = shift or die $usage;   # The output file format. Same properties as $infileformat

my $CurrentDir = cwd();                    # Get the current working directory

$SeqNum = 1;

# The sequence input file
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);

# The output file
my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
			      '-format' => $infileformat);

 # write each entry in the input file to the new output file
while (my $inseq = $seq_in->next_seq)
{

      $num = sprintf("%7d", $SeqNum);  # Pad number with zeros so that the total length of the string is 7 characaters
      $num =~ tr/ /0/;

      # --------------------------
      # CHECK SEQUENCE QUALITY 
      # --------------------------
      if ( $inseq->length < $MinLength )
      {
	  print "PROBLEM SEQUENCE RECORD: $num \n";
	  print "Primary ID: ".$inseq->primary_id."\n";
	  print "Length: ".$inseq->length."\n";
      }
      else
      {
	  # If the current sequence checks out, according
	  # to the quality checks above then the sequence
          # record is written to the outfile. 
	  $seq_out->write_seq($inseq); 
      }

      $SeqNum++;

#      # Place to kick out if this is a test run
#      if ($SeqNum == 100)
#      {
#	  exit;
#      }

}

# print "The files have all been placed in the directory: \n";

exit;
