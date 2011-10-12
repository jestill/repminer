#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# date_ltrs.pl - Estimate divergence times of LTRs          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill, Gina Baucom                     |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/11/2011                                       |
# UPDATED: 10/11/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
# Pipeline code to get dates for LTR retros based on code   |
# originally written by Gina Baucom.                        |
#                                                           |
#-----------------------------------------------------------+

package REPMINER;

# TO DO:
# accept list of conversions to use for converting
# from the distance to an estimated date since
# divergence.

# REQUIRED SOFTWARE
# clustalw
# - required for alignment
# baseml
# - part ofthe paml package

# Limitations
# There are limiations on the fasta file header sizes that the 
# programs used in the pipeline can make use of.
# If intput fasta files contain more than one sequence
# this will give results only for the first two sequences
# in the fasta file, and assume that these are the two
# LTRs
# Requires input as FASTA file with short headers
# should modfy thise to accept a tab delimited
# text file for input

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
my $indir;                      # Input directory of fasta files
my $outdir;                     # Output directory for baseml out
my $distout;                    # Output file for distance
my $tabin;                      # Tab delimited text file (name 5ltr  3ltr)

my $verbose;
my $show_usage;
my $show_version;
my $show_man;
my $show_help;

# Substitution rates to use to estimate date of divergence
my @sub_rates = (0.000000013,     # Grass rate (Ma and Bennetzen, 2008)
		 0.0000000065);   # Half the grass rate               

                               # Previously used divergence dates ...
                               # Rice --
                               # 1.3 x 10^-8 subs per site per year (Ma 
                               # 0.000000013
                               # Maize
                               # 6.5 x 10^-9 original SanMiguel 1998 rate
                               # 1.3 x 10^-8 used to get time (Years)


# It will also be a good idea to allow for the specification
# of the baseml model parameters from the command line
# so that the user does not need to modify the TMPCTL file
# by hand in the source code below.

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "t|tabin=s"    => \$tabin,
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "d|distout=s", => \$distout,
		    # ADDITIONAL OPTIONS
		    "verbose"      => \$verbose,
                     # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# CREATE THE OUTPUT DIR IF IT |
# DOES NOT ALREADY EXIST      |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "ERROR: Could not create the output directory:\n$outdir";
}


# Open the distance output file for writing
if ($distout) {
open (DISTOUT, ">".$distout) ||
    die "Can not open the distance file for output $distout \n";
}

#
# TO DO: ACCEPT INPUT AS TAB DELIMITED
#        TEXT FILE HERE
#        AND WRITE DIRECTORY OF FASTA FILES
#
# 
if ($tabin) {
    print STDERR "Generating fasta files\n";

    # Create the input dir that will hold the fasta files
    unless (-e $indir) {
	print "Creating input dir to hold fasta files ...\n" if $verbose;
	mkdir $indir ||
	    die "ERROR: Could not create the input directory:\n$indir";
    }

    open (TABIN, "<".$tabin) ||
	die "Can not open alignment input $tabin \n";
    my $tabline=0;
    while (<TABIN>) {
	$tabline++;
	chomp;
	my @tab_parts = split(/\t/,$_);
	my $num_tab_parts = @tab_parts;

	# If expected number of columns then write to the output
	# file otherwise report error. Note: no checking for input
	# as sequences.
	if ($num_tab_parts == 3) {
	    
	    my $fasta_file = $tab_parts[0].".fasta";
	    my $fasta_file_path = $indir.$fasta_file;
	    my $ltr5_name = $tab_parts[0]."_5";
	    my $ltr3_name = $tab_parts[0]."_3";
	    
	    # WRITE FASTA FILE
	    # CAN CHECK FOR LENGTH OF FASTA FILE HEADERS HERE
	    # TO MAKE SURE THAT THEY WILL WORK WITH THE DOWNSTREAM
	    # ALIGNMENT AND BASEML.
	    open (FASTOUT, ">".$fasta_file_path) ||
		die "can not open fasta file at $fasta_file_path\n";
	    print FASTOUT ">".$ltr5_name."\n".$tab_parts[1]."\n";
	    print FASTOUT ">".$ltr3_name."\n".$tab_parts[2]."\n";
	    close (FASTOUT);

	}
	else {
	    print STDERR "WARNING: Unexpected number of columns in line $tabline".
		"in input file $tabin\n";
	}
	
    } # End of while TABIN
    close (TABIN);
    
} # End of if tab delim file used as input


#-----------------------------+
# LOAD FASTA FILES            |
#-----------------------------+

# OPEN INPUT DIR TO GET FASTA FILES
opendir( DIR, $indir ) ||
    die "Can't open $indir";

# Look for all *fasta files in input dir 
my @ltr_files = grep /.fasta$/,readdir DIR ;
closedir( DIR );

#-----------------------------+
# PROCESS EACH FASTA FILE     |
#-----------------------------+
for my $ind_file (@ltr_files) {
    
    my $name_root;

    #-----------------------------+
    # Get the root name of the    |
    # fasta file to process       |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }

    #-----------------------------+
    # ALIGN WITH CLUSTALW         |
    #-----------------------------+

    my $ind_file_path = $indir.$ind_file;
    my $clustal_cmd = "clustalw ".$ind_file_path." -output=phylip";
    
    print STDERR "Aligning ".$ind_file."\n" if $verbose;
    
    # Run clustal command
    system ( $clustal_cmd );
    
    # Check for phy alignment file in the input dir
    my $clust_result = $indir.$name_root.".phy";
    if (-e $clust_result) {
	print STDERR "Found alignment file\n";
    }
    else {
	print "\a";
	print STDERR "WARNING: An alignment file was not found".
	    " for the input file: ".$ind_file."\n";
	# Skip to the next input file in the indir
	next;
    }

 
    #-----------------------------+
    # ADD I TO THE ALIGNMENT      |
    # RESULT                      |
    #-----------------------------+
    # I believe that the I tells the program
    # that this is an interleaved format
    # Gina's program.
    my $phy_in = $clust_result;
    my $phy_mod = $indir.$name_root."_i.phy";
    
    open (PHYIN, "<".$phy_in) ||
	die "Can not open alignment input $phy_in \n";
    open (PHYOUT, ">".$phy_mod) ||
	die "Can not open alignment output $phy_mod \n";
    
    my $line_num = 0;
    while (<PHYIN>) {
	chomp;
	$line_num++;

	if ($line_num == 1) {
	    print PHYOUT $_."\tI\n";
	}
	else {
	    print PHYOUT $_."\n";
	}

    }

    close (PHYIN);
    close (PHYOUT);

    # Check that the output file was actually created
    if (-e $clust_result) {
	print STDERR "Found modified alignment file\n" 
	    if $verbose;
    }
    else {
	print "\a";
	print STDERR "WARNING: A modified alignment file was not found".
	    " for the input file: ".$ind_file."\n";
	# Skip to the next input file in the indir
	next;
    }


    #-----------------------------+
    # RUN baseml ON THE MODIFIED  |       
    # ALIGNMENT RESULT            |
    #-----------------------------+

    my $treefile = $indir.$name_root.".dnd";
    my $seqfile = $phy_mod;
    my $baseml_out = $outdir."mlc.".$name_root;
    
    open (TMPCTL, ">baseml.ctl") ||
	die "Can not open the baseml.ctl file for output\n";

    # The following would be the place to allow for a modification of the
    # vars used by baseml
    print TMPCTL "seqfile = ".$seqfile."\n".
	"treefile = ".$treefile."\n".
#	"outfile = ".$outdir."mlc.".$name_root."\n\n".
	"outfile = ".$baseml_out."\n\n".
	"noisy = 0       * 0,1,2,3,9: how much rubbish on the screen\n".
	"verbose = 0     * 1: detailed output, 0: concise output\n".
	"runmode = 0     * 0: user tree;  1: semi-automatic;  2: automatic\n".
	"                * 3: StepwiseAddition; (4,5):PerturbationNNI\n\n".
	# * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
	"model = 5       * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n".
	"                * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n".
	"Mgene = 0       * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n".
#	"*ndata = 1\n".
	"clock = 0       * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n".
	"fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below\n".
	"kappa = 5       * initial or fixed kappa\n".
	"fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below\n".
	"alpha = 1       * initial or fixed alpha, 0:infinity (constant rate)\n".
	"Malpha = 0      * 1: different alpha's for genes, 0: one alpha\n".
	"ncatG = 5       * # of categories in the dG, AdG, or nparK models of rates\n".
	"nparK = 0       * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK\n".
	"nhomo = 0       * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n".
	"getSE = 1       * 0: don't want them, 1: want S.E.s of estimates\n".
	"RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states\n".
	"Small_Diff = 7e-6\n".
	"cleandata = 1   * remove sites with ambiguity data (1:yes, 0:no)?\n".
#	"icode = 0  * (with RateAncestor=1. try GC in data,model=4,Mgene=4)\n ".
	"fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n".
	"method = 0  * 0: simultaneous; 1: one branch at a time\n";
	
    close (TMPCTL);

    print STDERR "Running baseml for $ind_file\n" if $verbose;

     # this runs baseml, which will just blindly load and run 
     # the baseml.ctl file we just printed above
    system ("baseml");

    # Once complete, remove the temp files
#    system ("rm baseml.ctl") if (-e "baseml.ctl");
    system ("rm 2base.t*") if (-e "2base.t*");
    system ("rm lnf") if (-e "lnf");
    
    # Check that the output file exists
    if (-e $baseml_out) {
	print STDERR "Found baseml result file\n".$baseml_out."\n" 
	    if $verbose;
    }
    else {
	print "\a";
	print STDERR "WARNING: The baseml output file ".
	    " was not found\n: ".$baseml_out."\n";
	next;
    }

    #-----------------------------+
    # FETCH THE DISTANCE FROM THE |
    # BASEML RESULT FILE AND      |
    # PRINT TO DISTOUT FILE       |
    #-----------------------------+
    open (BASEML, "<".$baseml_out) ||
	die "Can not open baseml result $baseml_out \n";

    my $in_distance = 0;
    while (<BASEML>) {

	chomp;

	# Set boolean for being in distance 
	# part of the output file
	if (m/^Distances/) {
	    $in_distance = 1;
	}
	elsif (m/^TREE #/) {
	    $in_distance = 0;
	}

	# print intput line if in distance 
	# part of the output file
	if ($in_distance) {
	    print STDERR $_."\n";
	    # Use regexp to fetch the distance
	    if (m/\s(\S*)\(/) {
		my $distance = $1;

		# Print the distance out to the output fil
		print STDERR "======\n".$name_root."\t".$distance."\n"
		    if $verbose;

		#-----------------------------+
		# DETERMINE DIVERGENCE TIME   |
 		# AND PRINT OUTPUT            |
		#-----------------------------+
		if ($distance) {
		    print DISTOUT $name_root."\t".$distance."\t";

		    # Do the math for each substitution rate
		    for my $ind_rate (@sub_rates) {
			# Do division and convert to an integer
			my $div_time = int( ( $distance / (2 * $ind_rate) ) );
			
			print DISTOUT $div_time."\t";
		    }

		    print DISTOUT "\n";
		
		}
		
	    }
	}
	
    }

    close (BASEML);

}


close (DISTOUT);

exit;

