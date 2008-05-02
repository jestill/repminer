#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# glyph_test.pl - Testing the LTR Retrotransposon Glyph     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 04/16/2008                                       |
# UPDATED: 04/16/2008                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Test some simple glyph drawing of an LTR Retrotransposon |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO: Consider writing HTML output in show images
#       Fix how to deal with N for null values .. do not draw
#       .. currently tries to draw then at the 1 position with length 0.
#       .. this obscures the TSD letters.

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use lib "$ENV{HOME}/lib";
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outdir;
my $html_out;                 # HTML File linking to images

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
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
		    "t|html=s"    => \$html_out, # t for html 
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$infile) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$infile);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );

}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# FILE IO                     |
#-----------------------------+
if ($html_out) {
    open( HTML, ">$html_out") ||
	die "Can not open html output file: $html_out\n";
    print HTML "<HTML><HEAD><TITLE>Annotation Output</TITLE>\n".
	"</HEAD><BODY>\n";
}


open (INFILE, "<$infile") ||
    die "Can not open input file: $infile\n";

while (<INFILE>) {
    next if /^\#/;
    chomp;
    my @in_data = split;
    my $name = $in_data[0];
    my $rpn_start = $in_data[1];
    my $rpn_end = $in_data[2];
    my $ltr5_start = $in_data[3];
    my $ltr5_end = $in_data[4];
    my $ltr5_len = $ltr5_end - $ltr5_start;
    my $ltr5_dn_start = $in_data[5];
    my $has_ltr5_tg = 0;
    $has_ltr5_tg = 1 if $ltr5_dn_start =~ "TG";
    my $ltr5_dn_end = $in_data[6];
    my $has_ltr5_ca = 0;
    $has_ltr5_ca = 1 if $ltr5_dn_end =~ "CA";
    my $ltr3_start = $in_data[7];
    my $ltr3_end = $in_data[8];
    my $ltr3_dn_start = $in_data[9];
    my $has_ltr3_tg = 0;
    $has_ltr3_tg = 1 if $ltr3_dn_start =~ "TG";
    my $has_ltr3_ca = 0;
    $has_ltr3_ca = 1 if $ltr5_dn_end =~ "CA";
    my $ltr3_dn_end = $in_data[10];
    my $pbs_start = $in_data[11];
    my $pbs_end = $in_data[12];
    my $rr_start = $in_data[13];
    my $rr_end = $in_data[14];
    my $tsd_5_start = $in_data[15];
    my $tsd_5_end = $in_data[16];
    my $tsd_5_seq = $in_data[17];
    my $tsd_3_stat = $in_data[18];
    my $tsd_3_end = $in_data[19];
    my $tsd_3_seq = $in_data[20];
    my $gag_start = $in_data[21];
    my $gag_end = $in_data[22];
    my $zf_start = $in_data[23];
    my $zf_end = $in_data[24];
    my $rvp_start = $in_data[25];
    my $rvp_end = $in_data[26];
    my $rve_start = $in_data[27];
    my $rve_end = $in_data[28];
    my $rvt_start = $in_data[29];
    my $rvt_end = $in_data[30];
    my $rh_start = $in_data[31];
    my $rh_end = $in_data[32];
    my $chrom_start = $in_data[33];
    my $chrom_end = $in_data[34];
    my $env_start = $in_data[35];
    my $env_end = $in_data[36];
    my $flank_5_seq = $in_data[37];
    my $flank_3_seq = $in_data[38];
    
    print $name."\n";


   #-----------------------------+
   # DRAW LTR RETRO              |
   #-----------------------------+
   # CUrrently using tab delim input file
    
    my $img_out = $name."png";
    my $bsg = 'Bio::SeqFeature::Generic';
    my $span         = $bsg->new(-start=> $rpn_start,
				 -end  => $rpn_end );
    
   # The pad_left and pad_right values below can be used to
   # add padding to the left and right to allow for the
   # seq val and tsd to be drawn onto the canvas
    my $panel = Bio::Graphics::Panel->new(-width=>800,
					  -length=>$span->length,
					  -pad_left=>40,
					  -pad_right=>40);
    
   # Leave out the following to just print the feature
    $panel->add_track($span,-glyph=>'arrow',-double=>1,-tick=>2);

   # Test of LTR, give start and end of span here
   # further information goes below in drawing the track

   # All of the following scales are in base pairs
    my $test_feature = $bsg->new(-start        => $rpn_start,
				 -end          => $rpn_end,
				 -key_style    => 'between',
				 -display_name => $name,
				 -source_tag   => $name);
    
    $panel->add_track($test_feature,
		      -glyph          => 'ltr_retro',
		      -bgcolor        => 'gray',
		      -fgcolor        => 'gray',
		      -font2color     => 'blue',
		      -height         => 40,
		      -label          => 1,
		      -description    => 1,
		      -label_feat     => 1,  # Label biological features
		      # Default 0
		      # gag,pol etc
		      #-label_feat_color = 'black',
		      #-----------------------------+
		      # BASE LTR RETRO SPAN         |
		      #-----------------------------+
		      -span_fg_color  => 'gray',
		      -span_bg_color  => 'gray',
		      #-----------------------------+
		      # CODING REGIONS              |
		      #-----------------------------+
		      # PRIMER BINDING SITE
		      -pbs_start      => $pbs_start, # PBS start
		      -pbs_end        => $pbs_end,   # PBS end
		      -pbs_fg_color   => 'black',
		      -pbs_bg_color   => 'red',
		      # POLYPURINE TRACT
		      -rr_start       => $rr_start,  # Polypurine tract start
		      -rr_end         => $rr_end,    # Polypurine tract end
		      -rr_fg_color    => 'black',    # RR Foreground Color
		      -rr_bg_color    => 'red',      # RR Background Color
		      # ENVELOPE
		      -env_start      => $env_start, # GAG ORF Start
		      -env_end        => $env_end,   # GAG ORF End
		      -env_fg_color   => 'black',    # GAG Foreground Color
		      -env_bg_color   => 'purple',   # GAG Background Color
		      # GAG REGION
		      -gag_start      => $gag_start, # GAG ORF Start
		      -gag_end        => $gag_end,   # GAG ORF End
		      -gag_fg_color   => 'black',    # GAG Foreground Color
		      -gag_bg_color   => 'red',      # GAG Background Color
		      # ZF Knuckle, zf_cchc
		      -zf_start       => $zf_start,
		      -zf_end         => $zf_end,
		      -zf_fg_color    => 'black',
		      -zf_bg_color    => 'red',
#		      # POL
#		      -pol_start      => 5000,       # Pol Start
#		      -pol_end        => 6700,       # Pol End
#		      -pol_fg_color   => 'black',    # Pol Foreground Color
#		      -pol_bg_color   => 'red',      # Pol Background Color
		      # INTEGRASE
		      -int_start      => $rve_start, # Integrase Start
		      -int_end        => $rve_end,   # Integrase End
		      -int_fg_color   => 'black',    # Int Foreground Color
		      -int_bg_color   => 'red',      # Int Background Color
		      # REVERSE TRANSCRIPTASE
		      -rt_start       => $rvt_start, # RT start
		      -rt_end         => $rvt_end,   # RT end
		      -rt_fg_color    => 'black',    # RT Foreground Color
		      -rt_bg_color    => 'red',      # RT Background Color
		      # RNASEH
		      -rh_start       => $rh_start,  # RNAseH Start
		      -rh_end         => $rh_end,    # RNAseH End
		      -rh_fg_color    => 'black',    # RH Foreground Color
		      -rh_bg_color    => 'red',      # RH Background Color
		      # PROTEASE
		      -pro_start      => $rvp_start, # Protease start
		      -pro_end        => $rvp_end,   # Protease end
		      -pro_fg_color   => 'black',    # Pro Foreground Color
		      -pro_bg_color   => 'red',      # Pro Background Color
		      # CHRomatin Organisation MOdifier Domain
		      -chromo_start     => $chrom_start, # chromo start
		      -chromo_end       => $chrom_end,   # chromo end
		      -chromo_fg_color  => 'black',      # chromo fore color
		      -chromo_bg_color  => 'yellow',     # chromo back color
		      #-----------------------------+
		      # LONG TERMINAL REPEATS       |
		      #-----------------------------+
		      -ltr_fg_color   => 'red',      # LTR Foreground Color
		      -ltr_bg_color   => 'red',      # LRT Background Color
		      #-ltr_arrow      => 0,         # Draw LTR Arrow [1] 
		      -ltr_arrow_fg_color => 'black',  # Arrow Fore Color 
		      -ltr_arrow_bg_color => 'black',  # Arrow Back Color
		      -ltr5_start     => $ltr5_start+1,  # Start of the 5' LTR
		      -ltr5_end       => $ltr5_end,    # End of the 5' LTR
		      -ltr3_start     => $ltr3_start,  # 3' LTR Start
		      -ltr3_end       => $ltr3_end,    # 3' LTR End
		      #-----------------------------+
		      # TG/CA in LTR                |
		      #-----------------------------+
		      -tg_bg_color  => 'gold',       # TG/CA Background Color
		      -ca_bg_color  => 'gold',       # TG/CA Background Color
		      -has_ltr5_tg    => $has_ltr5_tg,  # 5' LTR Starts w tg
		      -has_ltr5_ca    => $has_ltr5_ca,  # 5' LTR Ends w ca
		      -has_ltr3_tg    => $has_ltr3_tg,  # 3' LTR starts w tg
		      -has_ltr3_ca    => $has_ltr3_ca,  # 3' LTR ends w ca
		      #-----------------------------+
		      # TARGET SITE DUPLICATION     |
		      #-----------------------------+
		      -seq_con_color  => 'black',    # Color of context text
		      -seq_con_ltr5   => $tsd_5_seq, # Sequence context 5' LTR
		      -seq_con_ltr3   => $tsd_3_seq, # Sequence context 3' LTR
		      -tsd_len        => 3,          # TSD Length, 0 is no TSD
		      -tsd_fg_color   => 'red',      # Default is feat fgcolor
		      -tsd_bg_color   => 'red',      # Default is feat bgcolor
		      );
    
   #-----------------------------+
   # CONVER THE OBJECT TO A PNG  |
   # FILE AND PRINT TO DISPLAY   |
   #-----------------------------+
    my $png_data = $panel->png;

    # MAKE PNG OUT FILE
    open (PNGOUT, ">".$outdir.$img_out.".png") ||
	die "Can not open png output test;\n";
    binmode PNGOUT;
    print PNGOUT $png_data;
    close PNGOUT;

    if ($html_out) {
	# Link to image created above
	print HTML "<a href=".$outdir.$img_out.".png>".$name."</a>|".$ltr5_len."<br>\n";
	unless ( $env_start=~"N") {
	    print HTML "HAS ENV<br>\n";
	}
    }
    
} # End of while infile

if ($html_out) {
    print HTML "</BODY>";
    print HTML "</HTML>";
    close (HTML);

}
exit;


#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 04/16/2008
# - Program started
# 
# 04/18/2008
# - Adding HTML output with -h option
