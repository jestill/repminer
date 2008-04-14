#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# glyph_test.pl - Testing the LTR Retrotransposon Glyph     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/12/2008                                       |
# UPDATED: 03/28/2008                                       |
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

use strict;
use lib "$ENV{HOME}/lib";
 
use Bio::Graphics;
use Bio::SeqFeature::Generic;

#-----------------------------+
# TEST OF LTR RETRO           |
#-----------------------------+

my $bsg = 'Bio::SeqFeature::Generic';
 
my $span         = $bsg->new(-start=>1,-end=>10000);

# The pad_left and pad_right values below can be used to
# add padding to the left and right to allow for the
# seq val and tsd to be drawn onto the canvas
my $panel = Bio::Graphics::Panel->new(-width=>1200,
				      -length=>$span->length,
				      -pad_left=>20,
				      -pad_right=>20);

# Leave out the following to just print the feature
$panel->add_track($span,-glyph=>'arrow',-double=>1,-tick=>2);

# Test of LTR, give start and end of span here
# further information goes below in drawing the track

# All of the following scales are in base pairs
my $test_feature = $bsg->new(-start=>2000,
			     -end=>8000,
			     -key_style => 'between',
                             -display_name=>' Huck',
                             -source_tag=>' This is a test Huck LTR Retro');

$panel->add_track($test_feature,
                  -glyph          => 'ltr_retro',
                  -bgcolor        => 'gray',
		  -fgcolor        => 'gray',
                  -font2color     => 'blue',
		  -height         => 40,
		  -label          => 1,
		  -description    => 1,
		  -label_feat     => 1,          # Label biological features
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
		  # PRIMER BINDING STE
		  -pbs_start      => 3000,       # Primer Binding Site start
		  -pbs_end        => 3020,       # Primer Binding Site end
		  -pbs_fg_color   => 'black',
		  -pbs_bg_color   => 'red',
		  # POLYPURINE TRACT
		  -rr_start       => 6980,       # Polypurine tract start
		  -rr_end         => 7000,       # Polypurine tract end
		  -rr_fg_color    => 'black',    # RR Foreground Color
		  -rr_bg_color    => 'red',      # RR Background Color
		  # ENVELOPE
		  -env_start      => 6700,       # GAG ORF Start
		  -env_end        => 6950,       # GAG ORF End
		  -env_fg_color   => 'black',    # GAG Foreground Color
		  -env_bg_color   => 'purple',      # GAG Background Color
		  # GAG REGION
		  -gag_start      => 3500,       # GAG ORF Start
		  -gag_end        => 4500,       # GAG ORF End
		  -gag_fg_color   => 'black',    # GAG Foreground Color
		  -gag_bg_color   => 'red',      # GAG Background Color
		  # ZF Knuckle, zf_cchc
		  -zf_start       => 4550,
		  -zf_end         => 4600,
		  -zf_fg_color    => 'black',
		  -zf_bg_color    => 'red',
		  # POL
		  -pol_start      => 5000,       # Pol Start
		  -pol_end        => 6700,       # Pol End
		  -pol_fg_color   => 'black',    # Pol Foreground Color
		  -pol_bg_color   => 'red',      # Pol Background Color
		  # INTEGRASE
		  -int_start      => 5700,       # Integrase Start
		  -int_end        => 6000,       # Integrase End
		  -int_fg_color   => 'black',    # Int Foreground Color
		  -int_bg_color   => 'red',      # Int Background Color
		  # REVERSE TRANSCRIPTASE
		  -rt_start       => 6000,       # Reverse transcriptase start
		  -rt_end         => 6300,       # Reverse transcriptase end
		  -rt_fg_color    => 'black',      # RT Foreground Color
		  -rt_bg_color    => 'red',      # RT Background Color
		  # RNASEH
		  -rh_start       => 6310,       # RNAseH Start
		  -rh_end         => 6600,       # RNAseH End
		  -rh_fg_color    => 'black',      # RH Foreground Color
		  -rh_bg_color    => 'red',      # RH Background Color
		  # PROTEASE
		  -pro_start      => 5000,       # Protease start
		  -pro_end        => 5680,       # Protease end
		  -pro_fg_color   => 'black',      # Pro Foreground Color
		  -pro_bg_color   => 'red',      # Pro Background Color
		  # CHRomatin Organisation MOdifier Domain
		  -chromo_start     => 6620,       # chromo start
		  -chromo_end       => 6670,       # chromo end
		  -chromo_fg_color  => 'black',    # chromo foreground color
		  -chromo_bg_color  => 'yellow',      # chromo background color
		  #-----------------------------+
		  # LONG TERMINAL REPEATS       |
		  #-----------------------------+
		  -ltr_fg_color   => 'red',      # LTR Foreground Color
		  -ltr_bg_color   => 'red',      # LRT Background Color
		  #-ltr_arrow      => 0,         # Draw the LTR Arrow, def 1
		  -ltr_arrow_fg_color => 'black',  # Arrow Foreground Color 
		  -ltr_arrow_bg_color => 'black',  # Arrow Background Color
		  -ltr5_start     => 2000,       # Start of the 5' LTR
		  -ltr5_end       => 3000,       # End of the 5' LTR
		  #-ltr5_len       => 1000,      # Can also just give length
		  -ltr3_start     => 7000,       # 3' LTR Start
		  -ltr3_end       => 8000,       # 3' LTR End
		  #-----------------------------+
		  # TG/CA in LTR                |
		  #-----------------------------+
		  -tg_bg_color  => 'gold',       # TG/CA Background Color
		  -ca_bg_color  => 'gold',       # TG/CA Background Color
		  -has_ltr5_tg    => 1,          # 5' LTR Starts with tg
		  -has_ltr5_ca    => 1,          # 5' LTR Starts with ca
		  -has_ltr3_tg    => 1,          # 3' LTR starts with tg
		  -has_ltr3_ca    => 1,          # 3' LTR ends with ca
		  #-ltr3_len       => 1000,
		  #-----------------------------+
		  # TARGET SITE DUPLICATION     |
		  #-----------------------------+
		  -seq_con_color  => 'black',    # Color of context text
		  -seq_con_ltr5   => 'left',     # Sequence context 5' LTR
		  -seq_con_ltr3   => 'right',    # Sequence context 3' LTR
		  -tsd_len        => 3,          # TSD Length, 0 is no TSD
		  -tsd_fg_color   => 'red',      # Default is feat fgcolor
		  -tsd_bg_color   => 'red',      # Default is feat bgcolor
		  );

#-----------------------------+
# CONVER THE OBJECT TO A PNG  |
# FILE AND PRINT TO DISPLAY   |
#-----------------------------+
my $png_data = $panel->png;


# Print to test output file
open (PNGOUT, ">test_image.png") ||
    die "Can not open png output test;\n";
binmode PNGOUT;
print PNGOUT $png_data;
close PNGOUT;

# Print to display
#open (DISPLAY,"| display -") || die;
#binmode DISPLAY;
#print DISPLAY $png_data;
#close DISPLAY;

