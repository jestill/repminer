#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# glyph_test.pl - Testing the bioperl glyphy function       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/12/2008                                       |
# UPDATED: 03/12/2008                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Test some simple glyph drawing                           |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

use GD;

my $max_x = 100;
my $max_y = 100;
my $pxs = 5;

my $img = new GD::Image($max_x, $max_y);

#-----------------------------+
# ALLOCATE COLORS FOR THE     |
# REPEAT CATEGORIES           |
#-----------------------------+
my $white = $img->colorAllocate      ( 255,  255,255);
my $colSep = $img->colorAllocate     (   0,   0,   0);
my $colUnk = $img->colorAllocate     (   0,   0,   0);
my $colLTR = $img->colorAllocate     ( 255,   0,   0);
my $colNonLTR = $img->colorAllocate  (   0,   0, 255);
my $colMITE = $img->colorAllocate    (   0, 255, 255);
my $colCACTA = $img->colorAllocate   (   0, 255,   0);
my $colSelf = $img->colorAllocate    ( 200, 200, 200);  
my $colOryza = $img->colorAllocate   (   0, 255,   0);
my $colWess = $img->colorAllocate    ( 255,   0,   0);
my $colSanMig = $img->colorAllocate  (   0,   0, 255);
my $colZea = $img->colorAllocate     ( 200, 200,   0);


my $x_pos = 25;
my $y_pos = 25;

$img->arc($x_pos, $y_pos, $pxs, $pxs, 0, 360, $colSelf);
$img->fill($x_pos, $y_pos, $colNonLTR);

my $png_data = $img->png;
open (DISPLAY,"| display -") || die;
binmode DISPLAY;
print DISPLAY $png_data;
close DISPLAY;
