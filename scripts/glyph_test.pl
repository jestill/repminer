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
use strict;
use lib "/Library/Perl/darwin/";
 
use Bio::Graphics;
use Bio::SeqFeature::Generic;
my $bsg = 'Bio::SeqFeature::Generic';
 
my $span         = $bsg->new(-start=>1,-end=>1000);
my $test_feature = $bsg->new(-start=>300,-end=>700,
                             -display_name=>'Test Feature',
                             -source_tag=>'This is only a test');
 
my $panel        = Bio::Graphics::Panel->new(-width=>600,-length=>$span->length,
                                             -pad_left=>12,-pad_right=>12);
$panel->add_track($span,-glyph=>'arrow',-double=>1,-tick=>2);
 
$panel->add_track($test_feature,
                  -glyph   => 'hourglass',
                  -bgcolor => 'orange',
                  -font2color => 'red',
                  -height  => 20,
                  -label   => 1,
                  -description => 1);
 
print $panel->png;
