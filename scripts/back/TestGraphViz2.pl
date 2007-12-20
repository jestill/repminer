#!/usr/bin/perl -w
# Jamie's short test of using graph viz to draw graphs
# This will eventually (hopefull very quickly) be added as an
# optional output to RepMiner
##
#Module Docs are at:
#http://search.cpan.org/~lbrocard/GraphViz-2.02/lib/GraphViz.pm
#
# COLORS
# The X11 Color scheme is here:
# http://www.graphviz.org/doc/info/colors.html
# Colors may also be hue, saturation, brightness
# values between 0 and 1

use strict;
use GraphViz;

print "The GraphViz test has started.\n";

#my $g = GraphViz->new();

my $g = GraphViz->new(
		      layout => 'neato', 
		      ratio => 'compress',
		      width => 8.5, 
		      height => 11,
		      directed => 0,
		      overlap => 'false'
		      );

# Test of using variable name
my $Lon = 'London';

# Test Random  nodes
# Num Nodes
my $NumNodes = 100;
#my $i = 1;
my $MaxVal = $NumNodes - 1;
for (my $i=1; $i <= $NumNodes; $i++)
{
    $g->add_node( $i, 
		  shape => 'circle',
		  style => 'filled',
		  color => 'blue',
		  fillcolor => 'red');
}

for (my $i=1; $i <= $MaxVal ; $i++)
{
    $g->add_edge( $i => $i+1);
}
$g->add_edge( 1 => $MaxVal+1);

# Test edges
#$g->add_edge( 1 => 2 );
#$g->add_edge( 1 => 1 );
# Test edges
#$g->add_edge('Augusta' => 'Paris', label => 'Infinite');
#$g->add_edge('London' => 'Paris');
#$g->add_edge('London' => 'New York', label => 'Far');
#$g->add_edge('Paris' => 'London');


# To print the dot file do as_cannon
# print $g->as_cannon;


# This will show the image using ImageMagic
my $png_data = $g->as_png;
open (DISPLAY,"| display -") || die;
binmode DISPLAY;
print DISPLAY $png_data;
close DISPLAY;

#print $g->as_png;

print "The GraphViz test has finished.\n";

exit;
