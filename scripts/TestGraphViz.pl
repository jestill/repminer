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

# Test nodes
$g->add_node($Lon, 
	     shape => 'circle',
	     color => 'blue',
	     fillcolor => 'red',
	     style => 'filled');
$g->add_node('Paris', 
	     label => 'City of\nlurve', 
	     shape => 'circle',
	     style => 'filled',
	     fillcolor => 'cyan3');
$g->add_node('New York');
$g->add_node('Augusta');






# Test edges
$g->add_edge('Augusta' => 'Paris', label => 'Infinite');
$g->add_edge('London' => 'Paris');
$g->add_edge('London' => 'New York', label => 'Far');
$g->add_edge('Paris' => 'London');

# This will show the image using ImageMagic
my $png_data = $g->as_png;
open (DISPLAY,"| display -") || die;
binmode DISPLAY;
print DISPLAY $png_data;
close DISPLAY;




#print $g->as_png;

print "The GraphViz test has finished.\n";

exit;
