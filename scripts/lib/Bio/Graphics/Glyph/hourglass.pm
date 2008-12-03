package Bio::Graphics::Glyph::hourglass;
 
use strict;
use base 'Bio::Graphics::Glyph::box';
 
sub draw_component {
  my $self = shift;
  my ($gd,$dx,$dy) = @_;
  my ($left,$top,$right,$bottom) = $self->bounds($dx,$dy);
 
  # draw the hourglass as a polygon
  my $poly = GD::Polygon->new;
  $poly->addPt($left,$top);
  $poly->addPt($right,$bottom);
  $poly->addPt($right,$top);
  $poly->addPt($left,$bottom);
  $poly->addPt($left,$top);
  $gd->filledPolygon($poly,$self->bgcolor);
  $gd->polygon($poly,$self->fgcolor);
}
 
1;
