package Bio::Graphics::Glyph::ltr_retro;

#-----------------------------------------------------------+
# LTR Retrotransposon Object
#-----------------------------------------------------------+
#
#  AUTHOR: James C. Estill
# STARTED: 03/13/2008
# UPDATED: 03/28/2008
# DESCRIPTION:
#  Glyph for a  highly detailed LTR Retrotransposon. This
#  is an object for use with the BioPERL drawing library.
#
#-----------------------------------------------------------+
#
# TO DO:
# - Modify the 0.05, 0.2 etc values to variables
#   that could be passed to the glyph. These could 
#   simply be height values ranging from 0 to 1 and
#   the values used to offset are $given_value/2
# - Switch half circles for tg/ac and arrows to point in
#   correct direction for minus strand orientation
# - Modify the glyph draw to take into account pixel height
#   when drawing LTR Retros, this will make sure that these
#   look okay all the way down to 3pixels in height. 
# - Modify the features to be levels if possible
# - When the named color passed to the feature is nonsense,
#   default to the base background and foreground color instead
#   of the odd gray matter color
# - Add sequence font as a variable accepting
#   tiny,small,medium,large,giant
#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
# Since box is for simple features, I will want to use
# segmets since ltr retros have multiple segments
#use base 'Bio::Graphics::Glyph::segments';
use base 'Bio::Graphics::Glyph::box';

#-----------------------------+
# CONSTANTS                   |
#-----------------------------+
use constant TSD_LEN => 0;

# CONTEXT FONT
# FONT CHOICES LIMITED TO
#   FONT NAME       WIDTH(pxs)    HEIGHT(pxs)
#   gdTinyFont           5            8
#   gdSmallFont          6           12
#   gdMediumBoldFont     7           12
#   gdLargeFont          8           16
#   gdGiantFont          9           15
# This assumes that we are using the GD image class
# SVG Would probably not work with the following font choice
# It would be more stable to get font_class first as is used in
# Bio::Graphics::Panel.pm. In Panel.pm this is done as
#  $gd->string($self->{key_font},$left,KEYPADTOP+$top,"KEY:",$text_color);
# However this works for me for now

my $con_font = GD::gdSmallFont;
#my $feat_font = GD::gdSmallFont;
my $feat_font = GD::gdMediumBoldFont;
#my $con_font = GD::gdTinyFont;
#my $con_font = GD::gdGiantFont;

#-----------------------------+
# PADDING                     |
#-----------------------------+
# TSDs are used to PAD since these are set to be outside
# of the LTR Retrotransposon feature
# FROM http://www.perlmonks.org/?node_id=66948
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}

sub pad_left {

    # Left needs to be padded by the drawn length of the
    # target site duplication as well as the length of 
    # the sequence context string
    my $self = shift;
    my $pad_left = 1;
    my $req_tsd_len = $self->option('tsd_len');
    my $req_seq_context = $self->seq_con_ltr5;
    
    if ($self->tsd_len) {
	# Since this can give fractional values, I will always
	# round this up using the internal roundup function
	$pad_left = $pad_left  + ($self->tsd_len*$self->scale);
	$pad_left = roundup($pad_left);
    }
    
    if ($self->seq_con_ltr5) {
	my $char_width = $con_font->width;
	my $seq_width = length($self->seq_con_ltr5) * $char_width;
	$pad_left = $pad_left + ($seq_width);
    }

    return $pad_left;
    
}

sub pad_right {
    my $self = shift;
    my $pad_right = 1;
    my $req_tsd_len = $self->option('tsd_len');
    my $req_seq_context = $self->seq_con_ltr3;
    
    if ($self->tsd_len) {
	# Since this can give fractional values, I will always
	# round this up using the internal roundup function
	$pad_right = $pad_right  + ($self->tsd_len*$self->scale);
	$pad_right = roundup($pad_right);
    }

    if ($self->seq_con_ltr3) {
	my $char_width = $con_font->width;
	my $seq_width = length($self->seq_con_ltr3) * $char_width;
	$pad_right = $pad_right + ($seq_width);
    }

    return $pad_right;

}

#-----------------------------+
# BASE LTR RETRO SPAN         |
#-----------------------------+
# This is the bare minimum that will be drawn
sub span_fg_color {
    my $self = shift;
    my $req_span_fg_color = $self->option('span_fg_color');
    if (defined $req_span_fg_color) {
	return $self->translate_color($req_span_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub span_bg_color {
    my $self = shift;
    my $req_span_bg_color = $self->option('span_bg_color');
    if (defined $req_span_bg_color) {
	return $self->translate_color($req_span_bg_color);
	
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# LABEL FEATURES              |
#-----------------------------+
sub label_feat {
    my $self = shift;
    my $req_label_feat = $self->option('label_feat');
    if (defined $req_label_feat) {
	return $req_label_feat;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# SEQUENCE CONTEXT            |
#-----------------------------+
# The purpose of printing the sequence context is to show the putative TSDs
# on either side of the LTR. This could also be used to draw a larger context
# sequence such as that returned by LTR_Struc
sub seq_con_ltr5 {
    my $self = shift;
    my $req_seq_con_ltr5 = $self->option('seq_con_ltr5');
    if (defined $req_seq_con_ltr5) {
	return $req_seq_con_ltr5;
    }
    else {
	return 0;
    }
}

sub seq_con_ltr3 {
    my $self = shift;
    my $req_seq_con_ltr3 = $self->option('seq_con_ltr3');
    if (defined $req_seq_con_ltr3) {
	return $req_seq_con_ltr3;
    }
    else {
	return 0;
    }
}

sub seq_con_color {
    my $self = shift;
    my $req_seq_con_color = $self->option('seq_con_color');
    if (defined $req_seq_con_color) {
	return $self->translate_color($req_seq_con_color);
    }
    else {
	return $self->fgcolor;
    }
}

#-----------------------------+
# TARGET SITE DUPLICATION     |
#-----------------------------+
sub tsd_len {
    my $self = shift;
    my $requested_gag = $self->option('tsd_len');
    if (defined $requested_gag) {
	return $requested_gag;
    }
    else {
	return TSD_LEN;
    }
}

sub tsd_fg_color {
    my $self = shift;
    my $req_tsd_fg_color = $self->option('tsd_fg_color');
    if (defined $req_tsd_fg_color) {
	return $self->translate_color($req_tsd_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub tsd_bg_color {
    my $self = shift;
    my $req_tsd_bg_color = $self->option('tsd_bg_color');
    if (defined $req_tsd_bg_color) {
	return $self->translate_color($req_tsd_bg_color);
	
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# GENERAL TG/CA COLORS        |
#-----------------------------+
sub tg_bg_color {
    my $self = shift;
    my $req_tg_bg_color = $self->option('tg_bg_color');
    if (defined $req_tg_bg_color) {
	return $self->translate_color($req_tg_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub ca_bg_color {
    my $self = shift;
    my $req_ca_bg_color = $self->option('ca_bg_color');
    if (defined $req_ca_bg_color) {
	return $self->translate_color($req_ca_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# 5' LONG TERMINAL REPEAT     |
#-----------------------------+
sub ltr_fg_color {
    my $self = shift;
    my $req_ltr_fg_color = $self->option('ltr_fg_color');
    if (defined $req_ltr_fg_color) {
	return $self->translate_color($req_ltr_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub ltr_bg_color {
    my $self = shift;
    my $req_ltr_bg_color = $self->option('ltr_bg_color');
    if (defined $req_ltr_bg_color) {
	return $self->translate_color($req_ltr_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub ltr5_start {
    my $self = shift;
    my $req_ltr5_start = $self->option('ltr5_start');
    if (defined $req_ltr5_start) {
#	if ($req_ltr5_start =~ "N") { # Treat N as NOT a value
#	    return 0;
#	}
#	else {
	    return $req_ltr5_start;
#	}
    }
    else {
	return 0;
    }
}

sub ltr5_end {
    my $self = shift;
    my $req_ltr5_end = $self->option('ltr5_end');
    if (defined $req_ltr5_end) {
	return $req_ltr5_end;
    }
    else {
	return 0;
    }
    
}

# Alternative to providing the start and stop of the ltr feature
sub ltr5_len {
    my $self = shift;
    my $req_ltr5_len = $self->option('ltr5_len');
    if (defined $req_ltr5_len) {
	return $req_ltr5_len;
    }
    else {
	return 0;
    }
    
}

sub has_ltr5_tg {
    my $self = shift;
    my $req_has_ltr5_tg = $self->option('has_ltr5_tg');
    if (defined $req_has_ltr5_tg) {
	return $req_has_ltr5_tg;
    }
    else {
	return 0;
    }
}

sub has_ltr5_ca {
    my $self = shift;
    my $req_has_ltr5_ca = $self->option('has_ltr5_ca');
    if (defined $req_has_ltr5_ca) {
	return $req_has_ltr5_ca;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# 3' LONG TERMINAL REPEAT     |
#-----------------------------+
sub ltr3_start {
    my $self = shift;
    my $req_ltr3_start = $self->option('ltr3_start');
    if (defined $req_ltr3_start) {
	return $req_ltr3_start;
    }
    else {
	return 0;
    }
}

sub ltr3_end {
    my $self = shift;
    my $req_ltr3_end = $self->option('ltr3_end');
    if (defined $req_ltr3_end) {
	return $req_ltr3_end;
    }
    else {
	return 0;
    }
    
}

# Altrenative is to just provide the length of the ltr
sub ltr3_len {
    my $self = shift;
    my $req_ltr3_len = $self->option('ltr3_len');
    if (defined $req_ltr3_len) {
	return $req_ltr3_len;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# LTR TG AND CA               |
#-----------------------------+
sub has_ltr3_tg {
    my $self = shift;
    my $req_has_ltr3_tg = $self->option('has_ltr3_tg');
    if (defined $req_has_ltr3_tg) {
	return $req_has_ltr3_tg;
    }
    else {
	return 0;
    }
}

sub has_ltr3_ca {
    my $self = shift;
    my $req_has_ltr3_ca = $self->option('has_ltr3_ca');
    if (defined $req_has_ltr3_ca) {
	return $req_has_ltr3_ca;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# LTR ARROW                   |
#-----------------------------+
sub ltr_arrow_fg_color {
    my $self = shift;
    my $req_ltr_arrow_fg_color = $self->option('ltr_arrow_fg_color');
    if (defined $req_ltr_arrow_fg_color) {
	return $self->translate_color($req_ltr_arrow_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub ltr_arrow_bg_color {
    my $self = shift;
    my $req_ltr_arrow_bg_color = $self->option('ltr_arrow_bg_color');
    if (defined $req_ltr_arrow_bg_color) {
	return $self->translate_color($req_ltr_arrow_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub ltr_arrow {
    my $self = shift;
    my $req_ltr_arrow = $self->option('ltr_arrow');
    if (defined $req_ltr_arrow) {
	return $req_ltr_arrow;
    }
    else {
	return 1;
    }
}

#-----------------------------+
# PRIMER BINDING SITE         |
#-----------------------------+
sub pbs_fg_color {
    my $self = shift;
    my $req_pbs_fg_color = $self->option('pbs_fg_color');
    if (defined $req_pbs_fg_color) {
	return $self->translate_color($req_pbs_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub pbs_bg_color {
    my $self = shift;
    my $req_pbs_bg_color = $self->option('pbs_bg_color');
    if (defined $req_pbs_bg_color) {
	return $self->translate_color($req_pbs_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub pbs_start {
    my $self = shift;
    my $req_pbs_start = $self->option('pbs_start');
    if (defined $req_pbs_start) {
	return $req_pbs_start;
    }
    else {
	return 0;
    }
}

sub pbs_end {
    my $self = shift;
    my $req_pbs_end = $self->option('pbs_end');
    if (defined $req_pbs_end) {
	return $req_pbs_end;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# POLYPURINE TRACT            |
#-----------------------------+
sub rr_fg_color {
    my $self = shift;
    my $req_rr_fg_color = $self->option('rr_fg_color');
    if (defined $req_rr_fg_color) {
	return $self->translate_color($req_rr_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub rr_bg_color {
    my $self = shift;
    my $req_rr_bg_color = $self->option('rr_bg_color');
    if (defined $req_rr_bg_color) {
	return $self->translate_color($req_rr_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub rr_start {
    my $self = shift;
    my $req_rr_start = $self->option('rr_start');
    if (defined $req_rr_start) {
	return $req_rr_start;
    }
    else {
	return 0;
    }
}

sub rr_end {
    my $self = shift;
    my $req_rr_end = $self->option('rr_end');
    if (defined $req_rr_end) {
	return $req_rr_end;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# ENVELOPE                    |
#-----------------------------+
sub env_fg_color {
    my $self = shift;
    my $req_env_fg_color = $self->option('env_fg_color');
    if (defined $req_env_fg_color) {
	return $self->translate_color($req_env_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub env_bg_color {
    my $self = shift;
    my $req_env_bg_color = $self->option('env_bg_color');
    if (defined $req_env_bg_color) {
	return $self->translate_color($req_env_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub env_start {
    my $self = shift;
    my $req_env_start = $self->option('env_start');
    if (defined $req_env_start) {
	if ($req_env_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_env_start;
	}
    }
    else {
	return 0;
    }
}

sub env_end {
    my $self = shift;
    my $req_env_end = $self->option('env_end');
    if (defined $req_env_end) {
	if ($req_env_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_env_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# GAG                         |
#-----------------------------+
sub gag_fg_color {
    my $self = shift;
    my $req_gag_fg_color = $self->option('gag_fg_color');
    if (defined $req_gag_fg_color) {
	return $self->translate_color($req_gag_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub gag_bg_color {
    my $self = shift;
    my $req_gag_bg_color = $self->option('gag_bg_color');
    if (defined $req_gag_bg_color) {
	return $self->translate_color($req_gag_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub gag_start {
    my $self = shift;
    my $req_gag_start = $self->option('gag_start');
    if (defined $req_gag_start) {
	if ($req_gag_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_gag_start;
	}
    }
    else {
	return 0;
    }
}

sub gag_end {
    my $self = shift;
    my $req_gag_end = $self->option('gag_end');
    if (defined $req_gag_end) {
	if ($req_gag_end =~ "N") {
	    return 0;
	} 
	else {
	    return $req_gag_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# Zinc Finger CCHC            |
#-----------------------------+
sub zf_fg_color {
    my $self = shift;
    my $req_zf_fg_color = $self->option('zf_fg_color');
    if (defined $req_zf_fg_color) {
	return $self->translate_color($req_zf_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub zf_bg_color {
    my $self = shift;
    my $req_zf_bg_color = $self->option('zf_bg_color');
    if (defined $req_zf_bg_color) {
	return $self->translate_color($req_zf_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub zf_start {
    my $self = shift;
    my $req_zf_start = $self->option('zf_start');
    if (defined $req_zf_start) {
	if ($req_zf_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_zf_start;
	}
    }
    else {
	return 0;
    }
}

sub zf_end {
    my $self = shift;
    my $req_zf_end = $self->option('zf_end');
    if (defined $req_zf_end) {
	if ($req_zf_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_zf_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# POL                         |
#-----------------------------+
sub pol_fg_color {
    my $self = shift;
    my $req_pol_fg_color = $self->option('pol_fg_color');
    if (defined $req_pol_fg_color) {
	return $self->translate_color($req_pol_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub pol_bg_color {
    my $self = shift;
    my $req_pol_bg_color = $self->option('pol_bg_color');
    if (defined $req_pol_bg_color) {
	return $self->translate_color($req_pol_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub pol_start {
    my $self = shift;
    my $req_pol_start = $self->option('pol_start');
    if (defined $req_pol_start) {
	if ($req_pol_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_pol_start;
	}
    }
    else {
	return 0;
    }
}

sub pol_end {
    my $self = shift;
    my $req_pol_end = $self->option('pol_end');
    if (defined $req_pol_end) {
	if ($req_pol_end =~ "N") {
	    return 0;
	}
	else {    
	    return $req_pol_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# PROTEASE                    |
#-----------------------------+
sub pro_fg_color {
    my $self = shift;
    my $req_pro_fg_color = $self->option('pro_fg_color');
    if (defined $req_pro_fg_color) {
	return $self->translate_color($req_pro_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub pro_bg_color {
    my $self = shift;
    my $req_pro_bg_color = $self->option('pro_bg_color');
    if (defined $req_pro_bg_color) {
	return $self->translate_color($req_pro_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub pro_start {
    my $self = shift;
    my $req_pro_start = $self->option('pro_start');
    if (defined $req_pro_start) {
	if ($req_pro_start =~ "N" ) {
	    return 0;
	} 
	else {
	    return $req_pro_start;
	}
    }
    else {
	return 0;
    }
}

sub pro_end {
    my $self = shift;
    my $req_pro_end = $self->option('pro_end');
    if (defined $req_pro_end) {
	if ($req_pro_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_pro_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# INTEGRASE                   |
#-----------------------------+
sub int_fg_color {
    my $self = shift;
    my $req_int_fg_color = $self->option('int_fg_color');
    if (defined $req_int_fg_color) {
	return $self->translate_color($req_int_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub int_bg_color {
    my $self = shift;
    my $req_int_bg_color = $self->option('int_bg_color');
    if (defined $req_int_bg_color) {
	return $self->translate_color($req_int_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub int_start {
    my $self = shift;
    my $req_int_start = $self->option('int_start');
    if (defined $req_int_start) {
	if ($req_int_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_int_start;
	}
    }
    else {
	return 0;
    }
}

sub int_end {
    my $self = shift;
    my $req_int_end = $self->option('int_end');
    if (defined $req_int_end) {
	if ($req_int_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_int_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# REVERSE TRANSCRIPTASE       |
#-----------------------------+
sub rt_fg_color {
    my $self = shift;
    my $req_rt_fg_color = $self->option('rt_fg_color');
    if (defined $req_rt_fg_color) {
	return $self->translate_color($req_rt_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub rt_bg_color {
    my $self = shift;
    my $req_rt_bg_color = $self->option('rt_bg_color');
    if (defined $req_rt_bg_color) {
	return $self->translate_color($req_rt_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub rt_start {
    my $self = shift;
    my $req_rt_start = $self->option('rt_start');
    if (defined $req_rt_start) {
	if ($req_rt_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_rt_start;
	}
    }
    else {
	return 0;
    }
}

sub rt_end {
    my $self = shift;
    my $req_rt_end = $self->option('rt_end');
    if (defined $req_rt_end) {
	if ($req_rt_end = "N") {
	    return 0;
	}
	else {
	    return $req_rt_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# CHRromatin Organiation      |
# MOdifier Domain             |
#-----------------------------+
sub chromo_fg_color {
    my $self = shift;
    my $req_chromo_fg_color = $self->option('chromo_fg_color');
    if (defined $req_chromo_fg_color) {
	return $self->translate_color($req_chromo_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub chromo_bg_color {
    my $self = shift;
    my $req_chromo_bg_color = $self->option('chromo_bg_color');
    if (defined $req_chromo_bg_color) {
	return $self->translate_color($req_chromo_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub chromo_start {
    my $self = shift;
    my $req_chromo_start = $self->option('chromo_start');
    if (defined $req_chromo_start) {
	if ($req_chromo_start =~ "N") {
	    return 0;
	}
	else {
	    return $req_chromo_start;
	}
    }
    else {
	return 0;
    }
}

sub chromo_end {
    my $self = shift;
    my $req_chromo_end = $self->option('chromo_end');
    if (defined $req_chromo_end) {
	if ($req_chromo_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_chromo_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# RNASEH                      |
#-----------------------------+
sub rh_fg_color {
    my $self = shift;
    my $req_rh_fg_color = $self->option('rh_fg_color');
    if (defined $req_rh_fg_color) {
	return $self->translate_color($req_rh_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub rh_bg_color {
    my $self = shift;
    my $req_rh_bg_color = $self->option('rh_bg_color');
    if (defined $req_rh_bg_color) {
	return $self->translate_color($req_rh_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

sub rh_start {
    my $self = shift;
    my $req_rh_start = $self->option('rh_start');
    if (defined $req_rh_start) {
	if ($req_rh_start =~ "N") {
	    return 0;
	} 
	else {
	    return $req_rh_start;
	}
    }
    else {
	return 0;
    }
}

sub rh_end {
    my $self = shift;
    my $req_rh_end = $self->option('rh_end');
    if (defined $req_rh_end) {
	if ($req_rh_end =~ "N") {
	    return 0;
	}
	else {
	    return $req_rh_end;
	}
    }
    else {
	return 0;
    }
}

#-----------------------------+
# DRAW THE LTR RETRO          |
#-----------------------------+
sub draw_component {

    my $self = shift;
    my ($gd,$dx,$dy) = @_;

    # Get the bounding box of the glyph
    my ($left,$top,$right,$bottom) = $self->bounds($dx,$dy);
    
    # Get the height
    my $height = $bottom - $top + 1;
    
    # Get the colors
    my $back_color = $self->bgcolor;
    my $fore_color = $self->fgcolor;

    #-----------------------------+
    # BASE LTR RETRO SPAN         |
    #-----------------------------+
    # This used the defauult background and foreground
    # color, it may be useful to make this a varaible
    # outside of the base.
    my $poly = GD::Polygon->new;
    my $span_top = $top+($height * 0.25);
    my $span_bottom = $bottom - ($height * 0.25);
    $poly->addPt($left,$span_top);
    $poly->addPt($right,$span_top);
    $poly->addPt($right,$span_bottom);
    $poly->addPt($left,$span_bottom);
    $poly->addPt($left,$span_top);

    $gd->filledPolygon($poly,$self->span_bg_color);
    $gd->polygon($poly,$self->span_fg_color);

    #-----------------------------+
    # 5' LTR                      |
    #-----------------------------+
    if ( ($self->ltr5_start) && ($self->ltr5_end) ) {
	
	my $ltr5_start_offset_len_bp = $self->ltr5_start - $self->start;
	my $ltr5_start_offset_len_pxs = 
	    $ltr5_start_offset_len_bp * $self->scale;
	my $ltr5_start_pxs = $left + $ltr5_start_offset_len_pxs;
	
	my $ltr5_end_offset_len_bp = $self->ltr5_end - $self->start;
	my $ltr5_end_offset_len_pxs = $ltr5_end_offset_len_bp * $self->scale;
	my $ltr5_end_pxs = $left + $ltr5_end_offset_len_pxs;
	
	$self->filled_box($gd,
			  $ltr5_start_pxs,
			  $top+($height * 0.05 ),
			  $ltr5_end_pxs,
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);
	
	#/////////////////////////////////////////////////////////////////
	# The following would only work for top strand LTRs
	# Test of draw triangle
	if ($self->ltr_arrow) {
	    my $ltr5_tri = GD::Polygon->new;
	    $ltr5_tri->addPt($ltr5_start_pxs, $top+($height * 0.05 ));  # PT A
	    $ltr5_tri->addPt($ltr5_end_pxs, $top+($height * 0.5) );     # PT B
	    $ltr5_tri->addPt($ltr5_start_pxs, $bottom-($height * 0.05 )); # PT C
	    $ltr5_tri->addPt($ltr5_start_pxs, $top+($height * 0.05  ));   # PT A
	    
	    $gd->filledPolygon($ltr5_tri,$self->ltr_arrow_bg_color);
	    $gd->polygon($ltr5_tri,$self->ltr_arrow_fg_color);

	    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	} # End of draw LTR Arrow

	# HAS 5' TA, DRAW HALF CIRCLE
	if ($self->has_ltr5_tg) {
	
	    my $arc_x = $ltr5_start_pxs;
	    my $arc_y = $top+($height * 0.5);
	    my $arc_width = $height * 0.25;
	    my $arc_height = $height * 0.25;
	    my $arc_start = 270;
	    my $arc_end = 90;

	    $gd->filledArc($arc_x,$arc_y,$arc_width,$arc_height,
			   $arc_start,$arc_end, $self->tg_bg_color);
	}

	# HAS 5' CA, DRAW HALF CIRCLE
	if ($self->has_ltr5_ca) {

	    my $arc_x = $ltr5_end_pxs;
	    my $arc_y = $top+($height * 0.5);
	    my $arc_width = $height * 0.25;
	    my $arc_height = $height * 0.25;
	    my $arc_start = 90;
	    my $arc_end = 270;

	    # For circle
	    $gd->filledArc($arc_x,$arc_y,$arc_width,$arc_height,
		     $arc_start,$arc_end,$self->ca_bg_color);
	}

    }
    elsif ( $self->ltr5_len) {
	
	# ALTERNATIVE IS TO JUST PASS THE LTR LENGTH
	# ASSUMING THAT LTR STARTS AT THE LEFT
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $left + ($self->ltr5_len * $self->scale),
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);
    }

    #-----------------------------+
    # 3' LTR                      |
    #-----------------------------+
    if ( ($self->ltr3_start) && ($self->ltr3_end) ) {
	
	my $ltr3_start_offset_len_bp = $self->ltr3_start - $self->start;
	my $ltr3_start_offset_len_pxs = 
	    $ltr3_start_offset_len_bp * $self->scale;
	my $ltr3_start_pxs = $left + $ltr3_start_offset_len_pxs;
	
	my $ltr3_end_offset_len_bp = $self->ltr3_end - $self->start;
	my $ltr3_end_offset_len_pxs = $ltr3_end_offset_len_bp * $self->scale;
	my $ltr3_end_pxs = $left + $ltr3_end_offset_len_pxs;
	
	$self->filled_box($gd,
			  $ltr3_start_pxs,
			  $top+($height * 0.05 ),
			  $ltr3_end_pxs,
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);
		
	#/////////////////////////////////////////////////////////////////
	# The following would only work for top strand LTRs
	# Test of draw triangle
	if ($self->ltr_arrow) {
	    my $ltr3_tri = GD::Polygon->new;
	    $ltr3_tri->addPt($ltr3_start_pxs, $top+($height * 0.05 ));  # PT A
	    $ltr3_tri->addPt($ltr3_end_pxs, $top+($height * 0.5) );     # PT B
	    $ltr3_tri->addPt($ltr3_start_pxs, $bottom-($height * 0.05 )); # PT C
	    $ltr3_tri->addPt($ltr3_start_pxs, $top+($height * 0.05  ));   # PT A
	    
	    $gd->filledPolygon($ltr3_tri,$self->ltr_arrow_bg_color);
	    $gd->polygon($ltr3_tri,$self->ltr_arrow_fg_color);
	    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	}

	if ($self->has_ltr3_tg) {
	
	    my $arc_x = $ltr3_start_pxs;
	    my $arc_y = $top+($height * 0.5);
	    my $arc_width = $height * 0.25;
	    my $arc_height = $height * 0.25;
	    my $arc_start = 270;
	    my $arc_end = 90;

	    $gd->filledArc($arc_x,$arc_y,$arc_width,$arc_height,
			   $arc_start,$arc_end,$self->tg_bg_color);

	}

	if ($self->has_ltr3_ca) {

	    my $arc_x = $ltr3_end_pxs;
	    my $arc_y = $top+($height * 0.5);
	    my $arc_width = $height * 0.25;
	    my $arc_height = $height * 0.25;
	    my $arc_start = 90;
	    my $arc_end = 270;

	    $gd->filledArc($arc_x,$arc_y,$arc_width,$arc_height,
		     $arc_start,$arc_end,$self->ca_bg_color);
	}



    }
    elsif ( $self->ltr3_len) {
	
	$self->filled_box($gd,
			  $right - ($self->ltr3_len * $self->scale),
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);
    }
    
    #-----------------------------+
    # PRIMER BINDING SITE         |
    #-----------------------------+
    # Will onlty draw a PBS when both a start and end location
    # are passed
    if ( ($self->pbs_start) && ($self->pbs_end) ) {

	# Given the start and stop of the PBS in base pairs
	# convert the feature locations to pixels and draw
	my $pbs_start_offset_len_bp = $self->pbs_start - $self->start;
	my $pbs_start_offset_len_pxs = $pbs_start_offset_len_bp * $self->scale;
	my $pbs_start_pxs = $left + $pbs_start_offset_len_pxs;
	
	my $pbs_end_offset_len_bp = $self->pbs_end - $self->start;
	my $pbs_end_offset_len_pxs = $pbs_end_offset_len_bp * $self->scale;
	my $pbs_end_pxs = $left + $pbs_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $pbs_start_pxs,
			  $top,
			  $pbs_end_pxs,
			  $bottom,
			  $self->pbs_bg_color,
			  $self->pbs_fg_color);
#			  $blue,
#			  $blue );
    }
    
    #-----------------------------+
    # POLYPURINE TRACT            |
    #-----------------------------+
    if ( ($self->rr_start) && ($self->rr_end) ) {

	my $rr_start_offset_len_bp = $self->rr_start - $self->start;
	my $rr_start_offset_len_pxs = $rr_start_offset_len_bp * $self->scale;
	my $rr_start_pxs = $left + $rr_start_offset_len_pxs;
	
	my $rr_end_offset_len_bp = $self->rr_end - $self->start;
	my $rr_end_offset_len_pxs = $rr_end_offset_len_bp * $self->scale;
	my $rr_end_pxs = $left + $rr_end_offset_len_pxs;

	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $rr_start_pxs,
			  $top,
			  $rr_end_pxs,
			  $bottom,
			  $self->rr_bg_color,
			  $self->rr_fg_color)
    }

    #-----------------------------+
    # ENVELOPE                    |
    #-----------------------------+
    if ( ($self->env_start) && ($self->env_end) ) {
	
	my $env_start_offset_len_bp = $self->env_start - $self->start;
	my $env_start_offset_len_pxs = $env_start_offset_len_bp * $self->scale;
	my $env_start_pxs = $left + $env_start_offset_len_pxs;
	
	my $env_end_offset_len_bp = $self->env_end - $self->start;
	my $env_end_offset_len_pxs = $env_end_offset_len_bp * $self->scale;
	my $env_end_pxs = $left + $env_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $env_start_pxs,
			  #$top+($height * 0.05 ),
			  $top,
			  $env_end_pxs,
			  $bottom,
			  #$bottom-($height * 0.05 ),
			  $self->env_bg_color,
			  $self->env_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "ENV";
	    my $feat_start_x = $env_start_pxs + 
		( ( ($env_end_pxs - $env_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);

	} # End of label feature gag
    }


    #-----------------------------+
    # GAG                         |
    #-----------------------------+
    if ( ($self->gag_start) && ($self->gag_end) ) {
	
	my $gag_start_offset_len_bp = $self->gag_start - $self->start;
	my $gag_start_offset_len_pxs = $gag_start_offset_len_bp * $self->scale;
	my $gag_start_pxs = $left + $gag_start_offset_len_pxs;
	
	my $gag_end_offset_len_bp = $self->gag_end - $self->start;
	my $gag_end_offset_len_pxs = $gag_end_offset_len_bp * $self->scale;
	my $gag_end_pxs = $left + $gag_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $gag_start_pxs,
			  $top+($height * 0.05 ),
			  #$top,
			  $gag_end_pxs,
			  #$bottom,
			  $bottom-($height * 0.05 ),
			  $self->gag_bg_color,
			  $self->gag_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "GAG";
	    my $feat_start_x = $gag_start_pxs + 
		( ( ($gag_end_pxs - $gag_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);

	} # End of label feature gag
    }

    #-----------------------------+
    # Zinc Finger CCHC            |
    #-----------------------------+
    if ( ($self->zf_start) && ($self->zf_end) ) {
	
	my $zf_start_offset_len_bp = $self->zf_start - $self->start;
	my $zf_start_offset_len_pxs = $zf_start_offset_len_bp * $self->scale;
	my $zf_start_pxs = $left + $zf_start_offset_len_pxs;
	
	my $zf_end_offset_len_bp = $self->zf_end - $self->start;
	my $zf_end_offset_len_pxs = $zf_end_offset_len_bp * $self->scale;
	my $zf_end_pxs = $left + $zf_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $zf_start_pxs,
			  $top+($height * 0.05 ),
			  #$top,
			  $zf_end_pxs,
			  #$bottom,
			  $bottom-($height * 0.05 ),
			  $self->zf_bg_color,
			  $self->zf_fg_color);
	
	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "Z";
	    my $feat_start_x = $zf_start_pxs + 
		( ( ($zf_end_pxs - $zf_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		    ) 
		  ) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature AP

    }

    #-----------------------------+
    # POL                         |
    #-----------------------------+
    # May want to offer option of drawing POL as cartoon and then
    # draw subcomponents on top of pol, these could be arrows within pol
    # These subcompoents will be drawn without the 5percent offset
    # and will be the same height as the 
    # In the same way it would be possible to annotate the u3/r/u5
    # of the LTRs
    if ( ($self->pol_start) && ($self->pol_end) ) {
	
	my $pol_start_offset_len_bp = $self->pol_start - $self->start;
	my $pol_start_offset_len_pxs = $pol_start_offset_len_bp * $self->scale;
	my $pol_start_pxs = $left + $pol_start_offset_len_pxs;
	
	my $pol_end_offset_len_bp = $self->pol_end - $self->start;
	my $pol_end_offset_len_pxs = $pol_end_offset_len_bp * $self->scale;
	my $pol_end_pxs = $left + $pol_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $pol_start_pxs,
			  $top+($height * 0.05 ),
			  #$top,
			  $pol_end_pxs,
			  #$bottom,
			  $bottom-($height * 0.05 ),
			  $self->pol_bg_color,
			  $self->pol_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "POL";
	    my $feat_start_x = $pol_start_pxs + 
		( ( ($pol_end_pxs - $pol_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature pol


    }

    #-----------------------------+
    # PROTEASE                    |
    #-----------------------------+
    if ( ($self->pro_start) && ($self->pro_end) ) {
	
	my $pro_start_offset_len_bp = $self->pro_start - $self->start;
	my $pro_start_offset_len_pxs = $pro_start_offset_len_bp * $self->scale;
	my $pro_start_pxs = $left + $pro_start_offset_len_pxs;
	
	my $pro_end_offset_len_bp = $self->pro_end - $self->start;
	my $pro_end_offset_len_pxs = $pro_end_offset_len_bp * $self->scale;
	my $pro_end_pxs = $left + $pro_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $pro_start_pxs,
			  $top,
			  $pro_end_pxs,
			  $bottom,
			  $self->pro_bg_color,
			  $self->pro_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "AP";
	    my $feat_start_x = $pro_start_pxs + 
		( ( ($pro_end_pxs - $pro_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);

	} # End of label feature AP

    }

    #-----------------------------+
    # INTEGRASE                   |
    #-----------------------------+
    if ( ($self->int_start) && ($self->int_end) ) {
	
	my $int_start_offset_len_bp = $self->int_start - $self->start;
	my $int_start_offset_len_pxs = $int_start_offset_len_bp * $self->scale;
	my $int_start_pxs = $left + $int_start_offset_len_pxs;
	
	my $int_end_offset_len_bp = $self->int_end - $self->start;
	my $int_end_offset_len_pxs = $int_end_offset_len_bp * $self->scale;
	my $int_end_pxs = $left + $int_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $int_start_pxs,
			  $top,
			  $int_end_pxs,
			  $bottom,
			  $self->int_bg_color,
			  $self->int_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "INT";
	    my $feat_start_x = $int_start_pxs + 
		( ( ($int_end_pxs - $int_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature integrase


    }
    
    #-----------------------------+
    # REVERSE TRANSCRIPTASE       |
    #-----------------------------+
    if ( ($self->rt_start) && ($self->rt_end) ) {
	
	my $rt_start_offset_len_bp = $self->rt_start - $self->start;
	my $rt_start_offset_len_pxs = $rt_start_offset_len_bp * $self->scale;
	my $rt_start_pxs = $left + $rt_start_offset_len_pxs;
	
	my $rt_end_offset_len_bp = $self->rt_end - $self->start;
	my $rt_end_offset_len_pxs = $rt_end_offset_len_bp * $self->scale;
	my $rt_end_pxs = $left + $rt_end_offset_len_pxs;
	
	# Using non black color for outline to make stand out
	$self->filled_box($gd,
			  $rt_start_pxs,
			  $top,
			  $rt_end_pxs,
			  $bottom,
			  $self->rt_bg_color,
			  $self->rt_fg_color);
	
	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "RT";
	    #my $feat_val = "T";
	    my $feat_start_x = $rt_start_pxs + 
		( ( ($rt_end_pxs - $rt_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature rt


    }
    
    #-----------------------------+
    # RNASEH                      |
    #-----------------------------+
    if ( ($self->rh_start) && ($self->rh_end) ) {
	
	my $rh_start_offset_len_bp = $self->rh_start - $self->start;
	my $rh_start_offset_len_pxs = $rh_start_offset_len_bp * $self->scale;
	my $rh_start_pxs = $left + $rh_start_offset_len_pxs;
	
	my $rh_end_offset_len_bp = $self->rh_end - $self->start;
	my $rh_end_offset_len_pxs = $rh_end_offset_len_bp * $self->scale;
	my $rh_end_pxs = $left + $rh_end_offset_len_pxs;
	
	$self->filled_box($gd,
			  $rh_start_pxs,
			  $top,
			  $rh_end_pxs,
			  $bottom,
			  $self->rh_bg_color,
			  $self->rh_fg_color);

	if ($self->label_feat) {
	    
	    # Change the following from seq context color to
	    # feature font color 
	    # This will not work if start > end
	    my $feat_val = "RH";
	    my $feat_start_x = $rh_start_pxs + 
		( ( ($rh_end_pxs - $rh_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		  ) 
		) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature integrase

    }

    #-----------------------------+
    # CHRromatin Organiation      |
    # MOdifier Domain             |
    #-----------------------------+
    if ( ($self->chromo_start) && ($self->chromo_end) ) {
	
	my $chromo_start_offset_len_bp = $self->chromo_start - $self->start;
	my $chromo_start_offset_len_pxs = $chromo_start_offset_len_bp * $self->scale;
	my $chromo_start_pxs = $left + $chromo_start_offset_len_pxs;
	
	my $chromo_end_offset_len_bp = $self->chromo_end - $self->start;
	my $chromo_end_offset_len_pxs = $chromo_end_offset_len_bp * $self->scale;
	my $chromo_end_pxs = $left + $chromo_end_offset_len_pxs;
	
	$self->filled_box($gd,
			  $chromo_start_pxs,
			  $top,
			  $chromo_end_pxs,
			  $bottom,
			  $self->chromo_bg_color,
			  $self->chromo_fg_color);

	if ($self->label_feat) {

	    #my $feat_val = "CHROMO";
	    my $feat_val = "C";
	    my $feat_start_x = $chromo_start_pxs + 
		( ( ($chromo_end_pxs - $chromo_start_pxs)/2 ) - 
		  ( (($feat_font->width * length($feat_val))/2) 
		    ) 
		  ) ;
	    my $char_height = $feat_font->height;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);
	    
	    $gd->string($feat_font, 
			$feat_start_x, 
			$feat_start_y, 
			$feat_val,
			$self->seq_con_color);
	    
	} # End of label feature integrase
	
    }

    #-----------------------------+
    # TARGET SITE DUPLICATION     |
    #-----------------------------+
    # Default length is zero, no TSDs are annotated
    # otherwise draw these just outside
    if ($self->tsd_len > 0) {

	# DRAW LEFT TSD
	$self->filled_box($gd,
			  $left-($self->tsd_len*$self->scale),
			  $top,
			  $left,
			  $bottom,
			  $self->tsd_bg_color,
			  $self->tsd_fg_color
	    );
	
	# DRAW RIGHT TSD
	$self->filled_box($gd,
			  $right,
			  $top,
			  $right+($self->tsd_len*$self->scale),
			  $bottom,
			  $self->tsd_bg_color,
			  $self->tsd_fg_color );
    }

    #-----------------------------+
    # 5' SEQUENCE CONTEXT         |
    #-----------------------------+
    if ($self->seq_con_ltr5) {

	my $char_width = $con_font->width;
	my $char_height = $con_font->height;
	my $seq_width = length($self->seq_con_ltr5) * $char_width;

	my $seq_val = $self->seq_con_ltr5;
	my $con_start_y =  $top + ($height * 0.5) - ($char_height/2);

	my $con_start_x = $left - $seq_width - 1;

	if ($self->tsd_len) {
	    my $offset = $self->tsd_len*$self->scale;
	    $offset = roundup($offset);
	    $con_start_x = $con_start_x - $offset;
	}
	
	$gd->string($con_font, 
		    $con_start_x, 
		    $con_start_y, 
		    $seq_val, 
		    $self->seq_con_color);
    }

    #-----------------------------+
    # 3' SEQUENCE CONTEXT         |
    #-----------------------------+
    if ($self->seq_con_ltr3) {
	
	my $char_width = $con_font->width;
	my $char_height = $con_font->height;
	my $seq_width = length($self->seq_con_ltr3) * $char_width;

	# Need to switch the following to not use pad_right since pad_right
	# will need to include the text string
	my $con_start_y =  $top + ($height * 0.5) - ($char_height/2);
	my $seq_val = $self->seq_con_ltr3;

	my $con_start_x = $right + 1;
	
	if ($self->tsd_len) {
	    my $offset = $self->tsd_len*$self->scale;
	    $offset = roundup($offset);
	    $con_start_x = $con_start_x + $offset;
	}

	$gd->string($con_font, 
		    $con_start_x, 
		    $con_start_y, 
		    $seq_val, 
		    $self->seq_con_color);
    }

}
 
1;

#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 03/13/2008 -- 03/14/2008
#  - Main object started, includes
#    - LTR by start stopewith triangle
#    - LTR by Length
#    - TSD
#    - GAG
#    - Pol
#    - Pro
#    - INT
#    - RH
#    - RT
#    - PBS
#    - PPT (using RR to follow song complience)
#  - Single object that only would be interpreted correctly
#    on the plus strand for a full length LTR retro
#
# 03/14/2008
# - Added boolean for drawing LTR arrows, default is true 1
# - Added half circles for the tg/ca based at either end of 
#   the LTR
# - Adding color features as variables
#    
# 03/16/2008
# - Added color features as variables for 
#    Pol, Gag, RR, PBS, INT, RH, PRO, LTR Arrows
# - Started genomic context 
#
# 03/17/2008
# - Continuing work on genomic context
#    This buffers sequences by one pixel on either side to 
#    avoid overlap with drawn objects, 
# - Modifed left and right padding to take into account all
#   sequence features that can exist outside of the drawn 
#   object. This includes a 1 pixel offset at all times.
# - Added tgca foreground and background colors as variables
# - Cleaned up code, removed STDERR print statments
# - Added span bgcolor and fg color
# - Added Lablels for Gag, Pol, AP, INT, RT, RH
# 03/26/2008
# - Added chromo domain with
#    -chromo_start,chromo_end,chromo_bg_color, chromo_fg_color
#    -This domain is about 60 aa so will not be labeled for now
#
# 03/28/2008
# - Adding env domain
#
# 04/18/2008
# - Working on dealing with value 'N' being passed as a coordinate
#   I will treat this as a null/N value and not draw the feature
#   by returned 0/false when fetching feature locations
# - Added labels to zinc finger cchc, env, and chromo domains 
