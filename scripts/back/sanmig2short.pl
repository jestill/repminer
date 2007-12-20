#!/usr/bin/perl -w
#  AUTHOR: James C. Estill
# STARTED: 03/20/2007
# UPDATED: 03/20/2007

# Quick and dirty script to get the "short name"
# for the SanMiguel TE hits

my $InFile = shift or die;
my $OutFile = shift or die;

open (IN, "<$InFile") ||
    die "Can not open input file:\n$InFile\n";
open (OUT, ">$OutFile") ||
    die "Can not open output file:\n$OutFile";
print OUT "SrcShortName\n";


while (<IN>)
{
    my @EqSplit = split (/\=/);
    my @UndSplit = split (/\_/, $EqSplit[1]);
    print OUT $EqSplit[0]."=".$UndSplit[0]."\n";
}

exit;
