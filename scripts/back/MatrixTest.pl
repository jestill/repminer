#!/usr/bin/perl -w
# 
# 

my $MaxRow = 5;
my $MaxCol = 10;
my @FeatMatrix;

for (my $row=1; $row<=$MaxRow; $row++)
{
    for (my $col=1; $col<=$MaxCol; $col++)
    {
	$FeatMatrix[$row][$col] = $row.":".$col;
	print $FeatMatrix[$row][$col]."\t";
    }

    print "\n"
}


exit;


