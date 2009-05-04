#!/usr/bin/perl -w

# Test to get the log value of a number


my $num = shift or die "Need a value";

# When num =1 the neg log value returns -0
# When num=0 the neg log returns, can not take log of zero
# log 0 is undefined, so will need to figure out how to score these ...?
# Consider using the  Math::Complex library
#http://search.cpan.org/~jhi/Math-Complex-1.54/lib/Math/Complex.pm#ERRORS_DUE_TO_DIVISION_BY_ZERO_OR_LOGARITHM_OF_ZERO
#my $num=0;

if ($num ==0) {
#    $num = 1e-1000;
    -realmax;
}

my $log = log($num);
my $log_ten = log_10($num);
my $log_neg = neg_log($num);
my $neg_log_ten = -log_10($num);

# These do not work for e value = 0

print "ANSWERS:\n";
print "LOG:\t $log\n";
print "TEN:\t $log_ten\n";
print "NEG:\t $log_neg\n";
print "NEG10:\t $neg_log_ten\n";

exit;

# SUBFUNCTIONS

# LOG 10
sub log_10 {
    my $n = shift;
    return log($n)/log(10);
}

# NEGATIVE LOG
sub neg_log {
    my $n = shift;

#    if ($n == 0) {
#	return "null";
#    }
#    else {
	return -log($n);
#    }

}

#sub 
