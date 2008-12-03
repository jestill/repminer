#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# extract_blast.pl - Simple extract of m8 BLAST to text     |
#                                                           |
#-----------------------------------------------------------+
#
# Extract blast output from the table in the form of
# # mysql> create table tbl_ltr_ltr (
#                qry INT,
#                sub INT,
#                bit INT,
#                INDEX (qry),
#                INDEX (sub));
#
#
#

# INCLUDES
use strict;
use DBI;
use Getopt::Long;

# VARIABLE SCOPE AND DEFAULTS
my $usrname = $ENV{DBI_USER};
my $pass=$ENV{DBI_PASSWORD};
my $dsn=$ENV{MAIZE_DSN};      # DSN FOR THE MAIZE DB

my $outfile;
my $prefout;                  # not required, Preference vector outfile

my $blast_tbl="tbl_ltr_ltr";

# BOOLEANS
my $verbose = 0;

# GET VARS FROM COMMAND LINE
my $ok = GetOptions ("d|dsn=s"      => \$dsn,
		     "u|dbuser=s"   => \$usrname,
		     "p|dbpass=s"   => \$pass,
		     "o|outfile=s"  => \$outfile,
		     # Additional options
		     "r|pref=s"     => \$prefout,
		     "verbose"      => \$verbose,
		     "t|table=s"    => \$blast_tbl,
		     );

print "Starting $0 ..\n" if $verbose;

# GET USER PASSWORD
unless ($pass) {
    print "\nEnter password for the user $usrname\n";
    system('stty', '-echo') == 0 or die "can't turn off echo: $?";
    $pass = <STDIN>;
    system('stty', 'echo') == 0 or die "can't turn on echo: $?";
    chomp $pass;
}


#-----------------------------+
# PREFERENCE
#-----------------------------+
if ($prefout) {
    open (PREFOUT, ">$prefout") ||
	die "Can not open preference outfile $prefout\n"
}


# Only supporting MYSQL at the moment
my $dbh = &connect_to_db($dsn, $usrname, $pass);

# THIS ASSUMES THAT i to j is without gaps
# This is probably and oversimplification
#

# GET THE MAX VALUE FOR THE qry SEQUENCE ID
my $sql = "SELECT MAX(qry) FROM $blast_tbl";
my $cur = $dbh->prepare($sql);
$cur->execute();
my @row=$cur->fetchrow;
my $max_qry=$row[0];
$cur->finish();

my $max_val = $max_qry+1;

for (my $i=1;$i<$max_val;$i++ ) {

    # Get median score value if preference outfile given
    if ($prefout) {
	# Fetch the array of hits for each qry ID
	my $qry = "SELECT bit FROM $blast_tbl".
	    " WHERE qry=$i";
	my $res = $dbh->prepare($qry);
	$res->execute();
	my @vals = @{$res->fetchall_arrayref([0])};
	my $num_vals = @vals;
	
	# Get median bit score from this array
	#my $median = get_median(@vals);
	# Get median as integer value
	my $median = int (get_median(@vals));
	
	# Print to median val output file
	print PREFOUT "$median\n";
	print STDERR "$i\t$median\n" if $verbose;
    }


} # END OF j

# FINISHED
# Shut down the dbh connection
$dbh->disconnect();

print "$0 finished\n" if $verbose;

close PREFOUT if ($prefout);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub connect_to_db {
    my ($cstr) = @_;
    return connect_to_mysql(@_) if $cstr =~ /:mysql:/i;
    return connect_to_pg(@_) if $cstr =~ /:pg:/i;
    die "can't understand driver in connection string: $cstr\n";
}

sub connect_to_pg {

    my ($cstr, $user, $pass) = @_;

        my $dbh = DBI->connect($cstr, $user, $pass,
                               {PrintError => 0,
                                RaiseError => 1,
                                AutoCommit => 0});
    $dbh || &error("DBI connect failed : ",$dbh->errstr);

    return($dbh);
} # End of ConnectToPG subfunction   


sub connect_to_mysql {

    my ($cstr, $user, $pass) = @_;

    my $dbh = DBI->connect($cstr,
                           $user,
                           $pass,
                           {PrintError => 0,
                            RaiseError => 1,
                            AutoCommit => 0});

    $dbh || &error("DBI connect failed : ",$dbh->errstr);

    return($dbh);
}

sub get_median {
    # from http://www.sourcesnip.com/2006/11/13/median-value-in-perl/
    my $rpole = shift;
    my @pole = @$rpole;

    my $ret;

    sort(@pole);

    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }

    return $ret;
}
