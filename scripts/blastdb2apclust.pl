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

my $blast_tbl="tbl_ltr_ltr";

# BOOLEANS
my $verbose = 0;

# GET VARS FROM COMMAND LINE
my $ok = GetOptions ("d|dsn=s"      => \$dsn,
		     "u|dbuser=s"   => \$usrname,
		     "p|dbpass=s"   => \$pass,
		     "o|outfile=s"  => \$outfile,
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


open (APOUT, ">$outfile") ||
    die "Can not open ap outfile.\n";

# PRINT THE MAX VALUE
print STDERR "MAX VALUE IS: $max_qry\n" if $verbose;


#
# SELECT BIT SCORES FOR ALL I and J 
# Could use BIT SCORE of I=J TO BE PREFERENCE VALUE
# THIS WOULD GIVE PREFERENCE TO THE LONGER SEQUENCES

# THE FOLLOWING IS A TEMP MAX VAL FOR DEBUG
#my $max_val = 5;

my $max_val = $max_qry+1;

for (my $i=1;$i<$max_val;$i++ ) {
    for (my $j=1;$j<$max_val;$j++) {
	
	print STDERR "processing query $i\n" if $verbose;

	my $qry = "SELECT * FROM $blast_tbl".
	    " WHERE qry=$i AND sub=$j";
	my $res = $dbh->prepare($qry);
	$res->execute();

	#////////////////////////////////
	# THIS IS THE INCORRECT WAY TO EXTRACT A SET OF INFORMATION
	my @vals=$res->fetchrow;
	my $ans = $vals[2] || "-Inf";
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	print APOUT "$i\t$j\t$ans\n";

	print STDERR "$i\t$j\t$ans\n" if $verbose;

	$res->finish();
	    
    } # END OF i

    # TEMP EXIT FOR DEBUG
    #if ($i==5) {$dbh->disconnect();exit;}
} # END OF j

# FINISHED
# Shut down the dbh connection
$dbh->disconnect();

print "$0 finished\n" if $verbose;

close APOUT;

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
