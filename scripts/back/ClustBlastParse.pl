#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# BLAST PARSER ::                                           |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/27/2006                                       |
# UPDATED: 09/22/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Parse a BLAST file, and report some contents based on    |
#  criteria set by the user.                                | 
#                                                           |
# USAGE:                                                    |
#  BLParser.pl                                              |
#
# DEPENDENCIES:                                             |
#  -BioPerl
#                                                           |
#-----------------------------------------------------------+

#print "The program has started\n";

#-----------------------------+ 
# INCLUDES                    |
#-----------------------------+
use Bio::SearchIO;            # Needed to parse BLAST output
use DBI();                    # Needed for database connections

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $first = "true";
#my $first = "false";


#-----------------------------+
# BLASTN VALUES               |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "AllBlastn_20060925.txt";
#
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_os_rep.blo";
#    "asgr_con_zm_rep.blo";
#    "asgr_con_Wessler.blo";
#    "asgr_con_TREP_8.blo";
#    "asgr_con_SanMiguel.blo";
#    "asgr_con_RB_pln.blo";
#    "asgr_con_RB_ory.blo";
#    "asgr_con_RB_ath.blo";
#    "asgr_con_gram_rep.blo";

#-----------------------------+
# TBLASTX VALUES              |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "AlltBlastx_20060925.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/blastx/".
#    "asgr_con_os_rep.blx";
#    "asgr_con_zm_rep.blx";
#    "asgr_con_Wessler.blx";
#    "asgr_con_TREP_8.blx";
#    "asgr_con_SanMiguel.blx";
#    "asgr_con_RB_pln.blx";
#    "asgr_con_RB_ory.blx";
#    "asgr_con_RB_ath.blx";
#    "asgr_con_gram_rep.blx";

#-----------------------------+
# TBLASTX MERGED              |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "AlltBlastx_20060925_merge.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/blastx/".
#    "asgr_con_merge.blx";

#-----------------------------+
# BLASTN MERGED               |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "AllBlastn_20060925_merge.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_merge.blo";

#-----------------------------+
# TREP 9 TOTAL                |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_TREP9_total_Blastn.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_TREP9_total.blo";
#    "asgr_con_merge.blo";

#-----------------------------+
# NEW SANMIGUEL RUN           |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_SanMiguel_200610.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_SanMiguel_200610.blo";
#    "asgr_con_merge.blo";

#-----------------------------+
# SAMI DATABASE               |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_SAMI_1.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_SAMI_1.blo";

#-----------------------------+
# PTREP PARSER                |
# 10/16/2006                  |
#-----------------------------+
#my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_PTRE9.txt";
#my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
#    "asgr_con_PTREP_9.blo";

#-----------------------------+
# MIPS REdat                  |
#-----------------------------+
my $SummaryOut = "/home/jestill/projects/asgr/200609/asgr_con/".
    "asgr_con_mips_REdata_4_3_TESTData.txt";
my $InFilePath = "/home/jestill/projects/asgr/200609/asgr_con/".
    "asgr_con_mips_REdat_4_3.blo";

my $MinQryLen = "50";
my $MinHitLen = "50";  # The minimum hit length to consider
my $MinScore = "50";
my $MaxE = "1.0e-06";
my $BlastResult;      # The BLAST Result 
my $BlastHit;         # Individual BLAST hit
my $HSP;              # BLAST HSP Result

my $BlastDB;          # The blast database that was queried
my $HitId;            # Unique ID of the BLAST hit
my $SpHInfo;          # Information related to the BLAST hit
my @HitDesc;          # Initially split the hit description
my @HitInfo;          # Array to hold the split of the hit info
my ($Class, $Subclass, $Superfamily, $Family, $Name ); 
                      # Information related to the repeat using Cedric's
                      # nomenclature from Nature review paper
                      # Class like Class 1 or Class 2
                      # Subclass like LTR/NonLTR/DNA Transposons
                      # Superfamily like copia, gypsy, lines, since
my $ResultCount = 0;
my $HitLength;

my $HSPMinStart;
my $HSPMaxEnd;


#-----------------------------+
# FILE I/O CONNECTIONS        |
#-----------------------------+

# OPEN THE BLAST REPROT FOR INTPUT
my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
				      '-file'   => $InFilePath, 
				      '-signif' => $MaxE, 
				      '-min_query_len' => $MinQryLen,
				      '-score' => $MinScore ) ||
    die "Could not open BLAST input file:\n$InFilePath.\n";

# OPEN THE OUTPUT FILE FOR OVERWRITE
#open (OUT, ">".$SummaryOut) ||

#OPEN THE FILE TO APPEND
open (OUT, ">>".$SummaryOut) ||
    die "Could not open out file\n$SummaryOut\n";

if ($first =~ "true")
{
    # COMMENT OUT FOR FIRST
    print OUT "QUERY\t";          # Name of the query seq
    print OUT "ALG\t";            # BLAST algorithm
    print OUT "DB\t";             # Database blasted against
    print OUT "LENGTH\t";         # Lenth of the hit
    print OUT "START\t";          # Minimum start of HSPs
    print OUT "END\t";            # Maximum end of HSPS
    print OUT "HSP-SPAN\t";       # Span of the HSP
    print OUT "SIG\t";            # Significance of the hit
    print OUT "BIT\t";            # Bit Score of the hit
    print OUT "NumHSP\t";         # Number of HSPS
    print OUT "NAME\t";           # Name of the repeat hit
    print OUT "CLASS\t";          # Class of the element
    print OUT "\n";
}

#-----------------------------------------------------------+
# DATABASE CONNECTION                                       |
#-----------------------------------------------------------+

#-----------------------------+
# DB VARIABLES                |
#-----------------------------+
my $DbName = "dbASGR";             # Database name to connect to
my $tblSrcTable = "tblSeqData";    # Table with src sequence data
my $tblDatCat = "tblDatCat";       # Table with information about data
                                   # categories. Can include color
                                   # and symbol for drawing.
my $DbUserName = "jestill";        # User name for the db connection

#-----------------------------+
# GET DB USER PASSWORD        |
#-----------------------------+
print "\nPassword for $DbUserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
$DbUserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $DbUserPassword;

#-----------------------------+
# CONNECT TO THE REPEAT DB    |
#-----------------------------+
$RepDbName = "dbRep";
my $RepDB = DBI->connect("DBI:mysql:database=$RepDbName;host=localhost",
			 $DbUserName, $DbUserPassword,
			 {'RaiseError' => 1});

#-----------------------------------------------------------+
# WORK WITH THE BLAST RESULT                                |
#-----------------------------------------------------------+
while ($BlastResult = $BlastReport->next_result())
{
    
    $ResultCount++;
    print $ResultCount."\n";

    #-----------------------------+
    # TEMP EXIT TO CHECK OUT HOW  |
    # THE PROGRAM IS WORKING      |
    #-----------------------------+
    #if ($ResultCount == 30) { exit; }

    while ( $BlastHit = $BlastResult->next_hit())
    {

	$HitLength = $BlastHit->length();

	#-----------------------------------------------------------+
	# EXTRACT REPEAT CLASS AND NAME FROM THE REPEAT DATABASES   |
	# USING FUNCTIONS SPECIFIC TO THE DATABASES THAT ARE        |
	# USED                                                      |
	#-----------------------------------------------------------+
	$BlastDB = $BlastResult->database_name;
	
	#-----------------------------+
	# TREP REPEAT DATABASE        |
	#-----------------------------+
	if ($BlastDB =~ 'TREP_8')
	{
	    $HitCount++;	    
	    my $HitId = $BlastHit->accession()  || "UnkAcc";
	    $RepClass = &GetTREPClass($HitId);
	    $Name = &GetTREPName($HitId);    
	}
	
	#-----------------------------+
	# REPBASE PLANTS              |
	#-----------------------------+
	elsif ($BlastDB =~ 'RB_pln')
	{
	    
	    $HitCount++;
	    $Name = $BlastHit->name();
	    my @SpName = (split /\#/, $Name);
	    my $CatSearch = trim($SpName[1] || "NONE");
	    $RepClass = &GetRBClass($CatSearch); 
	    
	}

	#-----------------------------+
	# REPBASE ARABIDOPSIS         |
	#-----------------------------+
	elsif ($BlastDB =~ 'RB_ath')
	{
	    
	    $HitCount++;
	    $Name = $BlastHit->name();
	    my @SpName = (split /\#/, $Name);
	    my $CatSearch = trim($SpName[1] || "NONE");
	    $RepClass = &GetRBClass($CatSearch); 
	    
	}


	#-----------------------------+
	# REPBASE ARABIDOPSIS         |
	#-----------------------------+
	elsif ($BlastDB =~ 'RB_ory')
	{
	    
	    $HitCount++;
	    $Name = $BlastHit->name();
	    my @SpName = (split /\#/, $Name);
	    my $CatSearch = trim($SpName[1] || "NONE");
	    $RepClass = &GetRBClass($CatSearch); 
	    
	}

	#-----------------------------+
	# TIGR ORYZA REPEAT DATABASE  |
	#-----------------------------+
	elsif ($BlastDB =~ 'os_rep')
	{
	    $HitCount++;
	    $Name = $BlastHit->name();
	    $RepClass = &GetTIGRClass($Name);
	}
	
	#-----------------------------+
	# WESSLER LAB REPEAT DATABASE |
	#-----------------------------+
	elsif ($BlastDB =~ 'Wessler')
	{
	    $HitCount++;
	    my $WesName = $BlastHit->name();
	    my @SpName = (split /\#/, $WesName);
	    $Name = $SpName[0];	    
	    my $WesClass = $SpName[1];
	    $RepClass = &GetWesClass($WesClass);	    
	}
	
	#-----------------------------+
	# SAN MIGUEL REPEAT DATABASE  |
	#-----------------------------+
	elsif ($BlastDB =~ 'SanMiguel')
	{
	    $HitCount++;
	    $Name = $BlastHit->name();
	    $RepClass = "UNK-SanMiguel";
	}
	
	#-----------------------------+
	# TIGR ZEA MAYS REPEAT        |
	# DATABASE                    |
	#-----------------------------+
	elsif ($BlastDB =~ 'zm_rep')
	{
	    $HitCount++;
	    $Name = $BlastHit->name();
	    $RepClass = &GetTIGRClass($Name); 
	}
	
	#-----------------------------+
	# TIGR GRAMINEAE REPEAT       |
	# DATABASE                    |
	#-----------------------------+
	elsif ($BlastDB =~ 'gram_rep')
	{
	    $HitCount++;
	    $Name = $BlastHit->name();
	    #$RepClass = "UNK-ZeaRepeat";
	    $RepClass = &GetTIGRClass($Name); 
	}

	

	#-----------------------------------------------------------+
	# 
	#
	# AT WORK
	# GETTING START AND STOP OF THE SET OF HSPS USING
	#
	#
	#-----------------------------------------------------------+
	# Reset HSP values to null
	$HSPMinStart = "99999999"; # Fake value
	$HSPMaxEnd = "0";          # Fake value	
	my $HitFracIdent = "";     # Fraction of identity for the HIT

	while ( $HSP = $BlastHit->next_hsp())
	{
	    #-----------------------------+
	    # GET THE START AND END
	    # COORDINATES OF THE HSPS    
	    #-----------------------------+
	    my $HSPStart = $HSP->start;
	    my $HSPEnd = $HSP->end;


	    #-----------------------------+
	    # SHOW SOME ADDITIONAL        |
	    # HSP INFORMATION             |
	    # 10/16/2006                  |
	    #-----------------------------+
	    # Working to incorporate this into the data used to 
	    # detemine 'signicant' hits
#	    $HitFracIdent = ($HSP->frac_identical('total') * 
#			     $HSP->length) + 
	    print "\tIDN:\t".$HSP->frac_identical('total')."\n";
	    print "\tLEN:\t".$HSP->length('total')."\n";
	    print "\tGAP:\t".$HSP->gaps."\n";



	    #-----------------------------+
	    # GET THE MAX AND MIN VALUES  |
	    #-----------------------------+

	    if ($HSPStart < $HSPEnd)
	    {
		#-----------------------------+
		# HSPStart IS THE MINIMUM VAL |
		#-----------------------------+

		# GET THE MIN VALUE
		if ( $HSPStart < $HSPMinStart){
		    $HSPMinStart = $HSPStart;}

		# GET THE MAX VALUE
		if ( $HSPEnd > $HSPMaxEnd){
		    $HSPMaxEnd = $HSPEnd;}
		
	    }else{
		#-----------------------------+
		# HSPStart IS THE MAXIMUM VAL |
		#-----------------------------+

		# GET THE MIN VALUE
		if ( $HSPEnd < $HSPMinStart){
		    $HSPMinStart = $HSPEnd;}

		# GET THE MAX VALUE
		if ( $HSPStart > $HSPMaxEnd){
		    $HSPMaxEnd = $HSPStart;}


	    }
	    

	}

	#-----------------------------------------------------------+
	# PRINT OUTPUT TO THE OUTPUT TAB DELIM TEXT FILE            |
	#-----------------------------------------------------------+
	my $QryName = $BlastResult->query_name || "NULL";
	my $HitName = $BlastHit->name()        || "NULL";
	my $BitScore = $BlastHit->bits()       || "NULL";
	my $HitSig = $BlastHit->significance() || "NULL";
	my $HitHsps = $BlastHit->num_hsps      || "NULL";
	my $BlAlg = $BlastHit->algorithm()     || "NULL";
	my $HitLen = $BlastHit->length()       || "NULL";
	my $HSPLen = $HSPMaxEnd - $HSPMinStart;

	# DEAL WITH MESSED UP VALUES
	if ($HSPLen ==  -99999999)
	{
	    $HSPMinStart = "NULL";
	}

	print OUT $QryName."\t";           # Name of the query seq
	print OUT $BlAlg."\t";             # BLAST algorithm
	print OUT $BlastDB."\t";
	print OUT $HitLen."\t";            # Length of the hit
	print OUT $HSPMinStart."\t";       # Min position of all HSPs
	print OUT $HSPMaxEnd."\t";         # Max position of all HSPs
	print OUT $HSPLen."\t";            # Len of the span of the HSP
	print OUT $HitSig."\t";            # Significance of the hit
	print OUT $BitScore."\t";          # Bit Score of the hit
	print OUT $HitHsps."\t";           # Number of HSPS
	print OUT $HitName."\t";           # Name of the repeat hit
	print OUT $RepClass."\t";          # Class of the element
	print OUT "\n";

    } # End of while BlastHit->next_hit
    
} # End of while BlastReport->next_result

close OUT;

#-----------------------------+
# WHEN FINISHED SOUND THE     |
# BELL                        |
#-----------------------------+
print "\a\a";

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub GetWesClass
{
#-----------------------------+
# GETS THE TERMINAL           |
# CLASSIFICATION NAME FOR HITS|
# TO THE WESSLER REPEAT       |
# DATABASE FROM MAGI SITE     |
#-----------------------------+
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    my %WesClass;

    %WesClass = (
		 # RETROTRANSPOSONS
		 'LTR/copia' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/Copia' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/COPIA' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/copia?' => 'Ty1-copia LTR-Retrotransposon',
		 'LTR/gypsy' => 'Ty3-gypsy LTR-Retrotransposon',
		 'LTR/Gypsy' => 'Ty3-gypsy LTR-Retrotransposon',
		 'LTR/?' => 'LTR-Retrotransposon',
		 'LTR/unknown' => 'LTR-Retrotransposon',
		 'LTR/Barley' => 'LTR-Retrotransposon',
		 'LTR/Maize' => 'LTR-Retrotransposon',
		 'SINE/12-15bp?' => 'SINE Retrotransposon',
		 'LINE' => 'LINE Retrotransposon',
		 'LTR/Tto1_NTlike' => 'LINE Retrotransposon',
		 'LTR/Zea' => 'LTR-Retrotransposon',
		 # TRANSPOSONS
		 'DNA/8bp' => 'Transposon',
		 'DNA/hAT' => 'hAT Transposon',
		 'DNA/CACTA' => 'CACTA, En/Spm Transposon',
		 'DNA/spm' => 'CACTA, En/Spm Transposon',
		 'DNA/CACTG' => 'CACTA, En/Spm Transposon',
		 'DNA/TAA' => 'ping/pong/SNOOPY Transposon',
		 # FOLDBACK ELEMENTS
		 "MITE/NNN" => 'MITE',
		 'MITE/TA' => 'MITE',
		 'MITE/TA(stw)' => 'MITE',
		 # OTHER
		 'vector' => 'Vector',
		 'Vector/Bao' => 'Vector',
		 'unknown' => 'Unknown',
		 'Novel/?' => 'Unknown',
		 'unknown/Bao' => 'Unknown'
		 );

    $cl = $WesClass{$qry} || "UNK";
    return $cl;   
}

sub GetTREPClass
{
#-----------------------------+
# GET THE CLASSIFICATION OF   |
# AN ELEMENT GIVEN THE TREP   |
# ACCESSION ID                |
#-----------------------------+
# This will need to use tblTREPMeta table in the database
    
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned

    my $QryTbl = "tblTREP";

    my $SelQry = "SELECT j_class FROM ".$QryTbl.
	" WHERE accession = '".$qry."'";

    $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    @row=$cur->fetchrow;
    $cl = $row[0] || "TREP CLASS UNK";
    $cur->finish();

    return $cl;

}


sub GetTREPName
{
    
    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    
    my $QryTbl = "tblTREP";

    my $SelQry = "SELECT name FROM ".$QryTbl.
	" WHERE accession = '".$qry."'";

    $cur = $RepDB->prepare($SelQry);
    $cur->execute();
    @row=$cur->fetchrow;
    $cl = $row[0] || "TREP NAME UNK";
    $cur->finish();

    #$cl = "Currently Unknown";
    return $cl;

}

sub GetTIGRClass
{

    #-----------------------------+
    # PARSE TIGR REPEAT CODES     |
    # Parse the codes used by TIGR|
    # for Repeat names. Returns   |
    # SuperClass-Class-SubClass   |
    #-----------------------------+


    my $qry = $_[0];        # The qry code sent to the subfunction
    my $cl;                 # The classification that is returned
    my %TIGRRep;            # Hash to translate from code to class
    my $code = substr ($qry,4,7);
    
    #-----------------------------+
    # HASH TO TRANSLATE FROM TIGR |
    # CODES TO REPEAT CLASSES     |
    #-----------------------------+
    %TIGRRep = (
		# RETROTRANSPOSONS
		TERT001 => 'Ty1-copia LTR-Retrotransposon',
		TERT002 => 'Ty3-gypsy LTR-Retrotransposon',
		TERT003 => 'LINE Retrotransposon',
		TERT004 => 'SINE Retrotransposon',
		TERTOOT => 'Unclassified Retrotransposon',
		# DNA TRANSPOSONS
		TETN001 => 'Ac/Ds Transposon',
		TETN002 => 'CACTA, En/Spm Transposon',
		TETN003 => 'Mutator (MULE) Transposon',
		TETN004 => 'Mariner (MLE) Transposon',
		TETN005 => 'ping/pong/SNOOPY Transposon',
		TETNOOT => 'Unclassified Transposon',
		# MITES
		TEMT001 => 'Tourist MITE',
		TEMT002 => 'Stowaway MITE',
		TEMT003 => 'Crackle MITE',
		TEMT004 => 'Explorer MITE',
		TEMT005 => 'Gaijin/Gaigin MITE',
		TEMT006 => 'Castaway MITE',
		TEMT007 => 'Snap MITE',
		TEMT008 => 'Amy/LTP MITE',
		TEMT009 => 'Ditto MITE',
		TEMT010 => 'Wanderer MITE',
		TEMT011 => 'p-SINE1 MITE',
		TEMT012 => 'Pop MITE',
		TEMT013 => 'Krispie MITE',
		TEMT014 => 'Snabo MITE',
		TEMT015 => 'MITE-adh, typeA MITE',
		TEMT016 => 'MITE-adh, typeB MITE',
		TEMT017 => 'MITE-adh, typeD MITE',
		TEMT018 => 'MITE-adh, typeG MITE',
		TEMT019 => 'MITE-adh, typeH MITE',
		TEMT020 => 'MITE-adh, typeI MITE',
		TEMT021 => 'MITE-adh, typeJ MITE',
		TEMT022 => 'MITE-adh, typeK MITE',
		TEMT023 => 'MITE-adh, typeL MITE',
		TEMT024 => 'MITE-adh, typeM MITE',
		TEMT025 => 'MITE-adh, tympN MITE',
		TEMT026 => 'MITE-adh-1',
		TEMT027 => 'MITE-adh-2',
		TEMT028 => 'MITE-adh-3',
		TEMT029 => 'MITE-adh-4',
		TEMT030 => 'MITE-adh-5',
		TEMT031 => 'MITE-adh-6',
		TEMT032 => 'MITE-adh-7',
		TEMT033 => 'MITE-adh-8',
		TEMT034 => 'MITE-adh-9',
		TEMT035 => 'MITE-adh-10',
		TEMT036 => 'MITE-adh-11',
		TEMT037 => 'MITE-adh-12',
		TEMT038 => 'Buhui MITE',
		TEMT039 => 'Casin MITE',
		TEMT040 => 'Centre MITE',
		TEMT041 => 'Delay MITE',
		TEMT042 => 'ECR MITE',
		TEMT043 => 'Helia MITE',
		TEMT044 => 'ID-2 MITE',
		TEMT045 => 'ID-3 MITE',
		TEMT046 => 'ID-4 MITE',
		TEMT047 => 'Lier MITE',
		TEMT048 => 'Stola MITE',
		TEMT049 => 'Stone MITE',
		TEMT050 => 'Susu MITE',
		TEMT051 => 'Wuji MITE',
		TEMT052 => 'Youren MITE',
		TEMT053 => 'Micron MITE',
		TEMT054 => 'Truncator MITE',
		TEMT055 => 'Heart Breaker MITE',
		TEMT056 => 'Frequent Flyer MITE',
		TEMT057 => 'Heart Healer MITE',
		TEMT058 => 'Acrobat MITE',
		TEMT059 => 'mPIF MITE',
		TEMT060 => 'Pangrangja MITE',
		TEMT061 => 'Kiddo MITE',
		TEMT062 => 'mPing/miniSNOOPY MITE',
		TEMT063 => 'MDM MITE',
		TEMTOOT => 'Unclassified MITE',
		# CENTROMERE RELATED SEQUENCES
		CMCM001 => 'Centromere Specific Retrotransposons',
		CMCM002 => 'Centromeric Satellite Repeats',
		CMCMOOT => 'Unclassified Centromere Sequence',
		# TELOMERE RELATED SEQUENCES
		TRTM000 => 'Telomere',
		TRTA000 => 'Telomere Associated',
		# RIBOSOMAL RNA GENES
		RGRR000 => '45s rDNA',
		RGRR005 => '5s rDNA',
		# UNCLASSIFIED
		OTOT000 => 'Unclassified'
		);

    $cl = $TIGRRep{$code} ||
	"UNKNOWN TIGR CODE";
    return $cl;

}

sub GetRBClass
{
#-----------------------------+
# GET REPEAT BASE CLASS       |
#-----------------------------+

# A standardized ontology for the class is returned given the name  
# used by RepBase for the repeat class.

    my $qry = $_[0];           # The qry code sent to the subfunction
    my $cl;                    # The classification that is returned
    my %RBClass;               # Hash to translate from code to class

    %RBClass = (
		#-----------------------------+
		# INTERSPERSED REPEATS        |
		#-----------------------------+
		Interspersed => 'Interspersed Repeat',
		#DNA TRANSPOSONS
		Mariner => 'Mariner (MLE) Transposon',
		hAT => 'DNA Transposon',
		MuDR => 'DNA Transposon',
		EnSpm => 'CACTA, En/Spm Transposon',
		piggyBac => 'DNA Transposon',
		P => 'DNA Transposon',
		Merlin => 'DNA Transposon',
		Harbinger => 'DNA Transposon',
		Transib => 'DNA Transposon',
		Novosib => 'DNA Transposon',
		Mirage => 'DNA Transposon',
		Helitron => 'Helitron',
		Polinton => 'DNA Transposon',
		Rehavkus => 'DNA Transposon',
		DNA => 'DNA Transposon',
		#LTR Retrotransposon
		Gypsy => 'Ty3-gypsy LTR-Retrotansposon',
		Copia => 'Ty1-copia LTR-Retrotransposon',
		LTR => 'LTR Retrotransposon',
		BEL => 'LTR Retrotransposon',
		DIRS => 'LTR Retrotransposon',
		#Endogenous Retrovirus
		ERV1 => 'Endogenous Retrovirus',
		ERV2 => 'Endogenous Retrovirus',
		ERV3 => 'Endogenous Retrovirus',
		#Non-LTR Retrotransposon
		SINE => 'SINE Retrotransposon',
		SINE1 => 'SINE Retrotransposon',
		SINE2 => 'SINE Retrotransposon',
		SINE3 => 'SINE Retrotransposon',
		CRE => 'Non-LTR Retrotransposon',
		NeSL => 'Non-LTR Retrotransposon',
		R4 => 'Non-LTR Retrotransposon',
		R2 => 'Non-LTR Retrotransposon',
		L1 => 'Non-LTR Retrotransposon',
		RTE => 'Non-LTR Retrotransposon',
		I => 'Non-LTR Retrotransposon',
		Jockey => 'Non-LTR Retrotransposon',
		CR1 => 'Non-LTR Retrotransposon',
		Rex1 => 'Non-LTR Retrotransposon',
		RandI => 'Non-LTR Retrotransposon',
		Penelope => 'Non-LTR Retrotransposon',
		# Caulimoviridae 
		Caulimoviridae => 'Caulimoviridae', 
		#-----------------------------+
		# SIMPLE REPEAT               |
		#-----------------------------+
		SAT => 'Satellite',
		MSAT => 'Satellite'
		);

    $cl = $RBClass{$qry} ||
	"UNKNOWN REPBASE NAME";
    return $cl;

}


 

sub trim($)
	 {
	     my $string = shift;
	     $string =~ s/^\s+//;
	     $string =~ s/\s+$//;
	     return $string;
	 }
# Left trim function to remove leading whitespace

sub ltrim($)
	  {
	      my $string = shift;
	      $string =~ s/^\s+//;
	      return $string;
	  }
# Right trim function to remove trailing whitespace

sub rtrim($)
	  {
	      my $string = shift;
	      $string =~ s/\s+$//;
	      return $string;
	  }
