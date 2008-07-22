#
# SHORT PROGRAM TO EVALUATE DISTIRBUTION OF PID IN ALL BY ALL BLAST
# RESULTS
# J. ESTILL
# 03/29/2007
# 

# Clean up
rm(list=ls(all=TRUE));

# BLAST COLS ARE
# [1] --> QryID
# [2] --> SubID
# [3] --> PID
# [4] --> AlignLen
# [5] --> Mismathces
# [6] --> GapOpen
# [7] --> QryStart
# [8] --> QryEnd
# [9] --> SubStart
# [10] -> SubEnd
# [11] -> EVal
# [12] -> BitScore

print ("Loading the BLAST output file.", quote=FALSE);
BlastFile = "/home/jestill/projects/SanMiguel/SanMig_SanMig_e10_m8.blo";


#AbA <- as.matrix(read.table(BlastFile, header=F, as.is=TRUE));
# The above code appears to load the PIDs as text instead of numbers
# The following appears to treat the 
#AbA <- as.matrix(read.table(BlastFile, header=F));

# Do scan test with 100 lines
AbA <- scan(BlastFile,
            nlines=10,
            what=c(
              "character", #[1]
              "character", #[2]
              "numeric",   #[3]
              "integer",   #[4]
              "integer",   #[5]
              "integer",   #[6]
              "integer",   #[7]
              "integer",   #[8]
              "integer",   #[9]
              "integer",   #[10]
              "numeric",   #[11]
              "integer"    #[12]
              ),
            flush=TRUE
            );

#NumBreaks = 100;

#HistPlot <-hist(AbA[,3],
#                breaks=NumBreaks,
#                col="blue",
#                main="PID Distribution for All by All BLAST",
#                xlab="PID",
#                ylab="FREQUENCY",
#                xlim=c(0,100));

