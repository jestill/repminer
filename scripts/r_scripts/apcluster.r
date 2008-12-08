#-----------------------------------------------------------+
#                                                           |
# AFFINITY PROPAGATION IN R                                 |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: JAMES C. ESTILL                                   |
# STARTED: 07/04/2008                                       |
# UPDATED: 07/08/2008                                       |
# CONTACT: JamesEstill_at_gmail.com                         |
#                                                           |
# DESCRIPTION:                                              |
#  An implementation of the affinity propagation clustering |
#  algorithm in R. Takes as its input a similarity file     |
#  of three columns representing a matrix of similarity     |
#  between elements. Input file is tab delimited text file  |
#  with no header and three columns:                        |
#   col 1. Row ID (Integer)                                 |
#   col 2. Column ID (Integer)                              |
#   col 3. Metric of similarity betwee row and column       |
#  The total number of elements in the similarity matrix    |
#  will be set to be the max of col 1 and col 2.            |
#                                                           |
#  Output is a text file listing integer values of the      |
#  exempars identified.                                     |
#                                                           |
# REFERENCES                                                |
#  Frey, B. J. and D. Dueck (2007). "Clustering by passing  |
#  messages between data points." Science 315(5814): 972-6. |
#                                                           |
#  Frey, B. J. and D. Dueck (2007).                         |
#  Mixture Modeling by Affinity Propagation.                |
#                                                           |
#  http://www.psi.toronto.edu/affinitypropagation/          |
#                                                           |
#-----------------------------------------------------------+
#
# TO DO:
#    and can check is.numeric for the first two columns to
#    decide. This can return the results in a way that is
#    easy to interpret.
#  - I am currently stopping this when the exemplar vector
#    does not change, may want to consider stopping this
#    when the index vector does not change, however this
#    may make for interesting options
#  - Add details option
#  - Add warning when a large N
#  - With the use of -Inf, will want to replace with -realmax
#  - Speed up calculation of S if possible
#  - Add option to send final plot to an external file
#  - Add memory test to determine approximate amount of memory
#    needed
#  - Add messages for the verbose mode
#  - Sparse matrix implementation=
#  - Make a global median value the default preference
#  - do object.size() test and report for error tracking
#    this can take the object.size of S and compare it to
#    the available memory
#  - Perhaps do this by connected compents for large matrices
#  - Pass data labels to the function
#  - Allow set designation by an external classifier column
#    (ie. connected components),
#    and sort out these and run apclust for each set ()
#  - Get vector of exemplars and K, send the results to
#    a K mediods program.
#  - Replace NA with -realmax or some other value to represent
#    not connected
#
#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
#inFile = "ToyProblemSimilarities.txt";
inFile = "ToyProblemSimilarities_nas.txt";
#inFile = "ltr_ltr.txt";
test_simData <- read.table(inFile, header=F);
#prefFile = "ToyProblemPreferences.txt";

#inFile = "TravelRoutingSimilarities.txt";
#inFile = "short_test.txt";
outFile = "ToyProblemOutput.txt";

#test_lam <- 0.7;                       # DAMPING FACTOR
#test_lam <- 0.5;                       # DAMPING FACTOR
test_lam <- 0.7;                       # DAMPING FACTOR
test_maxits <- 100;                   # NUMBER OF ITERATIONS
test_prefVal <- -16;            # Preference value for testing pref of length 1
#test_prefVector <- c (-16,-16,-16,-16,-16,
#                      -16,-16,-16,-16,-16,
#                      -16,-16,-16,-16,-16,
#                      -16,-16,-16,-16,-16,
#                      -16,-16,-16,-16,-16);



# the following test vector is the row medians
# for the toy data set.
test_prefVector <- c (-25.0191,
                      -15.97939,
                      -14.63548,
                      -16.84353,
                      -36.09047,
                      -24.74661,
                      -16.12498,
                      -27.14535,
                      -27.57591,
                      -15.41481,
                      -8.91651,
                      -12.70332,
                      -17.74917,
                      -11.47763,
                      -10.44418,
                      -11.42128,
                      -24.83111,
                      -30.56645,
                      -28.38697,
                      -16.98507,
                      -12.91641,
                      -11.12068,
                      -11.98069,
                      -16.09133,
                      -19.40381);


#test_simVector <- c(1,4,6,9,10,11,70,90,120);
#test_simData <- read.table(inFile, header=F);
#-----------------------------+
# TEST INSERTION DATE VECTOR  |
#-----------------------------+
inFile = "InsertDateShort.txt";
test_simVector <- as.vector(  (read.table(inFile, header=T))[,2] );

#-----------------------------------------------------------+
# AFFINITY PROPAGATION FUNCTION                             |
#-----------------------------------------------------------+
apcluster <- function (s, p, maxits=1000, dampfact=0.9, convits=100,
                       verbose=FALSE, noise=TRUE, do.plot=TRUE) {

  #-----------------------------+
  # MACHINE VARIABLES           |
  #-----------------------------+
  # These are named with the equivalent Matlab special varaible name
  start <- proc.time();
  eps <- .Machine$double.eps;
  realmin <- .Machine$double.xmin;
  realmax <- .Machine$double.xmax;
  # Rename variables for clarity when using below
  simData <- s;
  prefVal <- p;
  lam <- dampfact;

  #oldE <- NA;
  conCount <- 0; # The number of times exemplars are the same

  #-----------------------------+
  # FUNCTION TO REPLACE NAs     |
  #-----------------------------+
  subRealmax <- function () {
    returnVal <- -realmax;
    list(val=returnVal);
  };
  
  #-----------------------------------------------------------+
  # GENERATE THE MATRIX S                                     |
  #-----------------------------------------------------------+
  if (verbose == TRUE) {
    print ("Generating Similarity Matrix")
  }
  # a single vector that is used to extract pairwise distances
  # using basic Euclidian distance


  # simData is the similarity data file
  # three col file
 
  # Total number of vals is the max of col 1 and col 2
  #simDataDim <- dim(simData);
  #simDataCol <- simDataDim[2];
  
  #-----------------------------+
  # GET DIMENSIONS OF simData   |
  #-----------------------------+
  if (is.matrix(simData) == TRUE ) {
    simDataRow <- dim(simData)[1];
    simDataCol <- dim(simData)[2];
  }
  else if (is.data.frame(simData) == TRUE) {
    simDataRow <- dim(simData)[1];
    simDataCol <- dim(simData)[2];
  }
  else if (is.vector(simData == TRUE)) {
    simDataRow <- length(simData);
    simDataCol <- 1;
  }
  
  #-----------------------------+
  # THREE COLUMNS OF DATA       |
  #-----------------------------+
  # col 1 is row id num
  # col 2 is col id num
  # col 3 is similarity value
  if (simDataCol == 3) {

#    # ADD CHECK THAT THE VALUES ARE NUMBERS
#    
#    #-----------------------------+
#    # DETERMINE MAX ROW OR COL    |
#    #-----------------------------+
#    N <- max(simData[[1]], simData[[2]]);
#    #-----------------------------+
#    # LOAD SIMILARITY MATRIX      |
#    #-----------------------------+
#    S <- matrix( NA, nrow=N, ncol=N);
#
#    # TO DO: After getting median values, can put in the -INF values
#    #        -Inf it appears to be in R
#  
#    # LOAD INPUT DATA INTO THE MATRIX
#    numInputRows <- length(simData[[1]]);
#    for (i in 1:numInputRows) {
#      S[ simData[i,1], simData[i,2] ] = simData[i,3];
#    }

#
# USING as.matrix for xtabs still gives a slightly unexpected result   
#////////////////////////
# BEGIN TEST CODE
    
    # CAN ALSO USE THE XTABS FUNCTION
    # WILL XTABS ALWAYS CREATE A SQUARE MATRIX ?

    # This also gives 
    # The following will work with data frames that have been named
    # however, values that do not occured as counted as zeros
    #S <- as.matrix (xtabs(simData[,3] ~ simData[,1] + simData[,2]) )



#//////////////////////////////////
    S <- tapply(simData[,3], simData[,c(1,2)], c );
    dataLabels <- rownames(S);
    dimnames(S) <- NULL;
    N <- dim(S)[1];
    #print(S);
    #stop("test stop");
    # CHECK THAT THE MATRIX IS SQUARE
    if (dim(S)[1] != dim(S)[2]) {
      stop("Similarity matrix is not square.");
    }

    

    
#    # TEMP FOR DEBUG AND GIVE TIME
#    print ("SimilarityMatrix")
#    print (S);
##    print("TME")
##    totalTime <- proc.time() - start;
##    print (totalTime);
#    stop("Working on nas");

# END TEST CODE
#\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
  }



  
  #-----------------------------+
  # SQUARE MATRIX               |
  #-----------------------------+
  else if (simDataCol == simDataRow) {
    N <- simDataRow;
    S <- as.matrix(simData);
  }

  #-----------------------------+
  # SINGLE VECTOR               |
  #-----------------------------+
  else if (simDataCol == 1) {

    if (is.numeric(simData) == FALSE) {
      stop ("Input data vector is not numeric");
    }
    
    N <- simDataRow;

    # MAKE MATRIX TO HOLD SIMILARITES
    
    # DETERMINE PAIRWISE SIMILARITIES
    # NEED TO MAKE THIS FASTER
    # Using negative euclidian distance for similarities
    # CAN JUST DO j<i and then flip the matrix

    if (verbose == TRUE) {
      print ("Allocating Matrix")
    }
    S <-  matrix( NA, nrow=N, ncol=N);

    if (verbose == TRUE) {
      print ("Filling Matrix")
    }

    # THE FOLLOWING WORKS WITHOUT DUPLICATING
    # RESULTS, BUT
    # cannot allocate vector of size 1403440 Kb
    # This appears to be faster then using the dist funciton
#    for (i in 1:N) {
#      for (j in 1:N) {
#        S[i,j] <- -( abs(simData[i]-simData[j]) );
#      }
#    }


    #////////////////////////////////////
    # 07/09/2008 
    # TRY THE FOLLOWING
    S <- as.matrix( -1*dist(simData)  );
    dataLabels <- rownames(S);
    dimnames(S) <- NULL;
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\    
    
  }

  #-----------------------------+
  # NOT AN EXPECTED FORMAT      |
  #-----------------------------+
  else {
    stop("Similarity data is not in an expected format");
  }

  
  #-----------------------------+
  # REPLACE NA'S WITH -realmax? |
  #-----------------------------+
  # replace NA's with -realmax
  
  
  #-----------------------------------------------------------+
  # SET THE PREFERENCE VALUES                                 |
  #-----------------------------------------------------------+
  # If I can't make sense of the values passed, will use a
  # global median value. Should also throw a warning
  if (verbose == TRUE) {
    print ("Generating Preference Values")
  }
  
  if (is.numeric(prefVal) == TRUE) {

    prefLen <- length(prefVal);        # Length of preference Value
    
    #-----------------------------+
    # PREFERENCES AS SINGLE VALUE |
    #-----------------------------+
    if (prefLen == 1) {
    
      # Need to check that the preference value is a number
      # and not a text syntax error
      for (i in 1:N) {
        S[i,i] = prefVal;
      }
    }

    #-----------------------------+
    # PREFERENCE AS VECTOR        |
    #-----------------------------+
    else if (prefLen > 0) {
      # Need to add a check to make sure that the length of the preference
      # values is equal to the length of the
      for (i in 1:N) {
        S[i,i] = prefVal[i];
      }
    }
  }
  else {

    if (prefVal == "rowmedian") {
      
      prefVal <- vector(mode="numeric",length=N );
      #prefVal <- as.vector(prefVal);
      for (i in 1:N) {
        S[i,i] <- median (S[i,], na.rm=TRUE);
        prefVal[i] <- S[i,i];
        print (S[i,i]);
      }      
    }
    
#    #-----------------------------+
#    # PREFERENCES AS GLOBAL MEDIAN|
#    #-----------------------------+
#    if (prefVal == "median") {
    else {
      print ("GLOBAL MEDIAN");
      globalMedian <- median(S, na.rm=TRUE);
      for (i in 1:N) {
        S[i,i] <- globalMedian;
      }
      prefVal <- globalMedian;
      
    }
    
  } # END OF SET PREFERENCE VALUES


  #-----------------------------------------------------------+
  # REPLACES NA VALUES WITH -realmax .. very distant          |
  #-----------------------------------------------------------+
  # print the S values
  #print("S");
  #print(S);
  toReplace <- (is.na(S));
  S[toReplace] <- -realmax;
  #print("replaced S");
  #print(S);
  #stop("debug");

  #-----------------------------+
  # DEBUG:                      |
  # PRINT SIMILARITY MATRIX     |
  # AND PREFERENCE VALUES       |
  #-----------------------------+
  #print (prefVal);
  #print (S);
  #stop("peace out");
  
  #-----------------------------+
  # INITIALIZE MATRICES         |
  #-----------------------------+
  zeroMat <- matrix(0, nrow=N, ncol=N);
  Rp <- matrix(NA, nrow=N, ncol=N);

  #-----------------------------+
  # INITIALIZE MESSAGES         |
  #-----------------------------+
  A <- matrix( 0, nrow=N, ncol=N);
  R <- matrix( 0, nrow=N, ncol=N);

  # summaryVal matrix to hold summary values
  # col 1. iteration number
  # col 2. net similarity
  summaryVal <- matrix( NA, nrow=maxits, ncol=2);

  #-----------------------------+
  # REMOVE DEGENERACIES         |
  #-----------------------------+
  # Adds tiny bit of random numbers to matrix
  # For cases where there are many ties, may want to use a multiplication
  # greater then 100 here.
  if (noise == TRUE) {
    S <- S+( matrix(eps, nrow=N, ncol=N)%*%S   
            + matrix(realmin, nrow=N, ncol=N)%*%matrix(100,nrow=N,ncol=N)
            ) * matrix(runif(N^2),N);
  }
  
  for (iter in 1:maxits) {
    
    #-----------------------------+
    # COMPUTE RESPONSIBILITIES    |
    #-----------------------------+
    Rold <- R;
    AS <- A+S;
    
    I <- max.col(AS, ties.method=c("first"));  
    Y <- apply(AS,1,max);
    for (i in 1:N) {
      AS[i,I[i]]= -realmax;
    }
    
    I2 <- max.col(AS, ties.method=c("first"));
    Y2 <- apply(AS,1,max);
    
    R <- S-kronecker(matrix(1,1,N),Y);
    
    for (i in 1:N) {
      R[i,I[i]] <- S[i,I[i]] - Y2[i];
    }
    
    #-----------------------------+
    # DAMPEN RESPONSIBILITIES     |
    #-----------------------------+
    R <- matrix( (1-lam), nrow=N, ncol=N)*R + matrix(lam,nrow=N,ncol=N)*Rold;
  
    #-----------------------------+
    # COMPUTE AVAILABILITIES      |
    #-----------------------------+
    Aold <- A;
    
    Rp <- pmax(R,zeroMat);
    
    for (k in 1:N) {
      Rp[k,k] <- R[k,k];
    }
    
    A <- t(kronecker ( matrix (1,1,N), apply(Rp,2,sum) ) ) - Rp;
    
    dA <- diag(A);
    
    A <- pmin(A,zeroMat);
    for (k in 1:N) {
      A[k,k] <- dA[k];
    }
    
    #-----------------------------+
    # DAMPEN AVAILABILITIES       |
    #-----------------------------+
    A <- (1-lam)*A+lam*Aold;

    #-----------------------------------------------------------+
    # BEGIN NEW WORK
    #-----------------------------------------------------------+
    
    #-----------------------------+
    # CHECK FOR CONVERGENCE       |
    #-----------------------------+
    # E=((diag(A)+diag(R))>0);
    # The following gives true false values
    #new_E <- (diag(A) + diag(R)) > 0;
    # The following gives the indexes
    #new_E <- which((diag(A) + diag(R)) > 0);

    # Adding as.array to deal with the
    # cases when tmpE has only a single exemplar
    tmpE <- which((diag(A) + diag(R)) > 0);
    if (length(tmpE) == 1) {
      tmpEVal <- (A[tmpE,tmpE]) + R[tmpE,tmpE];
      tmpK <- sum(tmpEVal);
    }else {
      tmpEVal <- diag(A[tmpE]) + diag(R[tmpE]);
      tmpK <- sum(tmpEVal);
    }
    
    # e(:,mod(i-1,convits)+1)=E;
    # K=sum(E);
  
    #-----------------------------+
    # CURRENT MEASURE OF NET SIM  |
    #-----------------------------+
    #if new_K==0
    #          tmpnetsim=nan; tmpdpsim=nan; tmpexpref=nan; tmpidx=nan;
    if (tmpK == 0) {
      #cat ("NONE\n");
      tmpnetsim <- NA;
      tmpdpsim <- NA;
      tmpexpref <- NA;
      tmpidx <- NA;
      
      # Can stick in a value here if desired
      summaryVal[iter,1] = iter;
      #summaryVal[iter,2] = netSim;
      
      if (verbose == TRUE) {
        print ( c(iter, "NONE") );
      }
      
    }
    
    if (tmpK != 0) {
      # tmp message for debug
      #print ("SOME\n");
      #print ( c(iter, "SOME") );
  
      # NONE OF THE FOLLOWING COMMENTED
      # THIS IS THE SAME AS BELOW
      E <- R+A;
      I <- which(diag(E) > 0)
      K <- length(I);
      if ( K > 1) {
        tmp <- apply(S[,I],1,max);
        c_val<- max.col(S[,I]);
      }
      if ( K == 1) {
        tmp <- max(S[,I]);
        c_val <- max.col(S[,I]);
      }
      c_val[I] <- 1:K;
      idx <- I[c_val];
      netSim = 0;
      for ( i in 1:N) {
        indVal = S[ idx[i] , i ];
        netSim = netSim + indVal;
      }

      # The following deals with a problem when -realmax
      # is the value, this appears to cause problems plotting
      if (netSim == -realmax) {
        summaryVal[iter,1] = iter;
      }
      else {
        summaryVal[iter,1] = iter;
        summaryVal[iter,2] = netSim;
      }
 
      #-----------------------------------------------------------+
      # PLOT AS NEW PLOT OBJECT
      #-----------------------------------------------------------+
      if (do.plot == TRUE) {
        if (is.na(summaryVal[iter,2]) == FALSE) {
          plot(summaryVal, type="l", col="red",
               xlab="Iteration", ylab="Net Similarity",
               main = "Affinity Propagation Net Similarity Score");
        }
      }
      
      if (verbose == TRUE) {
        # iteration number and netsimilarity
        print ( c(iter, netSim) );
      }
      
    }


    #-----------------------------+
    # CHECK CONVERGENCE COUNTER   |
    #-----------------------------+
    # Convergence when the exemplars are the same, however
    # the net similarity and cluster membership
    # may be changing
    if (exists ("oldE") == TRUE) {
      if (identical (tmpE, oldE) == TRUE ) {
        iter <-  maxits;
        conCount <- conCount + 1;
 
        #print ( c("Con Count", conCount));
        
        if (conCount == convits) {
          break;
        }

      } else {
        conCount <- 0;
      }      
    }

    if (tmpK != 0) {
      oldE <- tmpE;
    }
    
    #-----------------------------------------------------------+
    # END NEW WORK                                              |
    #-----------------------------------------------------------+
    
    
  }

  #-----------------------------+
  # PSEUDOMARGINALS             |
  #-----------------------------+
  E <- R+A;

  #-----------------------------+
  # INDICIES OF EXEMPLARS       |
  #-----------------------------+
  I <- which(diag(E) > 0);
  # trying below to get rid of double information
  #I <- which( diag(E) > 0, arr.ind = TRUE );
  K <- length(I);
  #-----------------------------+
  # ASSIGNMENTS                 |
  #-----------------------------+
  # Assign seqs to their exemplars
  # There will be casees where I is a vector,
  # and there will be times when I is a single value.
  # Can use K to determine
  
  if ( K > 1) {
    tmp <- apply(S[,I],1,max);
    c_val<- max.col(S[,I]);
  }

   # TO DO: THE CASE WHERE THERE IS A SINGLE EXEMPLAR
  if ( K == 1) {
    tmp <- max(S[,I]);
    c_val <- max.col(S[,I]);
  }
  
  c_val[I] <- 1:K;
  idx <- I[c_val];

  #-----------------------------+
  # FINAL NET SIMILARITY        |
  #-----------------------------+
  # The following works, this should be sped up to
  # matrix based methods 
  netSim = 0;
  for ( i in 1:N) {
    indVal = S[ idx[i] , i ];
    netSim = netSim + indVal;
  }

  #-----------------------------+
  # TOTAL TIME                  |
  #-----------------------------+
  # Gives
  # [1] User Time
  # [2] System Time
  # [3] Total Elapsed Time
  # [4] Sum of child process user time
  # [5] Sum of child process system time
  totalTime <- proc.time() - start;
  print (totalTime);
  
  #-----------------------------+
  # FINAL PROGRAM STATUS        |
  #-----------------------------+
  # This will always print even when not verbose
  # TREATING CAT LIKE STD OUT

  # I seems to contains the names AND the index values
  print ("I IS:")
  print (I);

  #-----------------------------+
  # DEBUG INFO                  |
  #-----------------------------+
  
  # The following does not work, need to use
  #exemplarLabels<- rownames(S[c(I)]);
  print("LABELS");
  print(dataLabels);
  test_dim <- dim(I);
  print ("DIM");
  print (test_dim);
  print("K");
  print (K);
  #print("E");
  #print(E);
  print ("DONE EXEMPLARS");


  #-----------------------------+
  # RETURN ANSWERS              |
  #-----------------------------+
  # This will return data labels if they are available
  if (exists("dataLabels")) {
    
    # GET LABELS 
    # col 1. data label
    # col 2. label of exemplar data point
    exemplarLabels <- matrix(NA,nrow=N,ncol=2);
    exemplarLabels[,1] <- dataLabels;
    exemplarLabels[,2] <- dataLabels[idx];
    #print(exemplarLabels);

    # RETURN RESULTS OF THE FUNCTION
    list(summary=summaryVal,
         labels=exemplarLabels,
         netsimilarity=netSim,
         index=idx,
         preferences=prefVal,
         time=totalTime,
         exemplars=I,
         clusters=K);
    
  }
  else {

    # RETURN RESULTS OF THE FUNCTION
    list(netsimilarity=netSim,
         index=idx,
         preferences=prefVal,
         time=totalTime,
         exemplars=I,
         clusters=K);
        
  }
  
# The following also gives double sets for exemplars and index  
#  pairlist(similarity=netSim,
#           index=idx,
#           preferences=prefVal,
#           time=totalTime,
#           exemplars=I,
#           numClusters=K);

  
} # END OF apcluster function


#-----------------------------------------------------------+
# TEST THE FUNCTION
#-----------------------------------------------------------+
#apcluster (test_simData, test_prefVal,
#           test_maxits, test_lam, verbose=TRUE);

#myAnswer <- apcluster (test_simData, test_prefVal,
#                       test_maxits, test_lam,
#                       #do.plot = FALSE,
#                       verbose=TRUE);
# Anothe way to write this is

# convits - if exemplars stay the same for this number
#           of iterations then do not change

myAnswer <- apcluster (s = test_simData,
                       #s = test_simVector,
                       dampfact = 0.98,
                       p = "median",
                       maxits = 1000,
                       #convits = 20,
                       #p = "rowmedian",
                       #p= test_prefVector,
                       #p = -16,
                       #p = test_prefVal,
                       #maxits = test_maxits,
                       #dampfact = test_lam,
                       noise=FALSE,
                       #do.plot = FALSE,
                       do.plot = TRUE,
                       verbose=TRUE);

#-----------------------------------------------------------+
# WRITE OUTPUT                                              |
#-----------------------------------------------------------+
# idx is the index value
# The web page provides 5 files
# 1. ap.stop
#     This is the exit code of the ap program
#     Contains 0 for normal ending
# 2. idx.txt
#     For N data points, this has N integers.
#     The ith integer is the index of the data point
#     that is the exemplar for the data point.
# 3. netsim.txt
#     Net similarity achieved after each iteration
# 4. netsim2.txt
#     Net similarity achieved after each iteration
#     This on includes the iteration index before each
#     net similiarity value.
# 5. summary.txt
#     The training conditions specified by the user
#      (ie. max number of iterations, damping factor)
#     as well as summary of the ap
#      (number of iterations, number of clusters,
#       final net similarity, run time)
#
# I want to add the following
# labeled, allow for users to add a forth column giving the
# label for that data point. This should be a unique identified
# and will allow for easier interpretation of the results.

# idx is the final answer, the list of vals that an
# exemplar points to .. however this does not seem to work

#-----------------------------------------------------------+
# CHANGELOG                                                 |
#-----------------------------------------------------------+
# 07/04/2008
# - Started program
# - Input of three col file into the matrix S
#
# 07/05/2008
# - Added program information to header with references
# - Switched to using the following as opposed to iterating across matrix
#   I <- max.col(AS);
#   Y <- apply(AS,1,max);
# - Finished first draft of translated code, need to debug
#
# 07/07/2008
# - Running line by line test along side the matlab program
#   using a test dataset
# - Switched max.col from the default of picking ties at
#   random to always picking the first. This gives results
#   that are consistent with MATLAB
# - Working version of the program completed
# - Added the ability to use different types of preferences:
#      - single value
#      - vector of values equal to number of rows,cols
#      - global median
#      - row median
# - Added calculation of netSimilarity as determined
#   by an iterative method
# - Added net similarity
# - Added plot of net similarity
#
# 07/08/2008
# - Added verbose mode
# - Writing as an R function: apcluster()
#      apcluster (s, p, maxits, lam,
#                       verbose=FALSE, addNoise=TRUE,
#                       do.plot=TRUE)
# - Added do.plot as variable
# - Added alternative to get tmpK and tmpEVal when
# - Added tracking of the time it takes to run the function
#  
# 07/09/2008
# - Added option for adding noise
# - Added option for s to be a square matrix or data column
# - Added optin for s to a a vector of numbers
# - Added convits
#    --- if the estimated exemplars stay the same for
#        convit iterations, then stop
# - Renamed numIter to maxits
# - Added code to determine dimensions of simData to see if it
#   is a matrix , three cols or a single col of data
# - Added ability to pass named objects in three colums
#   instead of numbered 1:N in three colums, this makes
#   use of the R xtabs function.
#
# 07/10/2008
# - Adding the ability to return list of labels and the exemplars
#   as labels when labels are available. Returned as $labels
# - Added summary data matrix to the returned results
# 07/11/2008
# - Switched from using xtab to using tapply to get
#   data from labeled three column table, This seems to
#   be the way to get NA instead of 0 for empty cells
# - The NAs are later changed to -realmax values, this
#   is to indicate that they are not connected, it would
#   perhaps be possible to use -Inf as well, but the
#   AP program in matlab used -realmax
# - Using -realmax values seems to causes and infinit
#   axis error with plot so I swith -realmax with NA
#   in the summaryVal matrix
# - 
