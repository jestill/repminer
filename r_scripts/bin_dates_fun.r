#-----------------------------------------------------------+
# bin_dates.r - classify insertion dates into bins          |
#-----------------------------------------------------------+
# AUTHOR: James C. Estill                                   |
# STARTED: 06/19/2008                                       |
# UPDATED: 06/24/2008                                       |
#-----------------------------------------------------------+
# TO DO: - WRITE THIS AS A FUNCTION THAT TAKES AN ARRAY OR
#          Function should take two column matrix and name as
#          seqData (twoCols, [[1]] is the seq_id, [[2]] variable value)
#        - Make a key with data cut offs

# Remove existing objects
rm(list=ls(all=TRUE));

#-----------------------------+
# LOAD LIBRARIES              |
#-----------------------------+
library(colorRamps);           # Get color ramps

#-----------------------------+
# VARS TO CHANGE              |
#-----------------------------+
numClass = 10;                     # The number of classes to split it into
numBreaks = 200;                   # The number of breaks for drawing the histogram
drawPlot = TRUE;
plotTitle = "Maize LTR Retrotransposon Insertion Date";
plotXLabel = "Insertion Date (YA)";
#plotMaxX = 10000000;              # A X-Axis Maximum value to use for plotting Maximum
plotMaxX = 5000000;              # A X-Axis Maximum value to use for plotting Maximum
#plotMinX = 0;                     # An X-Axis Minimum value used for plotting
#mclSubset = "mcl_1";

#classMethod = "quantile";          # quantile | bagged
#classMethod = "bagged";          # quantile | bagged
classMethod = "kmeans";          # quantile | bagged
#classMethod = "equal_interval";    # equal_interval method
#classMethod = "pam";    # equal_interval method

#-----------------------------+
# SET EXPORT PATHS            |
#-----------------------------+
outColFile = (paste("class_", classMethod, "_", numClass, "_rgb.txt",  sep=""));
outHeader = (paste("class_", classMethod, "_", numClass, sep=""));
outNAFile = (paste("class_", classMethod, "_", numClass, ".NA",  sep=""));
outPlotImage = (paste("class_", classMethod, "_", numClass, ".png",  sep=""));

#-----------------------------+
# COLOR PALETTE               |
#-----------------------------+
colRamp <- blue2red(numClass);
#colRamp <- green2red(numClass);

#-----------------------------+
# DATA IO                     |
#-----------------------------+
inFile = "hq_insertion_date_for_r.txt";
outFile = "hq_insertation_date_bin.txt";
seqData <- read.table(inFile, header=T);

#-----------------------------+
# SUBSET THE DATA             |
#-----------------------------+
if (exists ("mclSubset", mode="any" ) ) {
  sub_mcl <- subset(seqData, seqData$mcl_clust==mclSubset);
  seqData <- sub_mcl;  
};

#-----------------------------+
# GET INFO ABOUT DATA         |
#-----------------------------+
minVal = min(seqData[[2]]);
maxVal = max(seqData[[2]]);

#-----------------------------------------------------------+
# CLASSIFY THE DATA                                         |
#-----------------------------------------------------------+
# TO DO: STD DEV Classification
#        pam -- part of the cluster package .. more robust version of k means
#        ap --- affinity propagation, will need preference vector, will need to do
#-----------------------------+
# KMEANS CLUSTER CLASSIFY     |
#-----------------------------+
# Minimizes the variation of the distances from the mean within each class 
# This is included with R and can be used when the package e1071
# is not installed
if (classMethod == "kmeans") {

  groupK <- kmeans(seqData[[2]], numClass);
  classK <- groupK$cluster;
  
  cols <- match(groupK$cluster, order(groupK$centers));
  ncl <- nrow(groupK$centers);

  stbrks <- matrix(unlist(tapply(seqData[[2]], factor(cols), range)),
                   ncol = 2, byrow = TRUE);
  
}

#-----------------------------+
# BAGGED CLUSTERING           |
#-----------------------------+
if (classMethod == "bagged") {
  
  # blucst requires the e1071 package
  library(e1071);

  # Set a random number seed
  set.seed(20080624);

  newBreaks <- bclust(seqData[[2]],  iter.base=20, centers=numClass, verbose=FALSE);
  cols <- match(newBreaks$cluster, order(newBreaks$centers));
  stbrks <- matrix(unlist(tapply(seqData[[2]], factor(cols), range)),
                   ncol = 2, byrow = TRUE);
}

#-----------------------------+
# QUANTILES                   |
#-----------------------------+
if (classMethod == "quantile") {
  
  quantVal = 1/numClass;
  stbrks <- matrix (NA, ncol=2,nrow=numClass);
  quantRes <- quantile(seqData[[2]], probs = seq (0,1, quantVal), names=FALSE );

  for (i in 1:numClass) {
    # Since quantiles represent percent below the given value we get..
    stbrks [i,1] = quantRes[i]+1; 
    stbrks [i,2] = quantRes[i+1];
  }

  breaks <- matrix (NA, ncol=1,nrow=numClass+1);
  breaks[1] = minVal;
  for (i in 1:numClass) {
    breaks[i+1] = stbrks[i,2];
  }
  
  cols <- findInterval(seqData[[2]], breaks , all.inside=TRUE);
} # End of if classMethod is quantile

#-----------------------------+
# EQUAL INTERVAL              |
#-----------------------------+
if (classMethod == "equal_interval" ) {
  
  stbrks <- matrix (NA, ncol=2,nrow=numClass);
  stepSize = (maxVal - minVal)/numClass;
  breakTop = minVal;              # Bottom of the break
  breakBot = minVal;              # Top of the break
  
  for (i in 1:numClass) {
    breakTop = breakTop + stepSize;
    stbrks [i,1] = breakBot;
    stbrks [i,2] = breakTop;
    breakBot = breakTop + 1;
  }
  
  breaks <- matrix (NA, ncol=1,nrow=numClass+1);
  breaks[1] = minVal;
  for (i in 1:numClass) {
    breaks[i+1] = stbrks[i,2];
  }
  
  cols <- findInterval(seqData[[2]], breaks , all.inside=TRUE);
  
}

#-----------------------------+
# PARTITION AROUND MEDIODS    | 
#-----------------------------+
# Very slow method for a large number of data points
if (classMethod == "pam" ) {

  library("cluster");
  newBreaks <- pam(seqData[[2]],  numClass, diss = FALSE,
                   metric = "euclidean", cluster.only = TRUE
                   );
  
}

#-----------------------------+
# STANDARD DEVIATIONS         |
#-----------------------------+
# Get the mean and do above and below number of std deviations

#-----------------------------------------------------------+
# DRAW THE PLOT                                             |
#-----------------------------------------------------------+
# Changing these to use the min value and max value for limits

if (drawPlot == TRUE ) {

  #-----------------------------+
  # GET X AXIS LIMITS           |
  #-----------------------------+
  if (!exists ("plotMaxX", mode="any" ) ) {
    plotMaxX = max(seqData[[2]]);
  };
  
  if (!exists ("plotMinX", mode="any" ) ) {
    plotMinX = min(seqData[[2]]);
  };

  subText =  (paste("class_", classMethod, "_", numClass,  sep=""));
  
  #-----------------------------+
  # SET UP THE PLOTTING OBJECT  |
  #-----------------------------+
  par <- par(mfrow = c(2,1), pty = "m", bg="white");

  #-----------------------------+
  # TOP HISTOGRAM PLOT          |
  #-----------------------------+
  # HISTOGRAM OBJECTS GIVES INFORMATION FOR EQUAL AREA
  # AGE CLASSES, FOR N BREAKS
  # CAN USE $breaks to delimit and use as a classifier
  # bottom, left, top, right
  par(mar = c (2,4,4,2) + .1);
  #plotTitle 
  histAges <-hist(seqData[[2]], numBreaks, col="blue", xlab=NULL, 
                  main=plotTitle, xlim=c(plotMinX,plotMaxX) );
#  main="Maize LTR Retrotransposon Age", xlim=c(plotMinX,plotMaxX) );

  #-----------------------------+
  # EMPIRICAL CUMULATIVE DISTN  |
  #-----------------------------+
  # Reset values for margin since this does not have a header
  # bottom, left, top, right
  # The following is the default numbers
  #par(mar = c (5,4,4,2) + .1)
  par(mar = c (5,4,1,2) + .1);
  plot.ecdf(seqData[[2]], verticals=TRUE , cex=0.1, 
            main = "", xlab=plotXLabel,
            ylim=c(-0.2,1), xlim=c(plotMinX,plotMaxX), sub = subText );
  rug(seqData[[2]]);
  
  #-----------------------------+
  # ADD COLOR CLASS BOXES       |
  #-----------------------------+
  par(new=T);

  # INITIALIZE MATRICES
  xVal <- matrix ( 1, nrow = numClass, ncol = 1);
  yVal <- matrix ( -0.05, nrow = numClass, ncol = 1);
  #recVals <- matrix(c(NA,NA),
  recVals <- matrix(c(100000,0.1),
                    nrow=numClass, ncol=2, byrow=TRUE);
  # SET VALS
  for (i in 1:numClass) {
    xVal[i] = stbrks[i,1] + ( (stbrks[i,2]- stbrks[i,1])/2 );
    recVals[i,1] =  stbrks[i,2]- stbrks[i,1];
  }

  # DRAW BOXES
  symbols ( xVal, yVal, rectangles = recVals, fg=colRamp , bg=colRamp,
           axes = FALSE, xlab="",ylab="", inches=FALSE,
           ylim=c(-0.2,1), xlim=c(plotMinX,plotMaxX) );

  #-----------------------------+
  # ADD LINES TO ECDF           |
  #-----------------------------+
  xLinesLeft <- matrix (stbrks[,1], nrow = numClass, ncol=2 );
  xLinesRight <- matrix (stbrks[,2], nrow = numClass, ncol=2 );
  yLines <- matrix( c(0,1), nrow=numClass, ncol=2, byrow=TRUE);
  
  # DRAW LINES ITERATIVELY
  for (i in 1:numClass) {
    
    # Left side of breaks
    lines ( c (xLinesLeft[i,1],xLinesLeft[i,2]),
           c (yLines[i,1], yLines[i,2]),
           col = "gray", lty=3);
  
    # Right side of breaks
    lines ( c (xLinesRight[i,1],xLinesRight[i,2]),
           c (yLines[i,1], yLines[i,2]),
           col = "gray", lty=3);
    
  }
  
  par(new=F);
  
} # End of if draw plot

#-----------------------------------------------------------+
# EXPORT THE CLASSIFIED DATA                                |
#-----------------------------------------------------------+

#-----------------------------+
# CREATE OUT MATRIX           |
#-----------------------------+
# This is the classified output matrix
numRows <- length(seqData[[1]]);
#nCols = 2;
classMat <- matrix (NA, ncol=2,nrow=numRows);
classMat[,1] = seqData[[1]];
classMat[,2] = cols;

# CONVERT MATRIX TO DATA FRAME
classFrame <- as.data.frame(classMat);

#-----------------------------+
# WRITE NA FILE FOR CYTOSCAPE |
#-----------------------------+

# WRITE THE HEADER
write ( outHeader, outNAFile, append=FALSE);

# WRITE THE CLASSIFIED DATA
write.table ( classFrame ,outNAFile , sep="=",
             quote=FALSE, row.names=FALSE,
             col.names=FALSE, append=TRUE);

#-----------------------------+
# WRITE RGB TEXT FILE         |
#-----------------------------+
# First need to convert from hex to rgb and then transpose
colRampRGB <- t(col2rgb(colRamp));
write.table (colRampRGB,outColFile );

#-----------------------------+
# WRITE THE PLOT TO PNG FILE  |
#-----------------------------+
if (drawPlot == TRUE ) {
  dev.copy (png, filename=outPlotImage,
            height=600, width=800, bg="white" );
  dev.off();
}

#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
#
# 06/19/2008
# - Started Program
# - Draw histogram of ages
#
# 06/20/2008
# - Added ability to draw ecdf
# - Working on the 
#
# 06/21/2008
# - Working on adding color blocks to ecdf
#
# 06/23/2008
# - Still working on adding color blocks, currently
#   working as vals, now switching to vectors
# - Added lines to ecdf plots
# - Made number of classes a variable for kmeans
#   Line Types
#    (0=blank, 1=solid, 2=dashed, 3=dotted, 4=dotdash,
#     5=longdash, 6=twodash)
#
# 06/24/2008
# - Covert R color codes to rgb codes for Cytoscape
# - Started adding pam .. too slow to work with existing dataset
# - Added option to not draw plot
# - Made plot labels variables
# - Making this a function
# - Filled space up in plot area better making use of 
#   of the par(mar) margin option
# - Added variables for x and y limits on plot
# - Added option to subset the data, example with mcl
#   cluster
#-----------------------------------------------------------+
# JUNKYARD
#-----------------------------------------------------------+
