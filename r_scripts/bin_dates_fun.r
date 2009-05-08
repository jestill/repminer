#!/usr/bin/r
#-----------------------------------------------------------+
# bin_dates.r - classify insertion dates into bins          |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# STARTED: 06/19/2008                                       |
# UPDATED: 05/07/2009                                       |
# VERSION: $Rev$                                            |
#-----------------------------------------------------------+
#
# DESCRIPTION:
#  Draws a histogram and emprical cumulative distribution
#  function with color classificaiton along the x-axis.
#  The text files that are produces allow this color data
#  to be used on classified objects in Cytoscape. This
#  script is written to accept variables from the command
#  line and can therefore be used from a perl program.
#
# AVAILABLE ARGS (default vals) - [options]:
#  --infile      Path to the intput file, expects a header in input file (STDIN)
#  --outdir      Path to the outfile that is produced (STDOUT)
#  --name        The base name to name the oufiles with
#  --num-class   The number of classes to generate (10)
#  --num-breaks  The number of breaks to draw in the histogram (200)
#  --class       The classification method to use (kmeans):
#                [quantile | bagged | kmeans | equal_interval | pam ]
#  --height      Height of the that is created
#  --width       Width of the image that is created
#  --title       Title of the graph
#  
# BOOLEANS
#  --plot        Draw the plot
#  --inv-color   Invert the color (Red Smaller Number, Blue Larger Number)
#  --verbose     Run the script in verbose mode
#
#
# EXAMPLE USE:
# R < bin_dates_fun.r --no-save --slave
# R < bin_dates_fun.r --no-save --slave --args --name dude
# R < bin_dates_fun.r --no-save --slave --args --num-class 20 --num-breaks 100 --infile test_insertion_date.txt --outdir test --verbose --plot --inv-color
#R < bin_dates_fun.r --no-save --slave --args --infile test_insertion_date.txt --outdir test2 --plot --title "Monkey underpants weight" --verbose --labelx today --max-x 5000000 --min-x 200000

# Remove existing objects
# This may be a terrible idea in some circumstances

rm(list=ls(all=TRUE));

#-----------------------------+
# LOAD LIBRARIES              |
#-----------------------------+
library('colorRamps');         # Get color ramps

#-----------------------------------------------------------+
# GET ARGS FROM THE COMMAND LINE                            |
#-----------------------------------------------------------+

#-----------------------------+
# GET ARGS FUNCTION           |
#-----------------------------+
# The GetArg subfunction allows for arguments to be passed
# to the R script from the command line
# Returns null when string not found in command lines
GetArg <-function(strSwitch) {
  
  ComArgs <- commandArgs();
  
  numArgs <- (length(ComArgs));
  myVar <- NULL;
  for (i in 1:numArgs) {
    if (ComArgs[i] == strSwitch) {
      myVar <- ComArgs[i+1];
    }
  }
  
  GetArg <- myVar;

}

#-----------------------------+
# GET BOOLEANS FUNCTION       |
#-----------------------------+
GetBool <- function(bolFlag) {
  
  BoolArgs <- commandArgs();
  
  bolVar <- FALSE;
  numArgs <- (length(BoolArgs));
  for (i in 1:numArgs) {
    if (BoolArgs[i] == bolFlag) {
      bolVar <- TRUE;
    }
  }
  
  GetBool <- bolVar;
  
}

#-----------------------------+
# BOOLEANS                    |
#-----------------------------+
# Draw the plot
drawPlot <- GetBool('--plot');
# Inverse the color ramp
inverseColor <- GetBool('--inv-color');
# Run script in verbose mode
isVerbose <- GetBool('--verbose');

#-----------------------------+
# INPUT FILE                  |
#-----------------------------+
# If no infile given, expect input from STIN
#inFile <- "test_insertion_date.txt";
argInFile <- GetArg('--infile');
if (is.null(argInFile) == FALSE) {
  inFile <- argInFile;
  seqData <- read.table(inFile, header=T);
} else {
  seqData <- read.table("/tmp/Rd.txt",TRUE);
}

#-----------------------------+
# OUTPUT DIR                  |
#-----------------------------+
argOutDir <- GetArg('--outdir');
if (is.null(argOutDir) == FALSE) {
  outDir <- argOutDir;  
}

#-----------------------------+
# NUMBER OF CLASSES           |
#-----------------------------+
# The number of classes to split it into
numClass <- 10;
argNumClass <- GetArg('--num-class');
if (is.null(argNumClass) == FALSE) {
  numClass <- as.numeric(argNumClass);
}

#-----------------------------+
# NUMBER OF BREAKS            |
#-----------------------------+
# The number of breaks for drawing the histogram
numBreaks = 200;
argNumBreaks <- GetArg('--num-breaks');
if (is.null(argNumBreaks) == FALSE) {
  numBreaks <- as.numeric(argNumBreaks);

}

#-----------------------------+
# CLASSIFICATION METHOD       |
#-----------------------------+
classMethod <- "kmeans";
argClassMethod <- GetArg('--class');
if (is.null(argClassMethod) == FALSE) {
  classMethod <- argClassMethod;
}

#-----------------------------+
# BASE NAME                   |
#-----------------------------+
baseName <- "class";
argBaseName <- GetArg('--name');
if (is.null(argBaseName) == FALSE) {
  baseName <- argBaseName;
}

#-----------------------------+
# PLOT HEIGHT                 |
#-----------------------------+
plotHeight <- 600;
argPlotHeight <- GetArg('--height');
if (is.null(argPlotHeight)==FALSE) {
  plotHeight <- as.numeric(argPlotHeight);
}

#-----------------------------+
# PLOT WIDTH                  |
#-----------------------------+
plotWidth <- 800;
argPlotWidth <- GetArg('--width');
if (is.null(argPlotWidth)==FALSE) {
  plotWidth <- as.numeric(argPlotWidth);
}

#-----------------------------+
# PLOT TITLE                  |
#-----------------------------+
plotTitle = "Maize LTR Retrotransposon Insertion Date";
argPlotTitle <- GetArg('--title');
if (is.null(argPlotTitle)==FALSE) {
  plotTitle <- argPlotTitle;
}

#-----------------------------+
# PLOT X-AXIS LABEL           |
#-----------------------------+
plotXLabel <- "Insertion Date (YA)";
argPlotXLabel <- GetArg('--labelx')
if (is.null(argPlotXLabel)==FALSE) {
  plotXLabel <- argPlotXLabel;
}

#-----------------------------+
# PLOT X-AXIS MAX             |
#-----------------------------+
argPlotMaxX <- GetArg('--max-x');
if (is.null(argPlotMaxX)==FALSE) {
  plotMaxX <- as.numeric(argPlotMaxX);
}

#-----------------------------+
# PLOT X-AXIS MIN             |
#-----------------------------+
argPlotMinX <- GetArg('--min-x');
if (is.null(argPlotMinX)==FALSE) {
  plotMinX <- as.numeric(argPlotMinX);
}

# PRINT VAR VALUES IF VERBOSE
if (isVerbose == TRUE) {
  cat ("Infile Path:\t",inFile,"\n");
  cat ("Outdir Path:\t",outDir,"\n");
  cat ("Num Breaks:\t",numBreaks,"\n");
  cat ("Num Classes:\t",numClass,"\n");
  cat ("Class Method:\t",classMethod,"\n");
  cat ("Plot Height:\t",plotHeight,"\n");
  cat ("Plot Width:\t",plotWidth,"\n");
  cat ("Plot Title:\t",plotTitle,"\n");
  cat ("Plot X-Label:\t",plotXLabel,"\n");

  if (exists ("plotMaxX", mode="any" )) {
      cat ("Plot Max-x:\t",plotMaxX,"\n");
    }
  
  if (exists ("plotMinX", mode="any" )) {
    cat ("Plot Max-x:\t",plotMinX,"\n");
  }
      
}


# VARS 
#mclSubset = "mcl_1";


#-----------------------------+
# SET EXPORT PATHS            |
#-----------------------------+
outColFile = (paste(baseName,"_", classMethod, "_", numClass, "_rgb.txt",  sep=""));
outNAFile = (paste(baseName,"_", classMethod, "_", numClass, ".NA",  sep=""));
outPlotImage = (paste(baseName,"_", classMethod, "_", numClass, ".png",  sep=""));

# Prepend OutDir if given
if (exists ("outDir", mode="any" ) ) {

  # Create dir if it does not exists
  if (!file.exists(outDir)) {
    dir.create(outDir, recursive=TRUE);
  }

  # Using file.path will create the file path in a platform independent way
  outColFile <- file.path(outDir, outColFile);
  outNAFile <- file.path(outDir,outNAFile);
  outPlotImage <- file.path(outDir,outPlotImage);
  
}

#-----------------------------+
# COLOR PALETTE               |
#-----------------------------+
colRamp <- blue2red(numClass);
if (inverseColor == TRUE) {
  colRamp <- rev(colRamp);
}

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

  png(filename=outPlotImage, height=plotHeight, width=plotWidth, bg="white");
      
#  dev.copy (png, filename=outPlotImage,
#            height=plotHeight, width=plotWidth, bg="white" );
  
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
  # Add lines to the empirical cumulative distribution function
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

  dev.off();
  
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

#-----------------------------------------------------------+
# WRITE OUTPUT FILES                                        |
#-----------------------------------------------------------+

#-----------------------------+
# WRITE NA FILE FOR CYTOSCAPE |
#-----------------------------+

# WRITE THE HEADER
outHeader = (paste(baseName,"_", classMethod, "_", numClass, sep=""));
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
# - Convert R color codes to rgb codes for Cytoscape
# - Started adding pam .. too slow to work with existing dataset
# - Added option to not draw plot
# - Made plot labels variables
# - Making this a function
# - Filled space up in plot area better making use of 
#   of the par(mar) margin option
# - Added variables for x and y limits on plot
# - Added option to subset the data, example with mcl
#   cluster
#
# 05/07/2009
# - Added the ability to do reverse color ramp order
# - Adding command line arguments
#    --infile
#    --outfile
#    --class
#    --num-class
#    --num-breaks
#-----------------------------------------------------------+
# JUNKYARD
#-----------------------------------------------------------+
