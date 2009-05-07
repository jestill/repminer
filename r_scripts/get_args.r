# The GetArg subfunction allows for arguments to be passed
# to the R script from the command line
GetArg <-function(strSwitch) {
  # The ComArgs variable is used to fetch the data following comargs
  ComArgs = commandArgs()
    # The number of arguments is equal to the length of the commands that occur after the args
  numArgs = (length(ComArgs))               
  for (i in 1:numArgs) {
        # Find the strSwitch string that was passed to the function
    if (ComArgs[i] == strSwitch) {
          # The variable variable is assumed to be the next string occuring directly after the flag string
      myVar <- ComArgs[i+1]             
    }
  }
  
  GetArg <- myVar
  
}

