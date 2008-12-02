## R random number generator


# Sample, testing with set variables
N <- 13588;    # Total sample size
n <- 100;      # Subsample size
out_file_path <- paste("outfile_",n,".txt", sep = "");
my_sample <- sample(1:N, n);
out_info <- array(n,c(n,2));
out_info[,2] <- my_sample;
write.table( out_info, file = out_file_path, quote = FALSE, row.names = FALSE, col.names = FALSE);

# Doing the above with a range
N <- 13588;    # Total sample size

for (i in 1:13) {
  n <- i*100;
  
  out_file_path <- paste("seq_sample_",n,".txt", sep = "");
  my_sample <- sample(1:N, n);
  out_info <- array(n,c(n,2));
  out_info[,2] <- my_sample;
  write.table( out_info, file = out_file_path, quote = FALSE, row.names = FALSE, col.names = FALSE);
}

# Doing the above with a shorter range
N <- 13588;    # Total sample size

for (i in 1:9) {
  n <- i*10;
  
  out_file_path <- paste("seq_sample_",n,".txt", sep = "");
  my_sample <- sample(1:N, n);
  out_info <- array(n,c(n,2));
  out_info[,2] <- my_sample;
  write.table( out_info, file = out_file_path, quote = FALSE, row.names = FALSE, col.names = FALSE);
}
