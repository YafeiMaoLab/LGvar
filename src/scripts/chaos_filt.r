options(warn = -1)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  divergence_threshold <- as.numeric(args[6]) 

  filt2saffire <- function(pos, name) {
    if (nrow(pos) == 0) {
      cat("Warning: Dataframe is empty for", name, "\n")
      return(NULL)
    }
    pos_copy <- copy(pos)
    col_names <- colnames(pos_copy)
    col_names[c(1, 3, 4, 8, 9, 6, 5)] <- c("query_chr", "query_start", "query_end", "ref_start", "ref_end", "ref_chr", "orient")
    colnames(pos_copy) <- col_names
    
    chrpc <- fread(args[1])
    chrpc <- distinct(chrpc)
    colnames(chrpc) <- c("ref_chr", "ref_len", "query_chr", "query_len")
    
    saffire <- merge(pos_copy, chrpc, by = c("ref_chr", "query_chr"))
    saffire <- saffire[, c("ref_chr", "ref_start", "ref_end", "ref_len", "orient", "query_chr", "query_start", "query_end", "query_len")]
    
    saffire_name <- c("#reference_name", "reference_start", "reference_end", "reference_length", "strand", 
                      "query_name", "query_start", "query_end", "query_length", "perID_by_matches", 
                      "perID_by_events", "perID_by_all", "matches", "mismatches", "deletion_events", 
                      "insertion_events", "deletions", "insertions")
    
    saffire[, c("perID_by_matches", "perID_by_events", "perID_by_all", "matches", "mismatches", 
                "deletion_events", "insertion_events", "deletions", "insertions")] <- 0
    colnames(saffire) <- saffire_name
    saffire$matches <- 100
    saffire$mismatches <- 100
    
    write.table(saffire, file = paste(args[2], name, sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
  }

  pos <- fread(args[3], fill = TRUE)
  cat("init", dim(pos)[1], "\n")

  filt2saffire(pos, "original.saffire")

  chrnames <- unique(pos$V6)
  numeric_part <- as.numeric(gsub("\\D", "", chrnames))
  sorted_chrnames <- chrnames[order(numeric_part)]
  
  filtered_list <- list()
  
  for (chrid in sorted_chrnames) {
    pos.chr <- pos[pos$V6 == chrid, ]
    print(chrid)

    v21_numeric <- as.numeric(gsub("[a-zA-Z:]", "", pos.chr$V21))
    pos.chr <- pos.chr[!is.na(v21_numeric) & v21_numeric <= divergence_threshold, ]
    
    filtered_list[[chrid]] <- pos.chr
  }

  final <- do.call(rbind, filtered_list)

  cat("Removed：", dim(final)[1], "\n") 
  filt2saffire(final, "final.saffire")
  
  write.table(final, file = args[4], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  print("Complete!!!!")
})
