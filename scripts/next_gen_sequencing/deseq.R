#!/usr/bin/env Rscript --vanilla

## An executable script for running DESeq2 in an automated fashion with a simple
## argument parser. Expects output from HTSeq counter. Coule be used for mRNA or sRNA.

library(DESeq2)

#### ---- Get the command line arguments ---- ####

args <- commandArgs(trailingOnly = TRUE)

# Throw an error if there are more args than expected
if (length(args) > 6){
  stop(sprintf("Too many arguments given. Only --input-files and --sample-names are accepted.
        Files and names must be comma-separated without spaces. The bad arguments given were 

        %s

        Did you add a space in a list?", args[5:length(args)]))
} else if (length(args) < 6){
  stop(sprintf("Not enough arguments given. Only parsed %d arguments.", length(args)))
} else if (length(args) == 0){
  stop("No arguments given. The following arguments are accepted:
       
       --input-files <file1,file2,file3>
              A comma separated (no spaces) list of count files for the run.
      
       --sample-names <name1,name2,name3>
              A comma separated (no spaces) list of sample names formatted as

              sample_replicate_N 

              where:
                    sample = condition or group to compare
                    N = the replicate number
      
       --outfile-prefix <outfile>
              Name of the output files to write. These will be created:
                  1. Normalized count table of all samples
                  2. Differential gene expression table per comparison

       ")
}

# Grab all the arg information
arg_pos <- 1
while (arg_pos <= length(args)){
  
  if (args[arg_pos] == '--input-files'){
    files <- unlist(strsplit(args[arg_pos + 1], ','))
  } else if (args[arg_pos] == '--sample-names'){
    samples <- unlist(strsplit(args[arg_pos + 1], ','))
  } else if (args[arg_pos] == '--outfile-prefix'){
    out_pref <- args[arg_pos + 1]
  } else if (!(args[arg_pos - 1] %in% c('--input-files', '--sample-names', '--outfile-prefix'))){
    stop(sprintf("This argument %s is not accepted. Did you accidentally include a space in your file or
         sample names?", args[arg_pos]))
  }
  arg_pos <- arg_pos + 1
}

#### ---- Set up the parameters ---- ####

# Create a data frame matching the file, name, condition
condition <- rep('none', length(samples))
for (i in 1:length(samples)){
  sample <- as.character(samples[i])
  condition[i] <- strsplit(sample, '_replicate_')[[1]][1]
}

# Create the count table (need to make fractional counts ints for DESeq2)
for (f in files){
  cat(as.character(f))
}

# Create the deseqdataset
sample_table <- data.frame(sampleName=samples, fileName=files, condition=condition)
deseq_ds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = "/", design = ~ condition)

#### ---- Run DESeq2 & write outputs ---- ####

# Create the DESeq Object
deseq_run <- DESeq(deseq_ds)

# Get the normalized counts
deseq_counts <- counts(deseq_run, normalized=TRUE)
write.csv(deseq_counts, paste(out_pref, "norm_counts.csv", sep="_"))

# Create & retrieve all possible comparisons
all_comparisons <- t(combn(unique(condition), 2))
for (i in 1:nrow(all_comparisons)){
  comparison <- all_comparisons[i,] 
  deseq_res <- results(deseq_run, c("condition", comparison[1], comparison[2]))
  deseq_res <- deseq_res[order(deseq_res$padj),]
  write.csv(deseq_res, paste(out_pref, comparison[1], comparison[2], "deseq_table.csv", sep="_"))
}

