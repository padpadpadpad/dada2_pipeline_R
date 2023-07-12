#--------------------------------#
# infer sequences using dada2 ####
#--------------------------------#

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# read in config file
config <- read.csv(snakemake@input[["config_file"]], stringsAsFactors = FALSE)
# config <- read.csv('config.csv', stringsAsFactors = FALSE)

# read in extra functions
source(snakemake@input[["extra_functions"]])
#source('scripts/functions.R')

# setup paths and packages and error if previously estimated ####
# this will automatically set up the environment to run the analysis
dada2_pipeline_setup(packages = trimws(unlist(strsplit(subset(config, input == 'packages')$value, ','))),
                     raw_path = subset(config, input == 'raw_path')$value,
                     filt_path = subset(config, input == 'filt_path')$value,
                     plot_path = subset(config, input == 'plot_path')$value,
                     output_path = subset(config, input == 'output_path')$value,
                     progress_path = subset(config, input == 'progress_path')$value,
                     ref_fasta = subset(config, input == 'ref_fasta')$value,
                     ref_fasta_spp = subset(config, input == 'ref_fasta_spp')$value,
                     meta_data = subset(config, input == 'meta_data')$value,
                     run = subset(config, input == 'run')$value,
                     pool = subset(config, input == 'pool')$value)

cat(paste("\ninfer_sequences.R started at ", Sys.time()), file = progress_file, append = TRUE)

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast|fq', full.names = TRUE, recursive = T))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1|r1", fns)]
fnRs <- fns[grepl("R2|r1", fns)]

# name the files in the list ####
# need to be the same as in the metadata file
sample_namesF <- get_sample_names(fnFs, meta$sample_id)
sample_namesR <- get_sample_names(fnRs, meta$sample_id)

# get filt path ####
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# name files
names(filtFs) <- sample_namesF
names(filtRs) <- sample_namesR

# load in error rates
dd_learnF <- readRDS(paste(output_path, 'temp', 'fwd_error.rds', sep = '/'))
dd_learnR <- readRDS(paste(output_path, 'temp', 'rev_error.rds', sep = '/'))

# infer sequence data ####
dadaFs <- dada(filtFs, err = dd_learnF, pool = TRUE, multithread = TRUE)
cat(paste('\nForward sequences inferred', Sys.time()), file = progress_file, append = TRUE)
dadaRs <- dada(filtRs, err = dd_learnR, pool = TRUE, multithread = TRUE)
cat(paste('\nReverse sequences inferred', Sys.time()), file = progress_file, append = TRUE)

# save out the dada class objects
saveRDS(dadaFs, paste(output_path, 'temp', 'dadaFs.rds', sep = '/'))
saveRDS(dadaRs, paste(output_path, 'temp', 'dadaRs.rds', sep = '/'))

cat(paste("\ninfer_sequences.R ended at ", Sys.time()), file = progress_file, append = TRUE)
