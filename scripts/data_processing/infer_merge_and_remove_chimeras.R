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

cat(paste("\ninfer_merge_and_remove_chimeras.R started at ", Sys.time()), file = progress_file, append = TRUE)

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast|fq', full.names = TRUE))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1|r1", fns)]
fnRs <- fns[grepl("R2|r2", fns)]

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

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample_namesF))
names(mergers) <- sample_namesF
dadaFs <- vector("list", length(sample_namesF))
names(dadaFs) <- sample_namesF
dadaRs <- vector("list", length(sample_namesF))
names(dadaRs) <- sample_namesF

for(sam in sample_namesF) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=dd_learnF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=dd_learnR, multithread=TRUE)
  
  dadaFs[[sam]] <- ddF
  dadaRs[[sam]] <- ddR
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR); rm(ddF); rm(ddR)

# save out objects
saveRDS(mergers, paste(output_path, 'temp', 'mergers.rds', sep = '/'))
saveRDS(dadaFs, paste(output_path, 'temp', 'dadaFs.rds', sep = '/'))
saveRDS(dadaRs, paste(output_path, 'temp', 'dadaRs.rds', sep = '/'))

rm(dadaFs); rm(dadaRs)

# construct sequence table ####
seqtab_all <- makeSequenceTable(mergers)
rm(mergers)

saveRDS(seqtab_all, paste(output_path, 'temp', 'seqtab_all.rds', sep = '/'))

# remove chimeric sequences ####
seqtab <- removeBimeraDenovo(seqtab_all, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab, paste(output_path, 'temp', 'seqtab.rds', sep = '/'))

cat(paste("\ninfer_merge_and_remove_chimeras.R ended at ", Sys.time()), file = progress_file, append = TRUE)
