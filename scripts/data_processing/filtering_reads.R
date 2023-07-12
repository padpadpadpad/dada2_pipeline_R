#-----------------#
# Filter reads ####
#-----------------#

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

cat(paste("\nfiltering_reads.R started at", Sys.time()), file = progress_file, append = TRUE)

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast|fq', full.names = TRUE))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1|r1", fns)]
fnRs <- fns[grepl("R2|r2", fns)]

# name the files in the list ####
# need to be the same as in the metadata file
sample_namesF <- get_sample_names(fnFs, meta$sample_id)
sample_namesR <- get_sample_names(fnRs, meta$sample_id)

# Trim and filter ####
filtFs <- file.path(filt_path, basename(fnFs)) 
filtRs <- file.path(filt_path, basename(fnRs))



# run filter parameters ####
# this can be based on the quality profiles in qual_plot_preFilt.pdf

trim_left <- c(0,0)
trunc_len=c(200,200)
max_n = 0
max_ee = Inf
trunc_q = 2

# set trimming parameters ####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     trimLeft = trim_left, 
                     truncLen=trunc_len,
                     maxN = max_n, 
                     maxEE = max_ee,
                     truncQ = trunc_q,
                     compress=TRUE, 
                     verbose=TRUE,
                     multithread=TRUE)

# check the number of reads post filtering
saveRDS(out, file = paste(output_path, 'temp', 'filtered_read_numbers.rds', sep ='/'))
# check the number of reads post filtering
write.csv(out, file = paste(output_path, 'temp', 'filtered_read_numbers.csv', sep ='/'))

cat(paste('\nfiltering parameters are:'), file = progress_file, append = TRUE)
cat(paste('\ntrimLeft =', trim_left), file = progress_file, append = TRUE)
cat(paste('\ntruncLen =', trunc_len), file = progress_file, append = TRUE)
cat(paste('\nmaxN =', max_n), file = progress_file, append = TRUE)
cat(paste('\nmaxEE =', max_ee), file = progress_file, append = TRUE)
cat(paste('\ntruncQ =', trunc_q), file = progress_file, append = TRUE)

# check quality of data ####
pdf(file.path(plot_path, 'qual_plot_preFilt.pdf'))
print(plotQualityProfile(fnFs, n = 2e6, aggregate = TRUE) +
        ggtitle('Fwd reads master quality profile'))
print(plotQualityProfile(fnRs, n = 2e6, aggregate = TRUE) +
        ggtitle('Rev reads master quality profile'))
dev.off()

pdf(file.path(plot_path, 'qual_plot_postFilt.pdf'))
print(plotQualityProfile(filtFs, n = 2e6, aggregate = TRUE) +
        ggtitle('Fwd reads master quality profile'))
print(plotQualityProfile(filtRs, n = 2e6, aggregate = TRUE) +
        ggtitle('Rev reads master quality profile'))
dev.off()

# add update to progress file
cat(paste('\nfiltering_reads.R ended at', Sys.time()), file = progress_file, append = TRUE)

