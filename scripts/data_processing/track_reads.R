#--------------------------------#
# merge sequences using dada2 ####
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

cat(paste("\ntrack_reads.R started at ", Sys.time()), file = progress_file, append = TRUE)

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

#--------------------------------------------------------------------#

out <- readRDS(paste(output_path, 'temp', 'filtered_read_numbers.rds', sep = '/'))
seqtab_all <- readRDS(paste(output_path, 'temp', 'seqtab_all.rds', sep = '/'))
seqtab <- readRDS(paste(output_path, 'temp', 'seqtab.rds', sep = '/'))
dadaFs <- readRDS(paste(output_path, 'temp', 'dadaFs.rds', sep = '/'))
mergers <- readRDS(paste(output_path, 'temp', 'mergers.rds', sep = '/'))

# track reads through pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab_all), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_namesF

saveRDS(track, paste(output_path, 'tracked_reads.rds', sep = '/'))

# plot out tracking of sample reads through stages ####
samps <- row.names(track)
track <- data.frame(track) %>%
  mutate(samps = samps) %>%
  gather(., 'stage', 'reads', c(input, filtered, denoised, merged, tabled, nonchim))

max <- max(log10(track$reads))

ggplot(track, aes(forcats::fct_relevel(stage, c('input', 'filtered', 'denoised', 'merged', 'tabled', 'nonchim')), reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2, height = 0), shape = 21, fill = 'white') +
  ylab('Number of reads') +
  xlab('Sequencing stage') +
  theme_bw(base_size = 16) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^0, 10^max),
                minor_breaks = NULL)

ggsave(file.path(plot_path, 'track_reads.pdf'), height = 5, width = 7)

cat(paste("\ntrack_reads.R ended at ", Sys.time()), file = progress_file, append = TRUE)
