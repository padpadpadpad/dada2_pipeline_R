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

cat(paste("\nasv_to_otus.R started at ", Sys.time()), file = progress_file, append = TRUE)

#--------------------------------------------------------------------#

percent_similarity <- as.numeric(subset(config, input == 'otu_percent_similarity')$value)
cut_off = (100-percent_similarity) / 100

# load in alignment file
alignment <- Biostrings::readDNAStringSet(file.path(output_path, paste('temp/', subset(config, input == 'alignment_file')$value, sep = '')))

ps <- readRDS(file.path(output_path, subset(config, input == 'ps_file')$value))

# calculate distance matrix for each sequence
dist_matrix <- readRDS(file.path(output_path, paste('temp/', subset(config, input == 'distance_matrix')$value, sep = '')))

clusters <- DECIPHER::IdClusters(
  dist_matrix, 
  method = "UPGMA",
  cutoff = cut_off,
  processors = NULL
)

ps0 <- merge_taxa_vec(
  ps, 
  group = clusters$cluster,
)

saveRDS(ps0, paste(output_path, '/ps_otu_', percent_similarity, 'percent.rds', sep = ''))

cat(paste("\nasv_to_otus.R.R ended at ", Sys.time()), file = progress_file, append = TRUE)
