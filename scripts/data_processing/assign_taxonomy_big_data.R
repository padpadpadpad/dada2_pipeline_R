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
#config <- read.csv('~/dada2_pipeline_myxorpoB/config.csv', stringsAsFactors = FALSE)

# read in extra functions
source(snakemake@input[["extra_functions"]])
# source('~/dada2_pipeline_myxorpoB/scripts/functions.R')

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

cat(paste("\nassign_taxonomy_big_data.R started at ", Sys.time()), file = progress_file, append = TRUE)

#--------------------------------------------------------------------#

seqtab <- readRDS(paste(output_path, 'temp', 'seqtab.rds', sep = '/'))

# assign taxonomy in batches of 50,000
to_split <- seq(1, ncol(seqtab), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, ncol(seqtab))

taxtab = NULL

for(i in 1:length(to_split)){
  seqtab2 <- seqtab[, to_split[i]:to_split2[i]]
  taxtab2 <- assignTaxonomy(seqtab2, refFasta = ref_fasta, multithread = TRUE)
  if(ref_fasta_spp != 'NULL'){taxtab2 <- addSpecies(taxtab2, refFasta = ref_fasta_spp, verbose = TRUE)}
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}
# save files in case phylogeny does not run
saveRDS(taxtab, paste(output_path, 'temp', 'taxtab.rds', sep = '/'))

# subset meta for just the samples present

# SampleID needs to match sample_namesF
rownames(meta) <- meta$sample_id
ps <- phyloseq(tax_table(taxtab), 
               sample_data(meta),
               otu_table(seqtab, taxa_are_rows = FALSE))

# save a phyloseq object without the phylogeny
saveRDS(ps, paste(output_path, 'ps_no_tree.rds', sep = '/'))
save(ps, file = paste(output_path, 'ps_no_tree.Rdata', sep = '/'))

cat(paste('\nEnd of processing without construction of phylogeny', Sys.time()), file = progress_file, append = TRUE)
cat(paste("\nassign_taxonomy_big_data.R ended at ", Sys.time()), file = progress_file, append = TRUE)

