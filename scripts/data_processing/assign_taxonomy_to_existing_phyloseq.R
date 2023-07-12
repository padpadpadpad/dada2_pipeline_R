#----------------------------------------#
# run assign taxonomy on renamed otus ####
#----------------------------------------#

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# set minBoot
min_boot = 50

# read in config file
config <- read.csv(snakemake@input[["config_file"]], stringsAsFactors = FALSE)

# load in necessary packages
packages = trimws(unlist(strsplit(subset(config, input == 'packages')$value, ',')))
lapply(packages, library, character.only = TRUE)

# assign ref_fasta and ref_fasta_spp
ref_fasta = subset(config, input == 'ref_fasta')$value
ref_fasta_spp = subset(config, input == 'ref_fasta_spp')$value

# read in spreadsheet of otu names and sequences
d <- read.csv('~/dada2_pipeline_myxorpoB/data/output/run_myxo_gtdbr202/myxo_seq_name_conversion.csv')
head(d)

seqs <- d$seq

# assign taxonomy as not over 50,000 OTUs
taxtab <- assignTaxonomy(seqs, refFasta = ref_fasta, multithread = TRUE)
if(ref_fasta_spp != 'NULL'){taxtab <- addSpecies(taxtab, refFasta = ref_fasta_spp, verbose = TRUE)}



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

