#----------------------#
# CHECK SETUP IS OK ####
#----------------------#


# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# read in config file
config <- read.csv(snakemake@input[["config_file"]], stringsAsFactors = FALSE)

# read in extra functions
source(snakemake@input[["extra_functions"]])
#source('scripts/extra_functions.R')

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

cat(paste("\ncheck_setup.R started at ", Sys.time()), file = progress_file, append = TRUE)

cat(paste('\nThis run was done using pool =', pool, 'during sample inference using dada2'), file = progress_file, append = TRUE)

# list files ####
fns <- sort(list.files(raw_path, pattern = 'fast|fq', full.names = TRUE, recursive = T))

# sort files for forward and reverse sequences ####
fnFs <- fns[grepl("R1|r1", fns)]
fnRs <- fns[grepl("R2|r1", fns)]

# name the files in the list ####
# need to be the same as in the metadata file
sample_namesF <- get_sample_names(fnFs, meta$sample_id)
sample_namesR <- get_sample_names(fnRs, meta$sample_id)

####################
##### CHECKS #######
####################

# check Fwd and Rev files are in the same order
if(!identical(sample_namesF, sample_namesR)) cat(paste('\nForward and reverse files do not match.'), file = progress_file, append = TRUE)

# check Fwd and Rev files are in the SampleID column of the meta data
if(all(sample_namesF %in% meta$sample_id) == FALSE) cat("\nForward and reverse file names are not present in metadata column Sample_id", file = progress_file, append = TRUE)

# check that ref_trainsets are present
if(file.exists(ref_fasta) == FALSE) cat("\nreference trainset, ref_fasta, is not in data/ref_trainset. Please make sure file has correct spelling and is present in the correct folder.", file = progress_file, append = TRUE)
if(!is.null(ref_fasta_spp)){
  if(file.exists(ref_fasta_spp) == FALSE) cat("\nreference trainset for species assignment, ref_fasta_spp, is not in data/ref_trainset. Please make sure file has correct spelling and is present in the correct folder.", file = progress_file, append = TRUE)
}

cat(paste('\ncheck_setup.R ended at' , Sys.time()), file = progress_file, append = TRUE)

write.csv(config, paste(file.path(output_path, 'config.csv')), row.names = FALSE, quote = FALSE)