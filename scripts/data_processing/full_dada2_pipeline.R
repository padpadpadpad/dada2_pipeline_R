#-----------------------------------------#
# dada2 amplicon sequencing workflow in R #
#-----------------------------------------#

#--------------------------#
# What this script does ####
#--------------------------#

# 1. installs the necessary packages for the workflow
# 2. set a bunch of paths for our different folders
# 3. renames raw files for easy parsing through the pipeline and so they match the sample names in the metadata file

#----------------------------------#
# 1. install necessary packages ####
#----------------------------------#

# first please check you have at least R 4.0.0 installed

# firstly check that BiocInstall and librarian are installed
if(length(find.package('BiocManager', quiet = TRUE)) == 0){install.packages('BiocManager')}
if(length(find.package('librarian', quiet = TRUE)) == 0){install.packages('librarian')}

# install and load packages used in the workflow using librarian::shelf
librarian::shelf(dada2, phyloseq, DECIPHER, phangorn, tidyverse)

#-----------------#
# 2. set paths ####
#-----------------#

# set run
run = '1'

# set path to raw files
raw_path = "data/raw_fastq"

# set path to filtered fastq files
filt_path = "data/filtered_fastq"

# set path to where to save plots
plot_path = "plots"

# set path to where to put the output
output_path = "data/output" 

# set path for where to store progress
progress_path = "data/progress" 

# set path for where reference databases are
ref_fasta = "data/ref_trainsets/rdp_train_set_16.fa"
ref_fasta_spp = "data/ref_trainsets/rdp_species_assignment_16.fa"

# set path for where metadata is
meta_data = "data/metadata.csv"

# specify whether to pool ASVs across samples - this will allow the finding of a few more rare variants, but will come at a computation cost
# find out more here: https://benjjneb.github.io/dada2/pool.html
# possible values, "pseudo", TRUE, or FALSE
pool = 'pseudo'

#--------------------#
# 3. rename files ####
#--------------------#

# this step is key because we want to rename the files to be a shorter name that links nicely to the metadata file. Names of sequencing files from sequencing centres are often quite long and contain extra information on the lane of the sequencer, the date of the run, and the adapter sequence for that file
# consequently I would STRONGLY recommend renaming the files BEFORE going any further. I would also recommend creating metadata for your samples BEFORE sending them for sequencing, and then naming files 1, 2, 3, 4, n to make renaming easier. This is what we have implemented in this example. With a sample number column in the metadata

# create dataframe of raw data files
d_rawfiles <- data.frame(full_name = list.files(raw_path, full.names = TRUE, recursive = TRUE, pattern = '.fastq.gz'))
d_rawfiles$full_name[1]

# we can then create some new columns to create the new names for the files, which we want to be sample_1_R1.fastq.gz sample_1_R2.fastq.gz etc
# I do not really know regex well, so I use chatGPT to help me write code to extract the bits of the file name that I need
# here I use strsplit() to keep the relevant parts between underscores in the filename
d_rawfiles <- mutate(d_rawfiles, filename = basename(tools::file_path_sans_ext(full_name)),
                     directory = dirname(full_name),
                     sample_number = map(filename, \(x){strsplit(x, '_') %>% unlist() %>% .[3]}),
                     fwd_rev = map(filename, \(x){strsplit(x, '_') %>% unlist() %>% .[4]}))

# we can then create a new filename
d_rawfiles <- mutate(d_rawfiles, new_filename = paste('sample_', sample_number, '_', fwd_rev, '.fastq.gz', sep = ''))

# we can then run a for loop to rename the files
# i will put the renamed files in a new folder here, but you can just overwrite them likely
dir.create(file.path(raw_path, 'renamed'))
for(i in 1:nrow(d_rawfiles)){
  # copy each file
  file.copy(d_rawfiles$full_name[i], paste(d_rawfiles$directory[i], 'renamed', d_rawfiles$new_filename[i], sep = '/'))
}
# this has copied each file, you can now delete the old named files if you so wish

# we will save out d_rawfiles so we can go back to the old names if need be
select(d_rawfiles, full_name, original_filename = filename, new_filename) %>%
  write.csv('data/filename_changes.csv', row.names = FALSE)

