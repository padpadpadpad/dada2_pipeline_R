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
# possible values, "pseudo", TRUE, or FALSE
pool = 'pseudo'

#--------------------#
# 3. rename files ####
#--------------------#

# this step is key because we want to rename the files to be a shorter name that links nicely to the metadata file. Names of sequencing files from sequencing centres are often quite long and contain extra information on the lane of the sequencer, the date of the run, and the adapter sequence for that file
# consequently I would STRONGLY recommend renaming the files BEFORE going any further. I would also recommend creating metadata for your samples BEFORE sending them for sequencing, and then naming files 1,2,3,4,n to make this renaming easier. This is what we have implemented in this example.



# 