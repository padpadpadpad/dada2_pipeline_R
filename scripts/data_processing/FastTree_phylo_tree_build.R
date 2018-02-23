# alternative phylogenetic tree with FastTree ####

# If the creation of a phylogenetic tree takes too long in R, you can attempt to use FastTree and bind the corresponding file with your existing sequencing run

# LINKS
# FastTree - http://www.microbesonline.org/fasttree/
# Borrows from a similar approach on this phyloseq GitHub Issue - https://github.com/padpadpadpad/bact_phage_temp

# load packages ####
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(dada2)

# where is your FastTree saved
FastTree_path <- '~/Desktop/FastTree/FastTree'

# find alignment fasta file ####
alignment <- list.files('data', pattern = 'alignment', full.names = TRUE, recursive = TRUE)

# stop if alignment is empty
if(length(alignment) == 0){
  stop("There is no alignment file. Please make sure there is an alignment file in one of the output folders of your analysis runs.")
}

# Choose alignment file if there are multiple files
alignment <- alignment[1]

# run using FastTree ####

# create FastTree call
FastTree_command <- paste(FastTree_path, ' -nt ', alignment, ' > ', dirname(alignment), '/tree', sep = '')

# run FastTree - this will not work on Windows!!!
system(FastTree_command)

# find files needed for phyloseq object 
seqtab <- list.files(dirname(alignment), pattern = 'seqtab', full.names = TRUE)
taxtab <- list.files(dirname(alignment), pattern = 'taxtab', full.names = TRUE)
tree_path <- list.files(dirname(alignment), pattern = 'tree', full.names = TRUE)

# load in previously run files
seqtab <- readRDS(seqtab)
meta <- read.csv('data/metadata_example.csv')
taxtab <- readRDS(taxtab)
tree <- read_tree(tree_path)

# change row names
rownames(meta) <- meta$SampleID

# make phyloseq object
ps <- phyloseq(tax_table(taxtab), 
              sample_data(meta),
              otu_table(seqtab, taxa_are_rows = FALSE), 
              phy_tree(tree))

# create new ps file path
ps_new <- paste(dirname(tree_path), '/', basename(dirname(tree_path)), '_ps.rds', sep = '')

# overwrite previous phyloseq object with new tree
saveRDS(ps, ps_new)
