# script for running dada2 pipeline on a server

# first install mamba as its faster for installing than conda
conda install -c conda-forge mamba

# create an environment
mamba create -n dada2_pipeline 

# open the environment
conda activate dada2_pipeline

# install needed dependencies for the pipeline
mamba install -c conda-forge snakemake
mamba install -c bioconda fasttree
mamba install -c r r
mamba install r-tidyverse r-phangorn
mamba install -c bioconda bioconductor-DECIPHER bioconductor-phyloseq bioconductor-dada2

# set working directory
wd=~/google_drive/work_googledrive/dada2_pipeline
cd $wd

snakemake --version

# available steps
# check_setup

# run each stage separately
# progress can be checked in the data/progress/run_*_progress.txt
snakemake --cores 2 -R --until check_setup
snakemake --cores 2 -R --until filter_reads
snakemake --cores 2 -R --until estimate_errors
snakemake --cores 2 -R --until infer_sequences
snakemake --cores 2 -R --until merge_sequences
snakemake --cores 2 -R --until make_sequence_table
snakemake --cores 2 -R --until remove_chimeras
snakemake --cores 2 -R --until assign_taxonomy
snakemake --cores 2 -R --until align_sequences
snakemake --cores 2 -R --until make_tree
snakemake --cores 2 -R --until make_phyloseq
snakemake --cores 2 -R --until track_reads

mamba search -c conda-forge r

mamba search -c bioconda -c conda-forge bioconductor-dada2
mamba search -c bioconda -c conda-forge liblapack
mamba install -c bioconda bioconductor-rhdf5filters
mamba install -c bioconda bioconductor-rhdf5
conda search -c bioconda -c conda-forge bioconductor-Rsamtools

# install necessary R packages
R

# from CRAN
install.packages('dplyr')
install.packages('tidyr')
install.packages("ggplot2")
install.packages('phangorn')
install.packages('BiocManager')

# from Bioconductor
BiocManager::install(version = "3.12")
BiocManager::install('DECIPHER')
BiocManager::install('phyloseq')
BiocManager::install('dada2')

quit()

git clone https://git@git.bioconductor.org/packages/IRanges

R CMD build IRanges
rm -r IRanges/vignettes
R CMD INSTALL rhdf5filters
