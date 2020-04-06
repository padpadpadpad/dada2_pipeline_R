#----------------------------------------------------------------------#
# install all the packages necessary to run the sequencing pipeline ####
#----------------------------------------------------------------------#

# install remotes if necessary
install.packages('remotes')

# install MicrobioUoE - a padpadpadpad R package
# https://github.com/padpadpadpad/MicrobioUoE
remotes::install_github('padpadpadpad/MicrobioUoE')

# some packages live on CRAN and some are only available through bioconductor

# we first need to check Bioconductor is installed and working correctly
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

# list packages to install ####
# cran
cran_packages  <-  c("ggplot2",
                     "dplyr",
                     'tidyr',
                     'magrittr',
                     'vegan',
                     "scales", 
                     "phangorn",
                     'forcats',
                     'quadprog')

# bioconductor
bioc_packages <- c("phyloseq", 
                   "genefilter", 
                   "impute", 
                   "DECIPHER",
                   'dada2')

# install all packages ####
MicrobioUoE::package_install_all(cran_packages = cran_packages, bioc_packages = bioc_packages)

# Huzzah all packages will be installed
