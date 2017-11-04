# install all the packages necessary to run everything in the sequencing pipeline ####

# install devtools if necessary
install.packages('devtools')

# install MicrobioUoE
devtools::install_github('padpadpadpad/MicrobioUoE')

# list packages ####
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

# install all packages
MicrobioUoE::package_install_all(cran_packages = cran_packages, bioc_packages = bioc_packages)

# Huzzah all packages will be installed

# to try update bioconductor to newest version
# This will update all packages
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
