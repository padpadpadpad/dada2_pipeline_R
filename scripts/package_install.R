# install all the packages necessary to run everything in the sequencing pipeline

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

# github
github_packages <- c('benjjneb/dada2')

# bioconductor
bioc_packages <- c("phyloseq", 
                   "genefilter", 
                   "impute", 
                   "DECIPHER")

# install all packages
MicrobioUoE::package_install_all(cran_packages = cran_packages, github_packages = github_packages, bioc_packages = bioc_packages)

# will give a list of the failed packages

# bingo - re-run to check
MicrobioUoE::package_install_all(cran_packages = cran_packages, github_packages = github_packages, bioc_packages = bioc_packages)
# Huzzah
