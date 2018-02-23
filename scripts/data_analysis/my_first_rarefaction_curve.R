# plotting rarefaction curves

# streamlined script for analysis thus far ####
rm(list = ls())

# load packages ####
library(phyloseq)
library(vegan)
# if not installed, install mctoolsr run devtools::install_github('leffj/mctoolsr')

# set seed
set.seed(42)

# figure path
path_fig <- 'plots'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('data/output/20171130_11h01m/ps_prevalence_filtered.rds')

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 50947. Woof.

# not going to rarefy those samples yet

# can plot rarefaction curves
# code from https://www.fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
# written by Gavin Simpson himself (an author of vegan)

# check rarefaction curves ####
ps_otu_table <- data.frame(otu_table(ps))
raremax <- min(rowSums(ps_otu_table))
col <- c('black', 'blue', 'yellow', 'red', 'orange', 'grey', 'hotpink', 'purple', 'green')
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, label = FALSE))

# save plot
pdf(file.path(path_fig, 'rarefaction_curve.pdf'), width = 9, height = 7)
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
out <- with(pars,
            rarecurve(ps_otu_table, step = 1000, sample = raremax, col = col,
                      lty = lty, label = FALSE))
dev.off()