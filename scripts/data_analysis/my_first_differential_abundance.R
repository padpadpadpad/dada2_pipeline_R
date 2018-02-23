# look at differential abundances of OTUs
# There are a couple of approaches that could be taken here
# 1. From Reese et al. 2017 Ecology and Evolution - take out the relative abundance of taxa (pooled at the Phylum level in this paper) and do frequentist statistics on this, lmer
# 2. Using Deseq2 to look at abundance changes based on treatment (would need to treat it as a block design)

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(lme4)
library(MicrobioUoE) # if not installed run devtools::install_github('padpadpadpad/MicrobioUoE')

# set seed
set.seed(42)

# figure path
path_fig <- 'plots'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('data/output/20171130_11h01m/ps_prevalence_filtered.rds')

# replace metadata with new metadata
meta_new <- read.csv('data/metadata.csv', stringsAsFactors = FALSE) %>%
  mutate(., time = as.character(time))
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# replace metadata
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at unique phyla
get_taxa_unique(ps, "Phylum")

# delete some samples based on taxonomy
ps <- subset_taxa(ps, Phylum != 'Cyanobacteria/Chloroplast')

# agglomerate results by phyla 
tax_group <- tax_glom(ps, taxrank = "Phylum")

# make counts proportions
tax_prop <- transform_sample_counts(tax_group, function(x){x / sum(x)})

# get data out of phyloseq object
d_group <- psmelt(tax_prop) %>%
  janitor::clean_names()

# filter things under 2%
# d_group <- filter(d_group, abundance > 0.02)

# add dummy data
# d_group2 <- tidyr::complete(., phylum, treat, time, fill = list(abundance = 0))

# plot
ggplot(d_group, aes(time, abundance, fill = treat, col = treat)) +
  geom_pretty_boxplot() +
  geom_point(fill = 'white', shape = 21, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1), size = 0.5) +
  facet_wrap(~ phylum, scales = 'free_y') +
  scale_color_grey() +
  scale_fill_grey() +
  theme_bw()

# save plot, other ways are available
ggsave(file.path(path_fig, 'diff_abund.pdf'), last_plot(), height = 9, width = 13)

# use a linear mixed effect model
mod <- lm(abundance ~ time*phylum, d_group)

# look at differences per phylum using lsmeans
lsmeans::lsmeans(mod, ~ time|phylum)

# no real difference in proportion across samples - quite like this approach though, easy to understand

# Have not implemented Deseq2 yet!!!