# my first clustering ####

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vegan)
library(viridis)
library(MicrobioUoE)
library(ggvegan)
# need to also install janitor from CRAN

# functions to use ####
# code stolen from phyloseq website
get_top_taxa <- function(ps, tax_rank, to_keep){
  temp <- tapply(phyloseq::taxa_sums(ps), phyloseq::tax_table(ps)[, tax_rank], sum, na.rm = TRUE)
  temp2 <-  names(sort(temp, TRUE))[1:to_keep]
  return(temp2)
}

# get distance from 00
dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}

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

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 50947. Woof.

# transform counts to relative abundances for ordination ####
ps_prop <- transform_sample_counts(ps, function(x){x / sum(x)})

# to change the ordination - see https://joey711.github.io/phyloseq/plot_ordination-examples.html and https://joey711.github.io/phyloseq/distance.html

#####################
# 1. Bray-Curtis dissimilarity ####

# do desired ordination. Here uses Bray-Curtis so use absolute abundances
ord_bray <- ordinate(ps, method = 'MDS', distance = 'bray')

evals <- ord_bray$values$Eigenvalues

# plot
plot_ordination(ps, ord_bray, col = 'time', shape = 'treat') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_point(size = 2) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('PCoA plot based on weighted Bray-Curtis distances')
# plot all on one panel by removing the facet_wrap command

# save plot, other ways are available
ggsave(file.path(path_fig, 'bray_ordination.pdf'), last_plot(), height = 6, width = 12)

# get the distance matrix out of the data
ps_bray <- phyloseq::distance(ps, method = 'bray')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp$treat<- as.factor(d_samp$treat)
d_samp$time <- as.factor(d_samp$time)

# run an Adonis test
mod1 <- vegan::adonis(ps_wunifrac ~ treat*time, data = d_samp)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_bray, d_samp, 'time', n_perm = 9999)
# loads of comparisons. Significant differences will be determined by p value correction.

################################################################
# 2. Do NMDS on Bray curtis to plot OTUs as well as effects ####
ord_nmds <- ordinate(ps, method = 'NMDS', distance = 'bray', try = 100)

plot_ordination(ps, ord_nmds, type = 'split', col = 'Phylum', shape = 'time')
# so many OTUs

# agglomerate at phylum level
tax_group <- tax_glom(ps, taxrank = "Phylum")
tax_prop <- transform_sample_counts(tax_group, function(x){x / sum(x)})

# get the distance matrix out of the data
ps_bray <- phyloseq::distance(tax_group, method = 'bray')
ps_table <- psmelt(tax_group) %>%
  select(SampleID, treat, time, Phylum, Abundance) %>%
  spread(., Phylum, Abundance) %>%
  janitor::clean_names() %>%
  select(., -cyanobacteria_chloroplast)

# abundance table 
d_abund <- select(ps_table, -c(treat, time)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(., 'sampleid')

# sample info table
d_samp <- select(ps_table, sampleid, treat, time)

# do nmds
ord_nmds <- metaMDS(d_abund, distance = 'bray')

# get data from nmds
d_nmds <- fortify(ord_nmds) %>%
  janitor::clean_names()

# wrangle sites
d_sites <- filter(d_nmds, score== 'sites') %>%
  rename(sampleid = label) %>%
  merge(., d_samp, by = 'sampleid')

# wrangle species
d_species <- filter(d_nmds, score == 'species') %>%
  mutate(dist = dist_from_00(nmds1, nmds2))

# plot
ggplot() +  
  geom_text(aes(nmds1, nmds2, label = label, hjust = 0.5*(1 - sign(nmds1)), vjust = 0.5*(1-sign(nmds2)), alpha = dist), col = 'red4', d_species) +
  geom_segment(aes(x = 0, y = 0, yend = nmds2, xend = nmds1, group = score, alpha = dist), d_species, arrow = arrow(length = unit(0.01, "npc")), col = 'red4') +
  geom_point(aes(nmds1, nmds2, col = time, shape = treat), d_sites) +
  scale_x_continuous(expand = c(.12, .12)) +
  scale_y_continuous(expand = c(.12, .12)) +
  scale_alpha(guide = FALSE) +
  ggtitle('NMDS of Bray-Curtis dissimilarity grouped at the Phylum level',
          subtitle = 'Loadings of each phylum are added as arrows') +
  theme_bw() +
  coord_equal()

# save plot, other ways are available
ggsave(file.path(path_fig, 'braycurtis_nmds.pdf'), last_plot(), height = 9, width = 11)

# test against the automatic method - LOOKS THE SAME - SUCCESS
ord_nmds2 <- ordinate(tax_group, method = 'NMDS', distance = 'bray')
plot_ordination(tax_prop, ord_nmds2, type = 'split', col = 'Phylum', shape = 'time')

#####################################################
# look at proportion of pseudomonas through time ####
# get data
d_ps <- psmelt(ps_prop)
unique(d_ps$Species)

d_Pseudomonas <- filter(d_ps, Genus == 'Pseudomonas') %>%
  select(., -OTU) %>%
  unite(., spp, c(Genus, Species), sep = '_', remove = FALSE)

d_Pseud_prop <- group_by(d_Pseudomonas, SampleID, rep, treat, time) %>%
  summarise(., prop = sum(Abundance)) %>%
  ungroup()

ggplot(d_Pseud_prop, aes(treat, prop, col = treat, fill = treat)) +
  geom_pretty_boxplot() +
  facet_wrap(~ time) +
  geom_point(position = position_jitter(height = 0, width = 0.1), fill = 'white', shape = 21)


################################################
# look at common taxa by treatment and time ####
# summarise phyloseq object at the Phylum level
tax_group <- tax_glom(ps, taxrank = "Phylum")

# filter for the 10 most common genus ####
tax_group_filt <- prune_taxa((tax_table(tax_group)[, "Phylum"] %in% get_top_taxa(tax_group, 'Phylum', 10)), tax_group)

# convert counts to proportions
tax_prop <- transform_sample_counts(tax_group, function(x){x / sum(x)})
tax_prop_filt <- transform_sample_counts(tax_group_filt, function(x){x / sum(x)})

# plot bar plot
plot_bar(tax_prop_filt, fill = "Phylum") +
  facet_wrap(~ time, scale = 'free_x')

# save plot, other ways are available
ggsave(file.path(path_fig, 'plot_bar1.pdf'), last_plot(), height = 7, width = 12)

# try to create a prettier bar plot

# get data
d_glom <- psmelt(tax_prop_filt)

# group by treatments
d_glom_group <- group_by(d_glom, treat, time) %>%
  do(., data.frame(prop = .$Abundance/sum(.$Abundance), Phylum = .$Phylum)) %>%
  ungroup()

# plot these
ggplot(d_glom_group, aes(interaction(treat, time), prop, fill = Phylum, col = Phylum)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Relative abundance of 10 most abundant genus groups') +
  ylab('Proportion') +
  xlab('Treatment by Time')

ggsave(file.path(path_fig, 'plot_bar2.pdf'), last_plot(), height = 5, width = 8)




