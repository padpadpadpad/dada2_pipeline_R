# create file for metadata
# rename raw files so they match with the metadata

library(readxl)
library(tidyverse)

# load in file
d <- read.csv('~/dada2_pipeline_myxo16S/data/metadata_complete.csv', stringsAsFactors = FALSE)

head(d)

d <- mutate(d, sample = paste('s', id, sep = ''))

# grab out IDs
seq_files <- tibble(file = list.files('~/dada2_pipeline_myxo16S/data/raw_fastq', pattern = 'fq.gz'),
                    id = map(file, function(x) paste(strsplit(x,"_")[[1]][1],collapse="_")) %>%
                      unlist()) %>%
  group_by(id) %>%
  mutate(., num = n()) %>%
  ungroup()

# check which samples are not present
filter(d, ! sample %in% unique(seq_files$id)) %>%
  write.csv(., 'dada2_pipeline_myxo16S/data/samples_not_present.csv', row.names = FALSE)
# 18, 29, 30, 45, 26
id2 <- d$id



# rename the original files in raw_fastq/original
# grab out IDs
seq_files <- tibble(file = list.files('~/dada2_pipeline_myxo16S/data/raw_fastq/original', pattern = 'fq.gz'),
                    id = map(file, function(x) paste(strsplit(x,"_")[[1]][2],collapse="_")) %>%
                      unlist()) %>%
  group_by(id) %>%
  mutate(., num = n()) %>%
  ungroup()

id2 <- c('18', '26', '29', '30')

# rename files
for(i in 1:length(id2)){
  temp <- filter(seq_files, id == id2[i]) %>% pull(file)
  fwd <- str_subset(temp, 'R1|r1')
  rev <- str_subset(temp, 'R2|r2')
  try(file.rename(file.path('~/dada2_pipeline_myxo16S/data/raw_fastq/original', fwd), file.path('~/dada2_pipeline_myxo16S/data/raw_fastq/original', paste('s', id2[i], '_R1.fq.gz', sep = ''))))
  try(file.rename(file.path('~/dada2_pipeline_myxo16S/data/raw_fastq/original', rev), file.path('~/dada2_pipeline_myxo16S/data/raw_fastq/original', paste('s', id2[i], '_R2.fq.gz', sep = ''))))
}
