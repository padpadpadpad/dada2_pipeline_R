# extra functions ####

# pipeline setup
dada2_pipeline_setup <- function (packages = c("ggplot2", "dada2", "phyloseq", "DECIPHER", 
                                               "tidyr", "dplyr"), 
                                  raw_path = "data/raw_fastq", 
                                  filt_path = "data/filtered_fastq", 
                                  plot_path = "plots", 
                                  output_path = "data/output", 
                                  progress_path = "data/progress", 
                                  ref_fasta = "data/ref_trainsets/rdp_train_set_16.fa", 
                                  ref_fasta_spp = "data/ref_trainsets/rdp_species_assignment_16.fa", 
                                  meta_data = "data/metadata.csv",
                                  run = '1',
                                  pool = 'pseudo') 
  {
  if (any(packages %in% utils::installed.packages()) == FALSE) {
    pkgs_not_pres <- packages[!packages %in% utils::installed.packages()]
    stop(paste("Some of the packages are currently not installed. Please install them before running this function again. Packages not present are:", 
               paste(pkgs_not_pres, collapse = " "), sep = " "), 
         call. = FALSE)
  }
  lapply(packages, library, character.only = TRUE)
  assign("raw_path", raw_path, envir = globalenv())
  assign("filt_path", filt_path, envir = globalenv())
  assign("output_path", output_path, envir = globalenv())
  assign("ref_fasta", ref_fasta, envir = globalenv())
  assign("ref_fasta_spp", ref_fasta_spp, envir = globalenv())
  run = paste('run', run, sep = '_')
  
  if(file.exists(paste(progress_path, "/", run, "_progress.txt", sep = "")) == FALSE)
    {file.create(paste(progress_path, "/", run, "_progress.txt", sep = ""))
    assign("progress_file", paste(progress_path, "/", run, "_progress.txt", 
                                  sep = ""), envir = globalenv())
    writeLines(paste("Run started at ", Sys.time()), progress_file)
  }
  
  assign("progress_file", paste(progress_path, "/", run, "_progress.txt", 
                                sep = ""), envir = globalenv())
  dir.create(file.path(plot_path, run))
  assign("plot_path", file.path(plot_path, run), envir = globalenv())
  dir.create(file.path(output_path, run))
  dir.create(paste(output_path, run, 'temp', sep ='/'))
  assign("output_path", file.path(output_path, run), envir = globalenv())
  
  assign("meta", utils::read.csv(meta_data, stringsAsFactors = FALSE), 
         envir = globalenv())
  assign('pool', pool, envir = globalenv())
}

# geom_pretty boxplot
geom_pretty_boxplot <- function (...) 
{
  list(ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.75), 
                             ...), ggplot2::stat_summary(geom = "crossbar", position = ggplot2::position_dodge(width = 0.75), 
                                                         fatten = 0, color = "white", width = 0.4, fun.data = function(x) {
                                                           return(c(y = stats::median(x), ymin = stats::median(x), 
                                                                    ymax = stats::median(x)))
                                                         }))
}

# track reads through pipeline
get_sample_names <- function(files, sample_ids){
  ids <- gsub('.*_', '', sample_ids)
  temp <- files
  for(i in 1:length(ids)){
    temp[stringr::str_detect(temp, ids[i])] <- paste('sample_', ids[i])
  }
  return(temp)
}
