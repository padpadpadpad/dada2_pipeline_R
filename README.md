# Outline

A pipeline to run amplicon sequencing analyses for the Buckling lab (and beyond) based on the dada2 workflow implemented in R.

This workflow is designed around the full stack dada2 and phyloseq workflow for microbiome analyses (Callahan _et al._ 2016, [link](https://f1000research.com/articles/5-1492/v2)). More documentation and information on dada2 can be found on the package's [website](https://benjjneb.github.io/dada2/index.html).

The workflow is a work in progress. I am currently trying to update it and aim to have it working with [targets](https://books.ropensci.org/targets/) but that might be a bit ambitious. At the very least we will aggregate all the scripts we already have and make a more sensible, and reuseable, workflow that can be followed.

## Updates

12/07/2023 - Dumped a bunch of new scripts into `scripts/data_processing.R`.

## Issues

Please report any issues to d.padfield@exeter.ac.uk or post in the [Issues tab](https://github.com/padpadpadpad/AB_dada2_pipeline_R/issues). Or if you are at the University of Exeter you can chat to me informally on Microsoft Teams.

## To run

### Data Processing

- Download folder using the big `Clone or download` button at the top right of the page. The folder can then be downloaded as as a ZIP file.
- Place folder in a suitable place but maintain the same folder structure.
- Place raw `.fastq` in the `data/raw_fastq` folder.
- Download [dada2 compatible reference](https://benjjneb.github.io/dada2/training.html) databases and place in `data/ref_trainsets`.
- Open the `.Rproj` file to set the working directory to the root directory.
- Run `scripts/package_install.R` to install all necessary packages needed for the current analyses (or I hope so)!
- Run `scripts/data_processing/check_quality_plots.R` to produce quality plots for each sample to decide on trimming parameters. The file will be saved in `figs`
- Run either `big_data_processing.R` or `raw_read_processing.R` to end up with data suitable for downstream analyses.
- For large datasets, the default way of making a phylogenetic tree in `phangorn` can be painfully slow. If this is the case and the script has created an alignment fasta file, you can run `scripts/data_processing/FastTree_phylo_tree_build.R` and it will make a phylogenetic tree using [FastTree](http://www.microbesonline.org/fasttree/). This is an R script but calls FastTree using __system()__ and will create a new phyloseq object with the tree from FastTree.
- At the end of data processing, it is likely there might be lots of folders in `output` that will not have very many files in, `data_processing/clean_up_folders.R` goes through the folders that contain output and deletes them if they are under a certain size.

### Data Analysis

As I had not done many analyses using amplicon sequence data and am learning all the time, very rudimentary and basic scripts for some of the usual analyses are available here. They are available as __.Rmd__, __html__ and __.R__ files.

- my_first_prevalence_filter:
- my_first_rarefaction_curve:
- my_first_clustering: a script to load in your phyloseq object, perform a clustering analysis and plot the data. It also runs through an example permutational anova using __vegan::adonis()__ and homeogeneity of variances using __vegan::betadisper()__.
- my_first_differential_abundance: a script to perform a very simple differential abundance analysis. The data is agglomerated at a specified taxonomic resolution/ The proportion of each taxa in each sample is then calculated and plotted against treatments. These can then be analysed in a linear model.

#### References

- Callahan, B.J., Sankaran, K., Fukuyama, J.A., McMurdle, P.J. (2016) Bioconductor workflow for microbiome data analysis: from raw reads to community analyses. F1000 Research