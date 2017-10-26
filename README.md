## Outline

A shell of the pipeline to run sequencing analyses for the Buckling lab (and beyond) based on the dada2 workflow implemented in R.

This workflow is designed around the full stack dada2 and phyloseq workflow for microbiome analyses (Callahan _et al._ 2016, [link](https://f1000research.com/articles/5-1492/v2)). More documentation and information on dada2 can be found on the package's [website](https://benjjneb.github.io/dada2/index.html).

The workflow produces and saves output, figures and a progress file separately and stand alone for each run, making it easy to re-run the pipeline with different parameters and keep previous runs. Jump to the example workthrough to what output is saved and an example progress file.

### Issues

Please report any issues to d.padfield@exeter.ac.uk or post in the [Issues tab](https://github.com/padpadpadpad/AB_dada2_pipeline_R/issues)

#### To run

- Download folder using the big `Clone or download` button at the top right of the page. The folder can then be downloaded as as a ZIP file.
- Place folder in a suitable place but maintain the same folder structure.
- Place raw `.fastq` in the `data/raw_fastq` folder.
- Download [dada2 compatible reference](https://benjjneb.github.io/dada2/training.html) databases and place in `data/ref_trainsets`.
- Open the `.Rproj` file to set the working directory to the root directory.
- Run `scripts/package_install.R` to install all necessary packages needed for the current analyses (or I hope so)!
- Run `scripts/data_processing/check_qual_plot.R` to produce quality plots for each sample to decide on trimming parameters. The file will be saved in `figs`
- Run either `big_data_processing.R` or `raw_read_processing.R` to end up with data suitable for downstream analyses.

#### Workthrough with example data 

#### References

- Callahan, B.J., Sankaran, K., Fukuyama, J.A., McMurdle, P.J. (2016) Bioconductor workflow for microbiome data analysis: from raw reads to community analyses. F1000 Research