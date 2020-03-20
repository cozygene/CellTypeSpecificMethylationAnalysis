# Cell-type-specific Methylation Analysis

This package provides functions that replicate the analyses reported in Rahmani et al. (2020).

## Installation

The package has a requirement of R version 3.4.0 or greater. The results found in Rahmani et al. (2020)
were run using R version 3.6.0. Other package requirements are listed in the Imports section of the
DESCRIPTION file.

With the `devtools` package, you can install this package using the following command:

```r
devtools::install_github("cozygene/CellTypeSpecificMethylationAnalysis")
```

## Usage

The main entrypoints for each analysis are as follows:

 * `compare_runtimes(results_dir, plot_type)`

This function will compare the runtime of TCA and CellDMC on various sample sizes and number of CpGs. Note that this function will take many hours to run.

* `parametric_simulations(data_path, results_path, plots_path, experiment_index=c(1,2,3), num_cores=1, image_format="tiff)`

This function will run one of three parametric simulations under different model assumptions. Specifically, 1 indicates phenotype affecting methylation levels, 2 indicates phenotype affecting methylation levels where effect directions are considered in power analyses, and 3 indicates methylation levels affecting the phenotype. 

* `refit_w_comparison(koestler, bbc, results_dir, plot_type)`

This function compares the cell type fraction estimation of different modes of TCA on two real datasets (provided in CellTypeSpecificMethylationData package).

* `smoking_analysis(liu, hannum, results_dir, plot_type)`

This function performs cell-type-specific methylation analyses of two real datasets (provided in CellTypeSpecificMethylationData package) in the context of smoking.

* `var_mle_gmm_comparison(results_dir, plot_type)`

This function compares TCA's MLE and GMM modes for optimization of parameters.
