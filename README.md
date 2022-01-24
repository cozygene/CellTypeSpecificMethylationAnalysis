# Cell-type-specific Methylation Analysis

This package provides functions that reproduce the analyses reported in Rahmani et al. 2021 (<a href="https://www.biorxiv.org/content/10.1101/2021.02.14.431168v1" target="_blank">bioRxiv</a>) comparing TCA (<a href="https://www.nature.com/articles/s41467-019-11052-9" target="_blank">Rahmani et al. 2019, Nature Communications</a>) and CellDMC (<a href="https://www.nature.com/articles/s41592-018-0213-x" target="_blank">Zheng et al. 2018, Nature Methods</a>). These functions were also used to generate the figuers in <a href="https://www.frontiersin.org/articles/10.3389/fbinf.2021.792605/" target="_blank">Rahmani et al. 2022, Frontiers in Bioinformatics</a>; see under "2.1 Software and computational tools" in the <a href="https://www.frontiersin.org/articles/file/downloadfile/792605_supplementary-materials_datasheets_1_pdf/octet-stream/Data%20Sheet%201.pdf/2/792605" target="_blank">Supplementary Information</a> of the latter publication for more details.



## Installation

The package has a requirement of R version 3.5.0 or greater. The results found in Rahmani et al. (2021)
were run using R version 3.6.0. Other package requirements are listed in the Imports section of the
DESCRIPTION file.

With the `devtools` package, you can install this package using the following command:

```r
devtools::install_github("cozygene/CellTypeSpecificMethylationAnalysis")
```

## Usage

The functions exported by this package will download most of the source data and output each figure found in Rahmani et al. (2021) into user-specified directories. Note that the data directory will contain just under 5GB of data and may temporarily max out at 13GB for temporary files.

### Parametric Simulations

Figure 1 and Supplementary Figure 7 evaluate the methods when a phenotype affects methylation (X|Y). They can be reproduced with the following

```r
parametric_simulations(data_dir, results_dir, plots_dir,
                       experiment_index=1)
```

which will save the simulation data, results, and plots (in tiff format) in the provided directories.

Figure 2 and Supplementary Figures 1 and 9 evaluate the methods when a phenotype is affected by methlation (Y|X), which can be reproduced with 

```r
parametric_simulations(data_dir, results_dir, plots_dir,
                       experiment_index=3)
```

Supplmentary Figures 3 and 4 perform the experiments in Figure 1 with smaller sample size and are generated with:


```r
parametric_simulations(data_dir, results_dir, plots_dir,
                       experiment_index=4)
```

Supplementary Figure 8 performs the experiment in Figure 1 but considers the estimated effect direction when determining performance. It can be generated as follows:


```r
parametric_simulations(data_dir, results_dir, plots_dir,
                       experiment_index=2)
```

### Smoking EWAS

The methods were benchmarked on two independent whole-blood datasets (Hannum et al., Liu et al.) .

Figure 3 visualizes these results and can be reproduced as follows:
```r
smoking.res <- smoking_analysis(hannum_smk_status_path="hannum_smoking_status.txt",
                                data_dir, results_dir, plot_dir)
```
where `hannum_smk_status_path` is the path to a CSV file that contains the smoking status of samples in the Hannum et al. dataset. This data is currently not publicly available but is required for this analysis. It contains a header indicating `ID` and `smk` columns that contain the sample IDs and smoking status, respectively. Smoking status is listed as either `Never`, `Past Smoker`, or `Current Smoker` in this file. 


### Cell Fraction Estimation

Supplementary Figure 2 evaluates estimation of cell fractions in two whole-blood datasets (Koestler et al., Jing et al.) with FACS-based fractions as ground truth. It can be reproduced with

```r
refit.w.res <- refit_w_comparison(bbc_pheno_path="PhenoTypesBBC.Rd",
                                  bbc_facs_path="facsBBC.Rd",
                                  data_dir, plot_dir)
```

The Koestler et al. data that is analyzed within this function is available as an object exported by this package and is also publicly available on GEO (GSE77797). Preprocessing of these data in the above package followed procedures described in Lehne et al. (Genome Biology 2015). The BBC data, which is available through NODE (OER035661), is downloaded and processed by this function. Note that the phenotype and FACS data are currently under restriced acces. The `bbc_pheno_path` and `bbc_facs_path` parameters point to the Rd files that contain these data and were provided by Jing et al.

### Runtime Comparison

The runtime of each method was also compared in Supplementary Figure 5. These results can be reproduced with the following function:

```r
runtime.res <- runtime_comparison(plot_dir)
```

Note that this function runs each simulation serially (took on the order of a day on a single core).

### TCA Optimization Comparison

Rahmani et al. (2021) also introduces an alternative procedure for fitting the TCA model in an efficient manner. While this generalized method of moments (GMM) approach provides different estimates than the original TCA maximum likelihood (MLE) model, we show in Supplementary Figure 6 that their estimates are highly concordant as sample sizes increase. It can be run with

```r
mle.gmm.comparison.res <- var_mle_gmm_comparison(plot_dir)
```
