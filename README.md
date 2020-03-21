# Cell-type-specific Methylation Analysis

This package provides functions that replicate the analyses reported in Rahmani et al. (2020) comparing TCA (Rahmani et al. 2019) and CellDMC (Zheng et al. 2018).

## Installation

The package has a requirement of R version 3.5.0 or greater (for loading the optional `CellTypeSpecificMethylationData` package). The results found in Rahmani et al. (2020)
were run using R version 3.6.0. Other package requirements are listed in the Imports section of the
DESCRIPTION file.

With the `devtools` package, you can install this package using the following command:

```r
devtools::install_github("cozygene/CellTypeSpecificMethylationAnalysis")
```

## Usage

The functions exported by this package will download most of the source data and output each figure found in Rahmani et al. (2020) into user-specified directories. Note that the data directory will contain just under 5GB of data and may temporarily max out at 13GB for temporary files.

### Parametric Simulations

Figure 1 and Supplementary Figure 7 evaluate the methods when a phenotype affects methylation (X|Y). They can be replicated with the following

```r
parametric_simulations(data_dir, results_dir, plots_dir,
                       experiment_index=1)
```

which will save the simulation data, results, and plots (in tiff format) in the provided directories

Figure 2 and Supplementary Figures 1 and 9 evaluate the methods when a phenotype is affected by methlation, which can be replicated with 

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

Figure 3 visualizes these results and can be replicated as follows:
```r
smoking.res <- smoking_analysis(data_dir, results_dir, plot_dir,
                                hannum_smk_status_path="hannum_smoking_status.txt")
```
where `hannum_smk_status_path` is the path to a CSV file that contains the smoking status of samples in the Hannum et al. dataset. This data is currently not publicly available but is required for this analysis. It contains a header indicating `ID` and `smk` columns that contain the sample IDs and smoking status, respectively. Smoking status is listed as either `Never`, `Past Smoker`, or `Current Smoker` in this file. 


### Cell Fraction Estimation

Supplementary Figure 2 evaluates estimation of cell fractions in two whole-blood datasets (Koestler et al, Jing et al) with FACS-based fractions as ground truth. It can be generated with

```r
refit.w.res <- refit_w_comparison(koestler, bbc, plot_dir)
```

where `koestler` is available in the `CellTypeSpecificMethylationData` analysis and is also publicly available on GEO (GSE77797). Preprocessing of these data in the above package followed procedures described in Lehne et al (Genome Biology 2015). The `bbc` object is also available in the data package. The raw data is available through NODE (OER035661) and was processed with the `epic_qc(idat_dir, data_dir)` function available in this package. Note that the phenotype data must be attached to the `bbc` in slots `W.facs` (Facs counts in matrix with cell types in columns, samples in rows) and `C2` (matrix with 1 column indicating plate number for each sample in `Plate` column name). These phenotype data are under restricted access and were acquired from Jing et al. 

### Runtime Comparison

The runtime of each method was also compared in Supplementary Figure 5. These results can be replicated with the following function:

```r
runtime.res <- runtime_comparison(plot_dir)
```

Note that this function runs each simulation serially (took on the order of a day on a single core).

### TCA Optimization Comparison

Rahmani et al. (2020) also introduces an alternative procedure for fitting the TCA model in an efficient manner. While this generalized method of moments (GMM) approach provides different estimates than the original TCA maximum likelihood (MLE) model, we show in Supplementary Figure 6 that their estimates are highly concordant as sample sizes increase. It can be run with

```r
mle.gmm.comparison.res <- var_mle_gmm_comparison(plot_dir)
```
