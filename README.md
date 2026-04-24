# Background
**ProteoMx**: An R package for analyzing proteomics data of Nanostring GeoMx spatial profiling. 

GeoMx protein sequencing is an antibody-based spatial profiling method. It relies on the IF-based Region of Interest (ROI) selection first, and then based on the classification within each ROI, up to 570 protein targets are sequenced for read counts. 

Different antibodies have different binding affinities and specificity to different tissue types, so the protein read values that are below the so-called "negative controls" in the protein panel cannot be treated as numeric zeros. Rather, "negative controls" are qualitative indicators. 

GeoMx does not reach single-cell resolution. Protein expression are measured with the minimal units as Areas of Interest (AOIs), which are subsets within ROIs. that being said, 

Moreover, for spatial analysis on human samples, usually no "standard" dataset paired compared to other proteomics approaches such as mass spec. 

Therefore, we will need a tool to address this issue, by applying our statistical model to 1) identify proteins that has true positive expression values in comparison to negative probe values; 2) predict cell subpopulation based on the components identified in our statistical model; 3) provide robust insights from the spatial data standing alone. 

**ProteoMx** provides a robust statistical framework for identifying expressed proteins in GeoMx DSP data. 

Instead of using arbitrary signal-to-noise thresholds, ProteoMx uses **Gaussian Mixture Models (GMM)** to mathematically distinguish true biological signal from background noise. This approach is data-driven, adaptive, and ideal for complex spatial proteomics data.

# Quick Start Workflow
This example demonstrates the standard ProteoMx pipeline:

**Normalize -> Fit Models -> Select Best Model -> Filter Noise -> Address Heteroscedasticity (Decay) -> QC & Visualize**.

## Installation

You can install the development version of ProteoMx from GitHub:

```r
# install.packages("devtools")
devtools::install_github("devin-qiu/ProteoMx")
library(ProteoMx)
library(Biobase)
```

# 1. Load Data
We use the tiny example dataset included in the package
```r
data("geomx_set")
```

# 2. Normalize (Q3)
A required preprocessing step. Performs 3rd quantile (Q3) normalization to correct for library size differences between AOIs.

`geomx_set`: A `NanoStringGeoMxSet` object.

*Output*: Adds the normalized matrix to the object under `assayDataElement(object, "q_norm")`.
```r
geomx_set <- Q3Normalize(geomx_set)
# This function is equivalent to 
# NanoStringNCTools::normalize(..., norm_method = "quant", desiredQuantile = 0.75, toElt = "q_norm")
```

# 3. Fit Mixture Models
Iterates through every protein and fits Gaussian mixture models for every combination of components (1 to ncomps) and variance structures (Equal and Unequal).

`geomx_set`: A Q3-normalized `NanoStringGeoMxSet` object.

`ncomp`: the number of components to be fit in the mixture models. 

*Output*: Safely stores all mathematical fits inside `experimentData(geomx_set)@other$MixModel`.

```r
geomx_set <- MixModelFit(geomx_set, ncomps = 3)
```

# 4. Optimize & Select Best Model
Runs a Bayesian tournament (via BIC and Bayes Factors) to determine the mathematically optimal component structure for every single protein, applying Occam's Razor to prevent overfitting.

`geomx_set`: The object processed by `MixModelFit`.

`ncomps`: Integer. Must match the maximum number of components run in the fitting step.

*Output*: Safely stores a summary table of the winning models in `experimentData(geomx_set)@other$Best_Model_Summary`.

```r
geomx_set <- BestMixModel(geomx_set, ncomps = 3)

# View the summary of winning models
best_fit <- experimentData(geomx_set)@other$Best_Model_Summary
head(best_fit)
```

# 5. Filter Noise
Removes proteins that never meaningfully exceed the background threshold in their statistically optimal model.

`geomx_set`: The object processed by `BestMixModel`.

`neg_ctrl`: Character. The specific negative control probe used to establish the baseline noise floor (e.g., `"Rb IgG"` or `"Rt IgG2a"`).

`n_sd`: Numeric. The number of standard deviations above the mean of the negative control to set the strict threshold. Default is `1`.

*Note*: Negative control probes are mathematically excluded from the biological "Pass" list but are silently preserved in the output object for downstream plotting.

```r
geomx_set_filtered <- FilterProteins(geomx_set, neg_ctrl = "Rb IgG", n_sd = 1)
```

# 6. Address Heteroscedasticity (Soft Thresholding)
Different antibodies possess varying binding affinities and background "stickiness", creating heteroscedasticity across targets. `DecayBackground` addresses this by applying an Empirical Cumulative Distribution Function (ECDF)-weighted exponential decay to sub-threshold reads. Data deep in the noise tail is aggressively squashed toward zero, while data near the threshold is spared, preserving natural variance for downstream differential expression.

`geomx_set`: The filtered `NanoStringGeoMxSet` object.

`neg_ctrl`: Character. Must match the negative control used during filtering.

`n_sd`: Numeric. The standard deviation multiplier for the threshold.

`base_decay`: Numeric. The baseline penalty multiplier. A higher number squashes the background noise more aggressively. Default is `2.0`.

*Output*: Creates a brand new, analysis-ready matrix inside the object called `"exp_decayed"`.

```r
geomx_set_filtered <- DecayBackground(geomx_set_filtered, neg_ctrl = "Rb IgG", n_sd = 1, base_decay = 2.0)
```

# 7. Quality Control (Ambiguity Flagging)
Cross-references the filtered data to flag highly expressed proteins that failed the detection threshold due to complex distribution overlaps (where the top 2 models were statistically indistinguishable).

`geomx_set`: The filtered object.

*Output*: Prints a summary to the console and saves a review table into the object metadata.

```r
geomx_set_filtered <- CheckBestModel(geomx_set_filtered)

# View proteins that may require manual inspection
review_table <- experimentData(geomx_set_filtered)@other$Ambiguous_Models_Review
```

# 8. Visualization
Generate density plots to visually assess the Gaussian fits and the threshold cutoffs. You can visualize both the original normalized data and your newly penalized data.

`geomx_set`: The processed object.

`protein`: Character. The target protein name to plot.

`assay`: Character. Which matrix to visualize (`"q_norm"` or `"exp_decayed"`). Default is `"q_norm"`.

`neg_ctrl`: Character or Vector. The threshold(s) to draw as red dashed lines. Accepts multiple controls.

`ncomp` & `ev`: (Optional) Force the plot to show a specific component number or variance assumption. If left as `NULL`, the function automatically plots the optimal model identified by `BestMixModel`.

```r
# Visualize the original mathematical fit
PlotMixModel(geomx_set_filtered, protein = "CD4", assay = "q_norm", neg_ctrl = c("Rb IgG", "Ms IgG1", "Ms IgG2a"))

# Visualize the data after Soft-Thresholding 
PlotMixModel(geomx_set_filtered, protein = "CD4", assay = "exp_decayed", neg_ctrl = c("Rb IgG", "Ms IgG1", "Ms IgG2a"))
```

# 9. Extracting Data for Downstream Analysis
Once the pipeline is complete, the `"exp_decayed"` assay is ready to be fed directly into your preferred differential expression tool (like `limma` or `standR`). You can extract the final matrix using standard Biobase commands:

```r
final_matrix <- Biobase::assayDataElement(geomx_set_filtered, "exp_decayed")
```
