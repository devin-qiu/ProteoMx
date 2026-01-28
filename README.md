# ProteoMx
An R package for analyzing proteomics data of Nanostring GeoMx spatial profiling. 

GeoMx protein sequencing is an antibody-based spatial profiling method. It relies on the IF-based ROI selection first, and then based on the classification within each ROI, 570 protein targets are sequenced for read counts. 

Different antibodies have different binding affinities, so the so-called "negative controls" cannot be treated as numeric zeros. Rather, they are qualitative indicators. 
For spatial analysis on human samples, usually no "standard" dataset paired

Therefore, we will need


# ProteoMx

[![R-CMD-check](https://github.com/USERNAME/ProteoMx/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/USERNAME/ProteoMx/actions/workflows/R-CMD-check.yaml)
**ProteoMx** provides a robust statistical framework for identifying expressed proteins in GeoMx DSP data. 

Instead of using arbitrary signal-to-noise thresholds, ProteoMx uses **Gaussian Mixture Models (GMM)** to mathematically distinguish true biological signal from background noise. This approach is data-driven, adaptive, and ideal for complex spatial proteomics data.

## Installation

You can install the development version of ProteoMx from GitHub:

```r
# install.packages("devtools")
devtools::install_github("devin-qiu/ProteoMx")
```

# Quick Start Workflow
This example demonstrates the standard ProteoMx pipeline: **Normalize** $\rightarrow$ **Fit** $\rightarrow$ **Visualize** $\rightarrow$ **Filter**.
```r
library(ProteoMx)
library(Biobase)
```

# 1. Load Data
We use the tiny example dataset included in the package
```r
data("geomx_set")
```

# 2. Normalize (Q3)
Corrects for library size differences between AOIs
```r
geomx_set <- Q3Normalize(geomx_set)
```

# 3. Fit Mixture Models
Fits Gaussian models (testing 1 to 3 components) for every protein
```r
geomx_set <- MixModelFit(geomx_set, ncomps = 3)
```

# 4. Visualization
Check model fit against the background threshold (Red Line) \n
Example: Plot a known negative control vs. a target protein
```r
PlotMixModel(geomx_set, protein = "Rt IgG2a")      # Negative Control
PlotMixModel(geomx_set, protein = "Cytokeratin 17") # Target Protein
```


# 5. Optimization
```r
best_fit <- BestMixModel(geomx_set, ncomps = 3)
# View results (Best_NComp = optimal number of components)
head(best_fit)
```

# 6. Filter Noise
Remove proteins that never exceed the background threshold (Mean + 1SD of Neg Control)
```r
geomx_filtered <- FilterProteins(geomx_set, neg_ctrl = "Rt IgG2a")
```




