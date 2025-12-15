# ProteoMx
An R package for analyzing proteomics data of Nanostring GeoMx spatial profiling. 

GeoMx protein sequencing is an antibody-based spatial profiling method. It relies on the IF-based ROI selection first, and then based on the classification within each ROI, 570 protein targets are sequenced for read counts. 

Different antibodies have different binding affinities, so the so-called "negative controls" cannot be treated as numeric zeros. Rather, they are qualitative indicators. 
For spatial analysis on human samples, usually no "standard" dataset paired

Therefore, we will need
