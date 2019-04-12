visual16S package demo
================
yeguanhua

Install and load package
========================

Note:

You need to run install from the parent working directory that contains the visual16S folder.

``` r
devtools::install('visual16S')
```

#### If you are can't install package in R Console, try this: open visual16S.Rproj first, then use RStudio → Build → Install and Restart.

Load package and demo data.
===========================

``` r
rm(list = ls())
library(visual16S)
demo_phyloseq_object <- visual16S::demo_phyloseq_object
demo_dada2_result <- visual16S::demo_dada2_result
```

Data status
===========

Primer:

CCTAYGGGRBGCASCAG ; GGACTACNNGGGTATCTAAT

DADA2 filter parameters:

dada2::filterAndTrim(truncLen=c(0,0), maxEE=c(2,2))

DADA2 taxonomy database:

silva\_nr\_v132

Metadata:

| SampleID  | diagnosis         |
|:----------|:------------------|
| s17118657 | healthy           |
| s17118661 | healthy           |
| s17118667 | healthy           |
| s17118714 | healthy           |
| s17118646 | intestinal cancer |
| s17118664 | intestinal cancer |
| s17118669 | intestinal cancer |
| s17118686 | intestinal cancer |
| s17118647 | liver cancer      |
| s17118684 | liver cancer      |
| s17118715 | liver cancer      |
| s17118730 | liver cancer      |
| s17118650 | lung cancer       |
| s17118680 | lung cancer       |
| s17118691 | lung cancer       |
| s17118703 | lung cancer       |

DADA2 workflow reads track
==========================

First use the track\_reads\_dada2 function to check reads drop associated with every step in DADA2.

Note:

You can plot reads track in either absolute abundance or relative abundance.

If the sample size is too large to show on legend, set legend\_position to "none".

``` r
track_reads_dada2(demo_dada2_result$reads_track, single_end = FALSE, relative_abundance = TRUE, 
                  legend_position = "top")
```

![](README_files/figure-markdown_github/reads%20track-1.png)

Stacked bar plot of phylogenetic composition
============================================

Use plot\_stacked\_bar function to plot the Order level abundance in every sample. You can change the 'level' argument to plot abundance in different level, or change the 'feature' to choose different feature you want to show in x-axis.

Note:

If the legend is too much to show, set legend\_position to "none".

### Order level

``` r
plot_stacked_bar(phyloseq = demo_phyloseq_object, level = "Order", feature = "diagnosis", 
                 x_size = 8, legend_position = "top", legend_size = 10)
```

![](README_files/figure-markdown_github/Stacked%20bar%20plot-1.png)

Alpha diversity
===============

Use plot\_alpha\_diversity to plot alpha diversity. Change 'measures' argument to use different measurement.

### Chao1

``` r
plot_alpha_diversity(phyloseq = demo_phyloseq_object, feature = "diagnosis", feature2 = NA, 
                     measures = "Chao1", p_test = "kruskal")
```

![](README_files/figure-markdown_github/Chao1-1.png)

Beta diversity
==============

Use plot\_beta\_diversity to plot beta diversity. Change 'method' to draw different beta diversity plot.

### Bray-Curtis

``` r
plot_beta_diversity(phyloseq = demo_phyloseq_object, feature = "diagnosis", feature2 = NA, 
                    method = "bray")
```

![](README_files/figure-markdown_github/Bray-Curtis-1.png)

Log2 fold change
================

Use log2fc function to show differential analysis result.

``` r
log2fc(phyloseq = demo_phyloseq_object, feature = "diagnosis", level = NA, p_value = 0.05)
```

    ## [1] "log2 fold change (MLE): diagnosis lung.cancer vs healthy"
    ##      OTU log2FoldChange         padj
    ## 1 OTU154       22.40152 2.806171e-11
    ## 4 OTU164       21.85581 5.133843e-05
    ## 5 OTU312       20.74025 1.658810e-04
    ## 8 OTU228       19.48117 4.637397e-04
    ## 6 OTU331      -18.93621 1.699453e-04
    ## 7 OTU282      -20.06961 2.784319e-04
    ## 3 OTU109      -22.19975 1.731358e-06
    ## 2 OTU152      -25.34458 1.312528e-08

![](README_files/figure-markdown_github/log2fc-1.png)
