visualization416S package demo
================
yeguanhua

Install and load package
========================

Note:

You need to run install from the parent working directory that contains the visualization416S folder.

If you are can't install package in R Console, try this: open visualization416S.Rproj first, then use RStudio → Build → Install and Restart.

``` r
devtools::install('visualization416S')
```

Load package and demo data.

``` r
rm(list = ls())
library(visualization416S)
demoPhyloseq <- visualization416S::demo_phyloseq_object
demoDADA2res <- visualization416S::demo_dada2_result
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

| subject\_id | diagnosis         |
|:------------|:------------------|
| s17118657   | healthy           |
| s17118661   | healthy           |
| s17118667   | healthy           |
| s17118714   | healthy           |
| s17118646   | intestinal cancer |
| s17118664   | intestinal cancer |
| s17118669   | intestinal cancer |
| s17118686   | intestinal cancer |
| s17118647   | liver cancer      |
| s17118684   | liver cancer      |
| s17118715   | liver cancer      |
| s17118730   | liver cancer      |
| s17118650   | lung cancer       |
| s17118680   | lung cancer       |
| s17118691   | lung cancer       |
| s17118703   | lung cancer       |

Track reads through DADA2 workflow
==================================

First use the dada2\_reads\_track function to check reads drop associated with every step in DADA2.

``` r
visualization416S::dada2_reads_track(demoDADA2res$reads_track, single_end = FALSE)
```

![](README_files/figure-markdown_github/reads%20track-1.png)

Stacked bar plot of phylogenetic composition
============================================

Use stacked\_bar\_plot function to plot the Order level abundance in every sample. You can change the 'level' argument to plot abundance in different level, or change the 'feature' to choose different feature you want to show in x-axis.

### Family level

``` r
visualization416S::stacked_bar_plot(phyloseq = demoPhyloseq, level = "Family", feature = "diagnosis")
```

![](README_files/figure-markdown_github/Stacked%20bar%20plot-1.png)

Alpha diversity
===============

Use alpha\_diversity\_plot to plot alpha diversity. Change 'measures' argument to use different measurement.

### Chao1

``` r
visualization416S::alpha_diversity_plot(phyloseq = demoPhyloseq, feature = "diagnosis", measures = "Chao1", 
                           p_test = "kruskal")
```

![](README_files/figure-markdown_github/Chao1-1.png)

Beta diversity
==============

Use beta\_diversity\_plot to plot beta diversity. Change method to draw different beta diversity plot.

### Bray-Curtis

``` r
visualization416S::beta_diversity_plot(phyloseq = demoPhyloseq, feature = "diagnosis", method = "bray")
```

![](README_files/figure-markdown_github/Bray-Curtis-1.png)

Log2 fold change
================

Use log2fc function to show differential analysis result.

``` r
visualization416S::log2fc(phyloseq = demoPhyloseq, feature = "diagnosis", level = NA, 
             p_value = 0.05, save_res = TRUE)
```

    ## [1] "DESeq2_result.rds has been saved to current working directory."
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

``` r
res <- readRDS('DESeq2_result.rds')
res
```

    ## log2 fold change (MLE): diagnosis lung.cancer vs healthy 
    ## Wald test p-value: diagnosis lung.cancer vs healthy 
    ## DataFrame with 1420 rows and 6 columns
    ##                baseMean     log2FoldChange            lfcSE
    ##               <numeric>          <numeric>        <numeric>
    ## OTU154 67.0226773751778   22.4015248306817 2.94452289002941
    ## OTU152 120.431649147706   -25.344577557431 3.80136683588675
    ## OTU109 121.962118499815  -22.1997536375473 3.80003022265851
    ## OTU164 52.9833807573347   21.8558119093275 4.20650098383088
    ## OTU312 19.8626227371063   20.7402498894252 4.20721933688743
    ## ...                 ...                ...              ...
    ## OTU652 3.01559579269091 -0.355409258234053 4.32794324774993
    ## OTU68  189.678929913974   1.51568862368933 1.61528457847523
    ## OTU347 14.2779775834037  -22.7173513024046 4.20748608897282
    ## OTU391 11.0604051702423  -22.4012225492369 4.20791261335844
    ## OTU283 11.7599955447353  -11.4582019101067 4.32794324771952
    ##                       stat               pvalue                 padj
    ##                  <numeric>            <numeric>            <numeric>
    ## OTU154    7.60786234895188 2.78666394581598e-14 2.80617059343669e-11
    ## OTU152   -6.66722751357903 2.60680829886671e-11 1.31252797847939e-08
    ## OTU109   -5.84199396762067 5.15796747152172e-09 1.73135774794079e-06
    ## OTU164    5.19572252409728 2.03926239000198e-07 5.13384306682998e-05
    ## OTU312    4.92968115723891 8.23639291399349e-07 0.000165880953287829
    ## ...                    ...                  ...                  ...
    ## OTU652 -0.0821196669847342                   NA                   NA
    ## OTU68    0.938341542962096                   NA                   NA
    ## OTU347   -5.39926949775147                   NA                   NA
    ## OTU391   -5.32359500007722                   NA                   NA
    ## OTU283    -2.6474935677921                   NA                   NA
