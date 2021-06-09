ml\_output\_notebook
================
Nick Lesniak
5/21/2021

## Hyperparameter Performance

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

## Model Performance

Features  
\- same day toxin presence - CFU, taxonomic abundances  
\- day 0 future toxin presence - taxonomic abundances  
\- day 0 future moribund - median cfu, day w/toxin presence (prior to
day 3), taxonomic abundances  
to avoid duplicating similar data, i sumarized cfu/toxin to be a single
feature to represent the amount or presence over the first two days of
challenge since all mice were present for at least 2 days - day 10
histology (high/mid/low) - cfu, toxin level, taxonomic abundances

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Models and features

#### Modeling future production of toxin based on initial community

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Using Logistic regression at the Genus level or RF at the Phylum level
the features with median differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

#### Modeling production of toxin with the community and CFU data from that same day

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Using Logistic regression at the Phylum level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

#### Modeling severe disease (moribund) from the initial community

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Using Logistic regression at the Class level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Modeling endpoint histology scores

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Using multiclass RF at the Phylum level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
