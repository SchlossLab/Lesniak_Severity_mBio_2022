ml\_output\_notebook
================
Nick Lesniak
5/21/2021

## Hyperparameter Performance

#### Same day toxin

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

#### Day 0 predict future toxin

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

#### Day 0 predict moribund outcome

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

#### Day 0 predict high/low clinical score on day 10

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

## Model Performance

Features  
\- same day toxin presence - CFU, taxonomic abundances  
\- day 0 future toxin presence - taxonomic abundances

  - day 0 future moribund - median cfu, day w/toxin presence (prior to
    day 3), taxonomic abundances  
  - to avoid duplicating similar data, i sumarized cfu/toxin to be a
    single feature to represent the amount or presence over the first
    two days of challenge since all mice were present for at least 2
    days  
  - day 0 histology (high vs low) - median cfu, days w/toxin, taxonomic
    abundances

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Models and features

#### Modeling future production of toxin based on initial community

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Using Logistic regression at the Genus level or RF at the Phylum level
the features with median differences greater than 0 are:

#### Modeling production of toxin with the community and CFU data from that same day

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Using Random Forest at the Phylum level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

#### Modeling severe disease (moribund) from the initial community

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Using random forest at the Class level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

#### Modeling endpoint histology scores

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Using Logistic regression at the Class level the features with median
differences greater than 0 are:

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

## What are the temporal dynamics of the taxa identified in the Lefse/RF models?

![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-27-4.png)<!-- -->![](ml_analysis_notebook_files/figure-gfm/unnamed-chunk-27-5.png)<!-- -->

#### identify temporal patterns differentiating high/low clinical score

potential last figure or supplemental to show significant temporal
trends
