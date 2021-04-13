initial\_analysis\_notebook\_030521
================
Nick Lesniak
4/1/2021

### Process sequencing data

  - removed unused fastq files

### run fastqs through mothur

    ## [1] "NP2_2630D7 (23 seqs), OP_931D1 (24 seqs), OP_931D3 (24 seqs), NP2_983D3 (25 seqs), NP2_983D1 (36 seqs), NP2_983D2 (49 seqs), CON3_NT_D5 (62 seqs), NP2_983D4 (76 seqs), NP2_2630D1 (90 seqs), NP1_2725D7 (103 seqs), OP_931D0 (116 seqs), NP2_2630D3 (126 seqs), 578_984D0 (130 seqs), 430_304D10 (403 seqs), NP2_979D7 (473 seqs), IN2_2673_D9 (598 seqs), NP2_2630D0 (921 seqs), OP_931D2 (990 seqs), OUTA_2063_D6 (1184 seqs), NP1_2597D6 (1215 seqs), 369_992D10 (1247 seqs), CON3_NT_D4 (1390 seqs), OUT2_2510_D2 (1476 seqs), CON2_2328_D2 (1541 seqs), 578_984D2 (1720 seqs)"

subsampled to 2017, eliminating the 25 samples listed above with their
counts error rate for the data was 0.19%

# Figure 1 - Human communities reporducibly colonize mice

Sampled spread of diverse donors for inoculating germ-free mice

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

> recalc nmds with just these samples?

Unique communities without conserved structure
![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

> does outcome matter at this point? in this figure?

Donor communities colonize mice and shift community, ~~but
consistently~~

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

# Figure 2 - Mice were colonized without any perturbation

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

![](initial_analysis_notebook_030521_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

> NAs on Day 0 for DA01134, DA10034, DA00884  
> Last samples for moribund mice is from cecum  
> Switch to cage color?  
> Recalc nmds?

# Figure 3 - Difference in severity

![](initial_analysis_notebook_030521_files/figure-gfm/histology-1.png)<!-- -->

Mice were challenged with a strain isolated that matched Clostridioides
difficile ribotype 027

> what is published about histology of 027 in C57/B6 mice?

![](initial_analysis_notebook_030521_files/figure-gfm/toxin-1.png)<!-- -->
