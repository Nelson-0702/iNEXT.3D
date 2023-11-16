<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.3D (R package)

<h5 align="right">
Latest version: 2023-11-16
</h5>
<font color="394CAE">
<h3 color="394CAE" style="font-weight: bold">
Introduction to iNEXT.3D (R package): Excerpt from iNEXT.3D User’s Guide
</h3>
</font> <br>
<h5>
<b>Anne Chao, Kai-Hsiang Hu</b> <br><br> <i>Institute of Statistics,
National Tsing Hua University, Hsin-Chu, Taiwan 30043</i>
</h5>

<br> iNEXT.3D (INterpolation and EXTrapolation for three dimensions of
diversity) is an R package, available in
[Github](https://github.com/AnneChao), for rarefaction and extrapolation
of species diversity (Hill numbers) for three dimensions. Here we
provide a quick introduction demonstrating how to run iNEXT.3D, and
showing three types of rarefaction/extrapolation sampling curves. See
Chao et al. (2021) for methodologies. An online version of [iNEXT.3D
Online](https://chao.shinyapps.io/iNEXT_3D/) is also available for users
without an R background. Detailed information about all functions in
iNEXT.3D is provided in the iNEXT.3D Manual in
[iNEXT.3D_vignettes](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/A%20Quick%20Introduction%20to%20iNEXT.3D%20via%20Examples.html),
which is also available from [Anne Chao’s
website](http://chao.stat.nthu.edu.tw/wordpress/software_download/).

`iNEXT.3D` is the extension of R package
[iNEXT](https://cran.r-project.org/web/packages/iNEXT/index.html) (Hsieh
et al., 2016). `iNEXT.3D` focuses on three measures of Hill numbers of
order q: species richness (`q = 0`), Shannon diversity (`q = 1`, the
exponential of Shannon entropy) and Simpson diversity (`q = 2`, the
inverse of Simpson concentration) and extend Hill number to three
dimensions: taxonomic diversity (TD), phylogenetic diversity (PD), and
functional diversity (FD) under Hill-Chao family frame work (Chao et
al., 2019). Besides, `iNEXT.3D` also provide statistic estimation for
three dimensions biodiversity (Chao et al., 2021). For each diversity
measure, `iNEXT.3D` uses the observed sample of abundance or incidence
data (called the “reference sample”) to compute diversity estimates and
the associated confidence intervals for the following two types of
rarefaction and extrapolation (R/E):

1.  Sample-size-based R/E sampling curves:`iNEXT3D` computes diversity
    estimates for rarefied and extrapolated samples up to an appropriate
    size. This type of sampling curve plots the diversity estimates with
    respect to sample size.  
2.  Coverage-based R/E sampling curves: `iNEXT3D` computes diversity
    estimates for rarefied and extrapolated samples with sample
    completeness (as measured by sample coverage) up to an appropriate
    coverage. This type of sampling curve plots the diversity estimates
    with respect to sample coverage.

`iNEXT.3D` also plots the above two types of sampling curves and a
sample completeness curve by `ggiNEXT3D`. The sample completeness curve
provides a bridge between these two types of curves.

### SOFTWARE NEEDED TO RUN INEXT.3D IN R

-   Required: [R](http://cran.rstudio.com/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### HOW TO RUN INEXT.3D:

The iNEXT.3D package is available on
[Github](https://github.com/AnneChao/iNEXT.3D) and can be downloaded
with a standard installation procedure using the commands shown below.
For a first-time installation, an additional visualization extension
package (ggplot2) must be installed and loaded.

``` r
## install iNEXT.3D package from CRAN
# install.packages("iNEXT.3D")  # coming soon

## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('AnneChao/iNEXT.3D')

## import packages
library(iNEXT.3D)
library(ggplot2)
```

In this document, here provide a quick introduction demonstrating how to
run the package `iNEXT.3D` (iNterpolation and EXTrapolation in three
Dimensions). `iNEXT.3D` has several main functions: `iNEXT3D`,
`ggiNEXT3D`, `ObsAsy3D`, `ggObsAsy3D`, `estimate3D`, and `DataInfo3D.`

### MAIN FUNCTION: iNEXT3D()

The main function iNEXT3D() with default arguments is described below:

``` r
iNEXT3D(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", 
        size = NULL, endpoint = NULL, knots = 40, nboot = 50, conf = 0.95, nT = NULL, 
        PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD', 
        FDdistM, FDtype = 'AUC', FDtau = NULL)
```

This main function computes diversity estimates of order q, the sample
coverage estimates and related statistics for K (if `knots = K`)
evenly-spaced knots (sample sizes) between size 1 and the `endpoint`,
where the endpoint is described below. Each knot represents a particular
sample size for which diversity estimates will be calculated. By
default, `endpoint` = double the reference sample size for abundance
data or double the total sampling units for incidence data. For example,
if `endpoint = 10`, `knot = 4`, diversity estimates will be computed for
a sequence of samples with sizes (1, 4, 7, 10).

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style="text-align: left;">
(a). For datatype = ‘abundance’, data can be input as a vector of
species abundances (for a single assemblage), matrix/data.frame (species
by assemblages), or a list of species abundance vectors. (b). For
datatype = ‘incidence_freq’, data can be input as a vector of incidence
frequencies (for a single assemblage), matrix/data.frame (species by
assemblages), or a list of incidence frequencies; the first entry in all
types of input must be the number of sampling units in each assemblage.
(c). For datatype = ‘incidence_raw’, data can be input as a list of
matrix/data.frame (species by sampling units); data can also be input as
a matrix/data.frame by merging all sampling units across assemblages
based on species identity; in this case, the number of sampling units
(nT, see below) must be input.
</td>
</tr>
<tr>
<td style="text-align: left;">
diversity
</td>
<td style="text-align: left;">
selection of diversity type: ‘TD’ = Taxonomic diversity, ‘PD’ =
Phylogenetic diversity, and ‘FD’ = Functional diversity.
</td>
</tr>
<tr>
<td style="text-align: left;">
q
</td>
<td style="text-align: left;">
a numerical vector specifying the diversity orders. Default is c(0, 1,
2).
</td>
</tr>
<tr>
<td style="text-align: left;">
datatype
</td>
<td style="text-align: left;">
data type of input data: individual-based abundance data (datatype =
‘abundance’), sampling-unit-based incidence frequencies data (datatype =
‘incidence_freq’), or species by sampling-units incidence matrix
(datatype = ‘incidence_raw’) with all entries being 0 (non-detection) or
1 (detection).
</td>
</tr>
<tr>
<td style="text-align: left;">
size
</td>
<td style="text-align: left;">
an integer vector of sample sizes for which diversity estimates will be
computed. If NULL, then diversity estimates will be calculated for those
sample sizes determined by the specified/default endpoint and knots;
</td>
</tr>
<tr>
<td style="text-align: left;">
endpoint
</td>
<td style="text-align: left;">
an integer specifying the sample size that is the endpoint for R/E
calculation; If NULL, then endpoint=double the reference sample size;
</td>
</tr>
<tr>
<td style="text-align: left;">
knots
</td>
<td style="text-align: left;">
an integer specifying the number of equally-spaced knots (say K, default
is 40) between size 1 and the endpoint;each knot represents a particular
sample size for which diversity estimate will be calculated. If the
endpoint is smaller than the reference sample size, then iNEXT3D()
computes only the rarefaction esimates for approximately K evenly spaced
knots. If the endpoint is larger than the reference sample size, then
iNEXT3D() computes rarefaction estimates for approximately K/2 evenly
spaced knots between sample size 1 and the reference sample size, and
computes extrapolation estimates for approximately K/2 evenly spaced
knots between the reference sample size and the endpoint.
</td>
</tr>
<tr>
<td style="text-align: left;">
nboot
</td>
<td style="text-align: left;">
a positive integer specifying the number of bootstrap replications when
assessing sampling uncertainty and constructing confidence intervals.
Enter 0 to skip the bootstrap procedures. Default is 50.
</td>
</tr>
<tr>
<td style="text-align: left;">
conf
</td>
<td style="text-align: left;">
a positive number \< 1 specifying the level of confidence interval.
Default is 0.95.
</td>
</tr>
<tr>
<td style="text-align: left;">
nT
</td>
<td style="text-align: left;">
(required only when datatype = ‘incidence_raw’ and input data is
matrix/data.frame) a vector of nonnegative integers specifying the
number of sampling units in each assemblage. If assemblage names are not
specified, then assemblages are automatically named as ‘assemblage1’,
‘assemblage2’,…, etc.
</td>
</tr>
<tr>
<td style="text-align: left;">
PDtree
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’), a phylogenetic tree in Newick
format for all observed species in the pooled assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
PDreftime
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’), a vector of numerical values
specifying reference times for PD. Default is NULL (i.e., the age of the
root of PDtree).
</td>
</tr>
<tr>
<td style="text-align: left;">
PDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’), select PD type: PDtype = ‘PD’
(effective total branch length) or PDtype = ‘meanPD’ (effective number
of equally divergent lineages). Default is ‘meanPD’, where meanPD =
PD/tree depth.
</td>
</tr>
<tr>
<td style="text-align: left;">
FDdistM
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’), select FD type: FDtype =
‘tau_values’ for FD under specified threshold values, or FDtype = ‘AUC’
(area under the curve of tau-profile) for an overall FD which integrates
all threshold values between zero and one. Default is ‘AUC’.
</td>
</tr>
<tr>
<td style="text-align: left;">
FDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’), select FD type: FDtype =
‘tau_values’ for FD under specified threshold values, or FDtype = ‘AUC’
(area under the curve of tau-profile) for an overall FD which integrates
all threshold values between zero and one. Default is ‘AUC’.
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
FDtau
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
(required only when diversity = ‘FD’ and FDtype = ‘tau_values’), a
numerical vector between 0 and 1 specifying tau values (threshold
levels). If NULL (default), then threshold is set to be the mean
distance between any two individuals randomly selected from the pooled
assemblage (i.e., quadratic entropy).
</td>
</tr>
</tbody>
</table>

This function returns an “iNEXT3D” object which can be further used to
make plots by the function `ggiNEXT3D()` to described later.

### DATA FORMAT/INFORMATION

Three types of data are supported:

1.  Individual-based abundance data (`datatype = "abundance"`): Input
    data for each assemblage/site include samples species abundances in
    an empirical sample of n individuals (“reference sample”). When
    there are N assemblages, input data consist of an S by N abundance
    matrix, or N lists of species abundances.

2.  Sampling-unit-based incidence data: There are two kinds of input
    data.  

<!-- -->

1.  Incidence-raw data (`datatype = "incidence_raw"`): for each
    assemblage, input data for a reference sample consist of a
    species-by-sampling-unit matrix; when there are N assemblages, input
    data consist of N lists of matrices, and each matrix is a
    species-by-sampling-unit matrix. The matrix of combined assemblage
    is allowed, but nT must be specified (see above description).

2.  Incidence-frequency data (`datatype = "incidence_freq"`): input data
    for each assemblage consist of species sample incidence frequencies
    (row sums of each incidence matrix). When there are N assemblages,
    input data consist of an (S + 1) by N matrix, or N lists of species
    incidence frequencies. The first entry of each list must be the
    total number of sampling units, followed by the species incidence
    frequencies.

``` r
data("Brazil_rainforest_data")

out.TD <- iNEXT3D(data = Brazil_rainforest_data$data, diversity = "TD", 
                  q = c(0, 1, 2), datatype = "abundance", 
                  nboot = 10)
out.PD <- iNEXT3D(data = Brazil_rainforest_data$data, diversity = "PD", 
                  q = c(0, 1, 2), datatype = "abundance", 
                  PDtree = Brazil_rainforest_data$tree, 
                  nboot = 10)
out.FD <- iNEXT3D(data = Brazil_rainforest_data$data, diversity = "FD", 
                  q = c(0, 1, 2), datatype = "abundance", 
                  FDdistM = Brazil_rainforest_data$dist,
                  nboot = 5)
```

``` r
out.TD$TDInfo
  Assemblage    n S.obs     SC(n)    SC(2n)  f1 f2 f3 f4 f5
1       Edge 1794   319 0.9387143 0.9744012 110 48 38 28 13
2   Interior 2074   356 0.9407166 0.9728436 123 48 41 32 19
```

Second part of output from function `iNEXT3D` is diversity estimates and
related statistics computed for these 40 knots by default (for example
in “EM” assemblage, corresponding to sample sizes m = 1, 20, 40, …, 372,
373, 374, …, 746), which locates the reference sample size at the
mid-point of the selected knots. The diversity can be based on
sample-size-based and sample coverage-based. The first data frame of
list `$TDiNextEst` (as shown below for ‘size_based’) includes the sample
size (`m`), the `Method` (`Rarefaction`, `Observed`, or `Extrapolation`,
depending on whether the size `m` is less than, equal to, or greater
than the reference sample size), the diversity order (`Order.q`), the
diversity estimate of order q, the lower and upper confidence limits of
diversity conditioning on sample size, and the sample coverage estimate
(`SC`) along with the lower and upper confidence limits of sample
coverage (`SC.LCL`, `SC.UCL`). These sample coverage estimates with
confidence intervals are used for plotting the sample completeness
curve. It is time consuming for `diversity = FD` and `FDtype = "AUC"`.
If the argument `nboot` is greater than zero, then the bootstrap method
is applied to obtain the confidence intervals for each diversity and
sample coverage estimates.

Here only show first six rows for taxonomic diversity:

``` r
head(out.TD$TDiNextEst$size_based)
  Assemblage Order.q   m      Method       qTD   qTD.LCL   qTD.UCL         SC     SC.LCL     SC.UCL
1       Edge       0   1 Rarefaction   1.00000   1.00000   1.00000 0.01164071 0.01092485 0.01235657
2       Edge       0  95 Rarefaction  66.30595  65.37725  67.23465 0.48409405 0.47156415 0.49662395
3       Edge       0 189 Rarefaction 106.74321 104.71921 108.76721 0.63758443 0.62552422 0.64964465
4       Edge       0 284 Rarefaction 137.02920 134.01159 140.04682 0.71845672 0.70738855 0.72952488
5       Edge       0 378 Rarefaction 161.01036 157.09083 164.92989 0.76846379 0.75807630 0.77885128
6       Edge       0 472 Rarefaction 181.07340 176.29921 185.84759 0.80314380 0.79318823 0.81309937
```

The second data frame of list `$TDiNextEst` (as shown below for
‘coverage_based’) includes the sample coverage estimate (‘SC’), the
sample size (`m`), the `Method` (`Rarefaction`, `Observed`, or
`Extrapolation`, depending on whether the size `m` is less than, equal
to, or greater than the reference sample size), the diversity order
(`Order.q`), the diversity estimate of order q, the lower and upper
confidence limits of diversity conditioning on sample coverage estimate.

Here only show first six rows for taxonomic diversity:

``` r
head(out.TD$TDiNextEst$coverage_based)
  Assemblage Order.q         SC         m      Method       qTD  qTD.LCL   qTD.UCL
1       Edge       0 0.01164071   1.00000 Rarefaction   1.00000   1.0000   1.00000
2       Edge       0 0.48409405  94.99999 Rarefaction  66.30594  63.0443  69.56758
3       Edge       0 0.63758443 189.00000 Rarefaction 106.74321 101.4235 112.06289
4       Edge       0 0.71845672 284.00000 Rarefaction 137.02920 130.0249 144.03355
5       Edge       0 0.76846379 378.00002 Rarefaction 161.01037 152.4969 169.52380
6       Edge       0 0.80314380 472.00002 Rarefaction 181.07340 171.1049 191.04195
```

The output `$TDAsyEst` lists the diversity labels, the observed
diversity, asymptotic diversity estimates, estimated bootstrap standard
error (`s.e.`) and confidence intervals for diversity with q = 0, 1, and
2 (`LCL`, `UCL`). The estimated asymptotic and observed diversity can
also be computed via the function `ObsAsy3D()`. The output are shown
below:

``` r
out.TD$TDAsyEst
  Assemblage Taxonomic Diversity Taxonomic Observed Taxonomic Estimator      s.e.       LCL       UCL
1       Edge    Species richness          319.00000           444.97141 27.631000 390.81565 499.12717
2       Edge   Shannon diversity          155.38611           177.99972  3.335778 171.46171 184.53772
3       Edge   Simpson diversity           82.02345            85.90541  4.622618  76.84524  94.96557
4   Interior    Species richness          356.00000           513.51776 21.914898 470.56535 556.47017
5   Interior   Shannon diversity          163.51447           186.98301  7.714017 171.86381 202.10220
6   Interior   Simpson diversity           72.15305            74.71763  5.472373  63.99198  85.44329
```

### BASIC GRAPHIC DISPLAYS: FUNCTION ggiNEXT3D()

The function `ggiNEXT3D()`, which extends `ggplot2` with default
arguments, is described as follows:

``` r
ggiNEXT3D(outcome, type = 1:3, se = TRUE,
          facet.var = "Assemblage", color.var = "Order.q")  
```

Here `outcome` is the object of `iNEXT3D()`’s output. Three types of
curves are allowed for different diversity dimensions:

1.  Sample-size-based R/E curve (`type = 1`): This curve plots diversity
    estimates with confidence intervals as a function of sample size.

2.  Sample completeness curve (`type = 2`): This curve plots the sample
    coverage with respect to sample size.

3.  Coverage-based R/E curve (`type = 3`): This curve plots the
    diversity estimates with confidence intervals as a function of
    sample coverage.

<br><br> The argument `facet.var = "Order.q"` or
`facet.var = "Assemblage"` is used to create a separate plot for each
value of the specified variable. For example, the following code
displays a separate plot of the diversity order q. The `ggiNEXT3D()`
function is a wrapper with package `ggplot2` to create a R/E curve in a
single line of code. The figure object is of class `"ggplot"`, so can be
manipulated by using the `ggplot2` tools.

When `facet.var = "Assemblage"` in `ggiNEXT3D` function, it creates a
separate plot for each assemblage and the different color lines
represent each diversity order. Sample-size-based R/E curve (`type = 1`)
as below:

``` r
# Sample-size-based R/E curves, separating by "assemblage""
ggiNEXT3D(out.TD, type = 1, facet.var = "Assemblage")
```

<img src="README/README-unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

When `facet.var = "Order.q"` in `ggiNEXT3D` function, it creates a
separate plot for each diversity order and the different color lines
represent each assemblage. Sample-size-based R/E curve (`type = 1`) as
below:

``` r
# Sample-size-based R/E curves, separating by "Order.q"
ggiNEXT3D(out.TD, type = 1, facet.var = "Order.q")
```

<img src="README/README-unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

The following command return the sample completeness (sample coverage)
curve (`type = 2`) in which different colors are used for the three
assemblages.

``` r
ggiNEXT3D(out.TD, type = 2, facet.var = "Order.q", color.var = "Assemblage")
```

<img src="README/README-unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the coverage-based R/E sampling curves in
which different colors are used for the three assemblages
(`facet.var = "Assemblage"`) and for three diversity orders
(`facet.var = "Order.q"`).

``` r
ggiNEXT3D(out.TD, type = 3, facet.var = "Assemblage")
```

<img src="README/README-unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.TD, type = 3, facet.var = "Order.q")
```

<img src="README/README-unnamed-chunk-14-1.png" width="672" style="display: block; margin: auto;" />

### EXAMPLE for INCIDENCE-RAW DATA

Incidence raw data is allowed for three diversity dimensions. For
illustration, use the Hinkley’s fish data (the dataset is included in
the package) at three time periods (1981-1985, 1987-1991 and 2015-2019).
This data set (`fish_incidence_data`) includes a list with three
matrices; each matrix is a species-by-sampling-unit data. Here only use
taxonomic diversity (TD) as demonstration below.

``` r
data("fish_incidence_data")
head(fish_incidence_data$data$`1981-1985`)
                    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61
Agonus_cataphractus 0 0 1 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  1  0
Alosa_fallax        0 0 0 0 0 0 0 0 1  0  0  0  1  0  0  0  0  0  1  1  0  1  0  1  1  0  1  0  0  1  1  1  1  1  1  0  1  0  0  1  1  1  0  1  1  1  1  1  0  0  1  1  1  1  1  0  0  0  0  1
Ammodytes_tobianus  0 0 0 0 0 0 1 0 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
Anguilla_anguilla   1 1 1 1 0 0 0 1 0  1  0  1  0  1  1  0  0  0  1  0  0  1  0  1  1  1  1  0  0  0  1  0  1  1  0  0  1  1  0  1  0  0  0  0  1  1  1  0  1  1  1  1  1  1  0  0  0  0  0  0
Aphia_minuta        0 0 0 0 1 1 0 0 0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  1  1  1  1  1  0  1  0  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  0  1  1  1  1  1  1  0  0  0  0
Atherina_boyeri     0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0
```

``` r
# fish.tree = fish_incidence_data$tree  ## for PD
# fish.dis = fish_incidence_data$dist   ## for FD

out.raw <- iNEXT3D(data = fish_incidence_data$data, diversity = "TD",
                   q = c(0, 1, 2), datatype = "incidence_raw", nboot = 30)
ggiNEXT3D(out.raw, type = 1)
```

<img src="README/README-unnamed-chunk-16-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.raw, type = 2)
```

<img src="README/README-unnamed-chunk-17-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.raw, type = 3)
```

<img src="README/README-unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

### EXAMPLE for INCIDENCE-FREQUENCY DATA

We transform Hinkley’s fish data (incidence-raw data) to
incidence-frequency data (`incidence_freq`). The first row must be the
total number of sampling units, followed by the species incidence
frequencies in each assemblage. Here only use taxonomic diversity (TD)
for demonstration below.

``` r
fish.freq = cbind(c( ncol(fish_incidence_data$data$`1981-1985`), rowSums(fish_incidence_data$data$`1981-1985`) ),
                  c( ncol(fish_incidence_data$data$`1987-1991`), rowSums(fish_incidence_data$data$`1987-1991`) ),
                  c( ncol(fish_incidence_data$data$`2015-2019`), rowSums(fish_incidence_data$data$`2015-2019`) ))
rownames(fish.freq)[1] = "sample units"
colnames(fish.freq) = names(fish_incidence_data$data)
```

``` r
head(fish.freq)
                    1981-1985 1987-1991 2015-2019
sample units               60        60        60
Agonus_cataphractus        14        19        21
Alosa_fallax               29        22        18
Ammodytes_tobianus          4         2         8
Anguilla_anguilla          30        29         8
Aphia_minuta               22        13        18
```

Note that incidence-frequency data (`datatype = "incidence_freq`) is
allowed only for diversity class: `"TD"`, `"FD"`. If `"PD"` required,
use incidence-raw data instead. The following commands return three R/E
sampling curves for fish data. The argument `color.var =  "Order.q"` is
used to display curves in different colors for diversity order.

``` r
out.incfreq <- iNEXT3D(data = fish.freq, diversity = "TD",
                       q = c(0, 1, 2), datatype = "incidence_freq",nboot = 30)


# Sample-size-based R/E curves
ggiNEXT3D(out.incfreq, type = 1, color.var = "Order.q")
```

<img src="README/README-unnamed-chunk-21-1.png" width="672" style="display: block; margin: auto;" />

``` r
# Sample completeness curves
ggiNEXT3D(out.incfreq, type = 2)
```

<img src="README/README-unnamed-chunk-22-1.png" width="672" style="display: block; margin: auto;" />

``` r
# Coverage-based R/E curves
ggiNEXT3D(out.incfreq, type = 3, color.var = "Order.q")     
```

<img src="README/README-unnamed-chunk-23-1.png" width="672" style="display: block; margin: auto;" />

### DATA INFORMATION FUNCTION: DataInfo3D()

``` r
DataInfo3D(data, diversity = "TD", datatype = "abundance", 
           nT = NULL, PDtree, PDreftime = NULL, 
           FDdistM, FDtype = "AUC", FDtau = NULL) 
```

Here provide the function `DataInfo3D` to compute three diversity
dimensions (‘TD’, ‘PD’, ‘FD’) data information, which including sample
size, observed species richness, sample coverage estimate, and the first
ten abundance/incidence frequency counts when `diversity = TD`. And so
on for PD, FD.

``` r
DataInfo3D(Brazil_rainforest_data$data, diversity = 'TD', datatype = "abundance")
  Assemblage    n S.obs     SC(n)    SC(2n)  f1 f2 f3 f4 f5
1       Edge 1794   319 0.9387143 0.9744012 110 48 38 28 13
2   Interior 2074   356 0.9407166 0.9728436 123 48 41 32 19
```

### POINT ESTIMATION FUNCTION: estimate3D()

``` r
estimate3D(data, diversity = "TD", q = c(0, 1, 2), datatype = "abundance", 
           base = "coverage", level = NULL, nboot = 50, conf = 0.95, 
           nT = NULL, PDtree, PDreftime = NULL, PDtype = "meanPD", 
           FDdistM, FDtype = "AUC", FDtau = NULL)
```

`estimate3D` is used to compute three diversity dimensions (TD, PD, FD)
estimates with q = 0, 1, 2 under any specified level of sample size
(when `base = "size"`) or sample coverage (when `base = "coverage"`) for
either abundance data (`datatype = "abundance"`) or incidence data
(`datatype = "incidence_freq"` or `"incidence_raw"`). If `level = NULL`,
this function computes the diversity estimates for the minimum sample
size among all samples extrapolated to double reference sizes (when
`base = "size"`) or the minimum sample coverage among all samples
extrapolated to double reference sizes (when `base = "coverage"`).

For example, the following command returns the taxonomic diversity
(‘TD’) with a specified level of sample coverage = 99.5% for the Brazil
rainforest data. For some assemblages, this coverage value corresponds
to the rarefaction part whereas the others correspond to extrapolation.

``` r
estimate3D(Brazil_rainforest_data$data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", 
           base = "coverage", level = 0.995)
  Assemblage Order.q    SC        m        Method       qTD      s.e.   qTD.LCL   qTD.UCL
1       Edge       0 0.995 6944.000 Extrapolation 434.69401 25.289606 385.12730 484.26073
2       Edge       1 0.995 6944.000 Extrapolation 176.09987  5.341759 165.63022 186.56953
3       Edge       2 0.995 6944.000 Extrapolation  84.86771  4.620648  75.81141  93.92402
4   Interior       0 0.995 8643.312 Extrapolation 500.23261 23.538804 454.09740 546.36782
5   Interior       1 0.995 8643.312 Extrapolation 185.36031  5.760025 174.07087 196.64975
6   Interior       2 0.995 8643.312 Extrapolation  74.08576  4.239065  65.77735  82.39418
```

### EMPIRICAL AND ASYMPTOTIC DIVERSITY FUNCTION: ObsAsy3D

``` r
ObsAsy3D(data,diversity = "TD",q = seq(0,2,0.2),datatype = "abundance",
         nboot = 50,conf = 0.95,nT = NULL,method = c("Asymptotic", "Observed"),
         PDtree,PDreftime = NULL,PDtype = "meanPD",
         FDdistM,FDtype = "AUC",FDtau = NULL)
```

The function `ObsAsy3D()` compute three diversity dimensions (TD, PD,
FD) for empirical (observed) diversity and estimated asymptotic
diversity with any diversity order. For example, the following commands
returns empirical and asymptotic taxonomic diversity (‘TD’) for Brazil
rainforest data, along with its confidence interval at diversity order q
from 0 to 2. Here only show the first ten rows.

``` r
out1 <- ObsAsy3D(Brazil_rainforest_data$data, diversity = 'TD', datatype = "abundance", 
                 method = c("Asymptotic", "Observed"), nboot = 30, conf = 0.95)

head(out1, 10)
   Assemblage Order.q       qTD      s.e.   qTD.LCL  qTD.UCL     Method
1        Edge     0.0 444.97141 34.206345 377.92821 512.0146 Asymptotic
2        Edge     0.2 375.27048 23.567431 329.07917 421.4618 Asymptotic
3        Edge     0.4 312.45184 15.615533 281.84595 343.0577 Asymptotic
4        Edge     0.6 258.37946 10.451587 237.89472 278.8642 Asymptotic
5        Edge     0.8 213.72995  7.541398 198.94909 228.5108 Asymptotic
6        Edge     1.0 177.99972  5.984010 166.27127 189.7282 Asymptotic
7        Edge     1.2 149.91442  5.078316 139.96111 159.8677 Asymptotic
8        Edge     1.4 127.94456  4.514335 119.09663 136.7925 Asymptotic
9        Edge     1.6 110.67216  4.183495 102.47266 118.8717 Asymptotic
10       Edge     1.8  96.94807  4.020864  89.06732 104.8288 Asymptotic
```

### GRAPHIC DISPLAYS FUNCTION: ggObsAsy3D()

Plots q-profile, time-profile, and tau-profile based on the outcome of
ObsAsy3D using the ggplot2 package.

The function ggObsAsy3D(), which extends ggplot2 with default arguments,
is described as follows:

``` r
ggObsAsy3D(outcome, profile = "q")
```

`ggObsAsy3D` plots q-profile, time-profile, and tau-profile based on
`ggplot2`. Here `outcome` is the object from the function `ObsAsy3D`,
and `profile` is a profile selection versus to diversity. Default is
`profile = "q"`. Note that `profile = "time"` is allowed for only when
diversity = “PD” and `profile = "tau"` profile is allowed for only when
`diversity = "FD"` and `FDtype = "tau_values"`.

``` r
# q-profile curves
ggObsAsy3D(out1,profile = "q")
```

<img src="README/README-unnamed-chunk-31-1.png" width="672" style="display: block; margin: auto;" />

The argument `profile = "time"` in `ggObsAsy3D` function creates a
separate plot for each diversity order q. Therefore the different
assemblages will be represented by different color lines.

``` r
# time-profile curves, separating by "Order.q"
data(fish_incidence_data)
data <- fish_incidence_data$data
tree <- fish_incidence_data$tree

out2 <- ObsAsy3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", 
                 nboot = 30, method = c("Asymptotic", "Observed"), 
                 PDtree = tree, PDreftime = seq(0.1, 82.8575, length.out = 40))
ggObsAsy3D(out2, profile = "time")
```

<img src="README/README-unnamed-chunk-32-1.png" width="672" style="display: block; margin: auto;" />

The argument `profile = "tau"` in `ggObsAsy3D` function creates a
separate plot for each diversity order q. Therefore the different
assemblages will be represented by different color lines.

``` r
# tau-profile curves, separating by "Order.q"
data <- Brazil_rainforest_data$data
distM <- Brazil_rainforest_data$dist

out3 <- ObsAsy3D(data, diversity = 'FD', q = c(0, 1, 2), datatype = "abundance", nboot = 30, 
                 method = c("Asymptotic", "Observed"), FDtau = seq(0, 1, 0.1), FDdistM = distM, FDtype = 'tau_values')
ggObsAsy3D(out3, profile = "tau")
```

<img src="README/README-unnamed-chunk-33-1.png" width="672" style="display: block; margin: auto;" />

### How to cite

If you publish your work based on results from iNEXT.3D package, you
should make references to the following Online reference:

-   Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H.,
    Dornelas, M and. Magurran, A. E. (2021). Measuring temporal change
    in alpha diversity: a framework integrating taxonomic, phylogenetic
    and functional diversity and the iNEXT.3D standardization. Methods
    in Ecology and Evolution, 12, 1926-1940.

### License

The iNEXT.3D package is licensed under the GPLv3. To help refine
`iNEXT.3D`, your comments or feedback would be welcome (please send them
to Anne Chao or report an issue on the iNEXT.3D github
[iNEXT.3D_github](https://github.com/AnneChao/iNEXT.3D).

### References

-   Chao, A., Chiu, C.-H., Villéger, S., Sun, I.-F., Thorn, S., Lin,
    Y.-C., Chiang, J. M. and Sherwin, W. B. (2019). An
    attribute-diversity approach to functional diversity, functional
    beta diversity, and related (dis)similarity measures. Ecological
    Monographs, 89, e01343. 10.1002/ecm.1343.

-   Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H.,
    Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in
    alpha diversity: a framework integrating taxonomic, phylogenetic and
    functional diversity and the iNEXT.3D standardization. Methods in
    Ecology and Evolution, 12, 1926-1940.

-   T.C. Hsieh, K. H. Ma, and Chao, A. (2016). iNEXT: An R package for
    rarefaction and extrapolation of species diversity (Hill numbers).
    Methods in Ecology and Evolution, 7, 1451-1456.
