<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.3D (R package)

<h5 align="right">
Latest version: 2023-04-10
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
showing three types of sampling curves. See Chao et al. (2021) for
methodologies. An online version of [iNEXT.3D
Online](https://chao.shinyapps.io/iNEXT_3D/) is also available for users
without an R background. Detailed information about all functions in
iNEXT.3D is provided in the iNEXT.3D Manual in
[iNEXT.3D_vignettes](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/A%20Quick%20Introduction%20to%20iNEXT.3D%20via%20Examples.html),which
is available from [Anne Chao’s
website](http://chao.stat.nthu.edu.tw/wordpress/software_download/).

`iNEXT3D` focuses on three measures of Hill numbers of order q: species
richness (`q = 0`), Shannon diversity (`q = 1`, the exponential of
Shannon entropy) and Simpson diversity (`q = 2`, the inverse of Simpson
concentration) among three classes: taxonomic diversity (‘TD’),
phylogenetic diversity (‘PD’), functional diversity (‘FD’). For each
diversity measure, `iNEXT.3D` uses the observed sample of abundance or
incidence data (called the “reference sample”) to compute diversity
estimates and the associated confidence intervals for the following two
types of rarefaction and extrapolation (R/E):

1.  Sample‐size‐based R/E sampling curves:`iNEXT3D` computes diversity
    estimates for rarefied and extrapolated samples up to an appropriate
    size. This type of sampling curve plots the diversity estimates with
    respect to sample size.  
2.  Coverage‐based R/E sampling curves: `iNEXT3D` computes diversity
    estimates for rarefied and extrapolated samples with sample
    completeness (as measured by sample coverage) up to an appropriate
    coverage. This type of sampling curve plots the diversity estimates
    with respect to sample coverage.

`iNEXT3D` also plots the above two types of sampling curves and a sample
completeness curve. The sample completeness curve provides a bridge
between these two types of curves.

### SOFTWARE NEEDED TO RUN INEXT.3D IN R

-   Required: [R](http://cran.rstudio.com/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### HOW TO RUN INEXT.3D:

The iNEXT.3D package is available on
[Github](https://github.com/AnneChao/iNEXT.3D) and can be downloaded
with a standard installation procedure using the commands shown
below.For a first-time installation, an additional visualization
extension package (ggplot2) must be installed and loaded.

``` r
## install the latest version from github
install.packages('devtools')
library(devtools)
# install_github('AnneChao/iNEXT.3D')
## import packages
library(iNEXT.3D)
library(ggplot2)
```

In this document, we provide a quick introduction demonstrating how to
run the package `iNEXT.3D`(iNterpolation and EXTrapolation in three
Dimensions). `iNEXT.3D` has several main function: `iNEXT3D`,
`ggiNEXT3D`, `AO3D`, `ggAO3D`, `estimate3D`, and `DataInfo3D.`

### MAIN FUNCTION: iNEXT3D()

The main function iNEXT3D() with default arguments is described below:
<br><br> iNEXT3D(data, diversity = ‘TD’, q = c(0,1,2), datatype =
“abundance”, size = NULL, endpoint = NULL, knots = 40, nboot = 50, conf
= 0.95, nT = NULL, PDtree = NULL, PDreftime = NULL, PDtype = ‘meanPD’,
FDdistM, FDtype = ‘AUC’, FDtau = NULL) <br><br> This main function
computes diversity estimates of order q, the sample coverage estimates
and related statistics for K (if knots=K) evenly-spaced knots (sample
sizes) between size 1 and the endpoint, where the endpoint is as
described below. Each knot represents a particular sample size for which
diversity estimates will be calculated. By default, endpoint is set to
be double the reference sample size.

This function returns an “iNEXT3D” object which can be further used to
make plots using the function ggiNEXT3D() to be described below.

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
make plots using the function ggiNEXT3D() to be described below.

### DATA FORMAT/INFORMATION

Three types of data are supported: (“abundance”, “incidence_raw”, or
“incidence_freq”):

1.  Individual‐based abundance data (datatype=“abundance”): Input data
    for each assemblage/site include samples species abundances in an
    empirical sample of n individuals (“reference sample”). When there
    are N assemblages, input data consist of an S by N abundance matrix,
    or N lists of species abundances.

2.  Sampling-unit-based incidence data: There are two kinds of input
    data.

(2a) Incidence‐raw data (datatype=“incidence_raw”): for each assemblage,
input data for a reference sample consist of a species‐by‐sampling‐unit
matrix; when there are N assemblages, input data consist of N lists of
matrices, and each matrix is a species‐by‐sampling‐unit matrix. If the
phylogenetic diversity required, the matrix of combined assemblage is
allowed, but nT must be specified (see above description).

(2b) Incidence‐frequency data (datatype=“incidence_freq”): input data
for each assemblage consist of species sample incidence frequencies (row
sums of each incidence matrix). When there are N assemblages, input data
consist of an S by N matrix, or N lists of species incidence
frequencies. The first entry of each list must be the total number of
sampling units, followed by the species incidence frequencies.

Two data sets (beetles for abundance data ) are included in iNEXT3D
package. The first list of beetles’ data consist of abundance data who
live on three treatments (“Control”, “Debarked” and “Scratched”) of
trees. The second list consist of the pylogenetic tree for every
beetles. And the third list consist of distance matrix for each pair of
beetles. For these data, the following commands display the sample
species abundances and run the iNEXT3D() function for three types of
diversty (“TD”, “PD”, “FD” with threshold dmean, “AUC” which integate FD
from threshold 0-1) in q = 0.

``` r
data("dunes")

out.TD <- iNEXT3D(data = dunes$data, diversity = "TD", 
               q = c(0, 1, 2), datatype = "abundance", 
               nboot = 10)
out.PD <- iNEXT3D(data = dunes$data, diversity = "PD", 
               q = c(0, 1, 2), datatype = "abundance", 
               PDtree = dunes$tree, 
               nboot = 10)
out.FD <- iNEXT3D(data = dunes$data, diversity = "FD", 
               q = c(0, 1, 2), datatype = "abundance", 
               FDdistM = dunes$dist,
               nboot = 5)
```

``` r
out.TD
Compare 3 assemblages with Hill number order q = 0, 1, 2.
$class: iNEXT3D

DataInfo: basic data information
  Assemblage    n S.obs     SC f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
1         EM  373    17 0.9920  3  1  1  0  0  1  0  1  1   1
2         MO 1490    39 0.9987  2  1  1  4  2  2  3  1  1   1
3         TR 1059    42 0.9962  4  2  4  2  2  2  0  1  2   2

iNextEst: diversity estimates with rarefied and extrapolated samples.
$size_based (LCL and UCL are obtained for fixed size.)

# A tibble: 45 x 10
   Assemblage     m Method        Order.q    qD qD.LCL qD.UCL    SC SC.LCL SC.UCL
   <chr>      <dbl> <chr>           <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl>
 1 EM             1 Rarefaction         0  1      1.00   1    0.136  0.128  0.145
 2 EM           196 Rarefaction         0 15.2   13.7   16.7  0.987  0.980  0.994
 3 EM           373 Observed            0 17     14.7   19.3  0.992  0.988  0.996
 4 EM           550 Extrapolation       0 18.2   15.3   21.2  0.994  0.990  0.998
 5 EM           746 Extrapolation       0 19.2   15.6   22.8  0.996  0.993  0.999
 6 EM             1 Rarefaction         1  1      1      1    0.136  0.128  0.145
 7 EM           196 Rarefaction         1  9.07   8.62   9.51 0.987  0.980  0.994
 8 EM           373 Observed            1  9.27   8.80   9.75 0.992  0.988  0.996
 9 EM           550 Extrapolation       1  9.37   8.88   9.85 0.994  0.990  0.998
10 EM           746 Extrapolation       1  9.43   8.93   9.93 0.996  0.993  0.999
# ... with 35 more rows
# i Use `print(n = ...)` to see more rows

NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$size_based to view complete output.

$coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)

# A tibble: 45 x 8
   Assemblage    SC     m Method        Order.q    qD qD.LCL qD.UCL
   <chr>      <dbl> <dbl> <chr>           <dbl> <dbl>  <dbl>  <dbl>
 1 EM         0.136    1  Rarefaction         0  1     0.967   1.03
 2 EM         0.987  196. Rarefaction         0 15.2  12.1    18.3 
 3 EM         0.992  373  Observed            0 17    12.6    21.4 
 4 EM         0.994  550. Extrapolation       0 18.2  13.3    23.1 
 5 EM         0.996  746. Extrapolation       0 19.2  13.8    24.5 
 6 EM         0.136    1  Rarefaction         1  1     0.968   1.03
 7 EM         0.987  196. Rarefaction         1  9.07  8.48    9.66
 8 EM         0.992  373  Observed            1  9.27  8.70    9.84
 9 EM         0.994  550. Extrapolation       1  9.37  8.80    9.93
10 EM         0.996  746. Extrapolation       1  9.43  8.87    9.99
# ... with 35 more rows
# i Use `print(n = ...)` to see more rows

NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$coverage_based to view complete output.

AsyEst: asymptotic diversity estimates along with related statistics.
  Assemblage         Diversity  Observed Estimator      s.e.       LCL       UCL
1         EM  Species richness 17.000000 21.487936 7.8075998 17.000000 36.790550
2         EM Shannon diversity  9.272102  9.541765 0.2735255  9.005665 10.077865
3         EM Simpson diversity  7.218106  7.340810 0.3037283  6.745514  7.936107
4         MO  Species richness 39.000000 40.998658 3.8025660 39.000000 48.451550
5         MO Shannon diversity 19.712107 19.987465 0.4717668 19.062819 20.912111
6         MO Simpson diversity 13.979950 14.102888 0.4343910 13.251498 14.954279
7         TR  Species richness 42.000000 45.996223 3.9371403 42.000000 53.712876
8         TR Shannon diversity 23.941938 24.481173 0.4528253 23.593652 25.368695
9         TR Simpson diversity 18.605455 18.920295 0.4792300 17.981021 19.859568
```

The second data frame of list $iNextEst (as shown below for
‘coverage_based’) under each class includes the reference sample
coverage (‘goalSC’), real sample coverage (‘SC’), sample size (m, i.e.,
each of the 40 knots), the method (Rarefaction, Observed, or
Extrapolation, depending on whether the size m is less than, equal to,
or greater than the reference sample size), the diversity order, the
diversity estimate of order q conditioning on ‘goalSC’, the 95% lower
and upper confidence limits of diversity. These sample coverage
estimates with confidence intervals are used for plotting the
coverage-based R/E curves.

$AsyEst lists the observed diversity, asymptotic estimates, estimated
bootstrap s.e. and 95% confidence intervals for Hill numbers with q = 0,
1, and 2. The estimated asymptotes and the observed diversity are
calculated via the functions AO3D() . The output for the dunes data is
shown below. All row and column variables are self‐explanatory.

### BASIC GRAPHIC DISPLAYS: FUNCTION ggiNEXT3D()

The ggiNEXT3D() function, which extends ggplot2 to the “iNEXT3D” object,
is described as follows with default arguments: <br><br>
ggiNEXT3D(outcome, type = 1:3, se = TRUE, facet.var = “Assemblage”,
color.var = “Order.q”)  
<br><br> Here outcome is the object of iNEXT3D()’s output. Three types
of curves are allowed for different diversity classes:

1.  Sample-size-based R/E curve (type=1): see Figs. 1a and 2a in the
    main text. This curve plots diversity estimates with confidence
    intervals (if se=TRUE) as a function of sample size up to double the
    reference sample size, by default, or a user‐specified endpoint.

2.  Sample completeness curve (type=2) with confidence intervals (if
    se=TRUE): see Figs. 1b and 2b in the main text. This curve plots the
    sample coverage with respect to sample size for the same range
    described in (1).

3.  Coverage-based R/E curve (type=3): see Figs. 1c and 2c in the main
    text. This curve plots the diversity estimates with confidence
    intervals (if se=TRUE) as a function of sample coverage up to the
    maximum coverage obtained from the maximum size described in (1).

<br><br> The argument facet.var=(“Order.q”, “Assemblage”) is used to
create a separate plot for each value of the specified variable. For
example, the following code displays a separate plot (in Figs 1a and 1c)
for each value of the diversity order q. The ggiNEXT3D() function is a
wrapper around ggplot2 package to create a R/E curve using a single line
of code. The resulting object is of class “ggplot”, so can be
manipulated using the ggplot2 tools.

The argument facet.var=“Assemblage” in ggiNEXT3D function creates a
separate plot for each assembalge, therefore the different Order.q will
seperated by different colours as shown below:

### RAREFACTION/EXTRAPOLATION FOR ABUNDANCE DATA

``` r
# Sample‐size‐based R/E curves, separating by "assemblage""
ggiNEXT3D(out.TD, type = 1, facet.var = "Assemblage")
```

<img src="README/README-unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

The argument facet.var=“Order.q” in ggiNEXT3D function creates a
separate plot for each order, therefore three assemblages will be
seperated by different colours as shown below:

``` r
# Sample-size-based R/E curves, separating plots by "Order.q"
ggiNEXT3D(out.TD, type = 1, facet.var = "Order.q")
```

<img src="README/README-unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the sample completeness curve in which
different colors are used for the three assemblages. Since the sample
completeness curve are same for differnet class of diversity, ggiNEXT3D
returns only one plot:

``` r
ggiNEXT3D(out.TD, type = 2, facet.var = "Order.q", color.var="Assemblage")
```

<img src="README/README-unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the coverage‐based R/E sampling curves in
which different colors are used for the three assemblages
(facet.var=“Assemblage”) and for three orders (facet.var=“Order.q”)

``` r
ggiNEXT3D(out.TD, type = 3, facet.var="Assemblage")
```

<img src="README/README-unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.TD, type = 3, facet.var="Order.q")
```

<img src="README/README-unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

### RAREFACTION/EXTRAPOLATION FOR RAW INCIDENCE DATA: (incidence_raw)

For illustration, we use the Hinkley’s fish data (in the dataset fish
included in the package) at three time periods (1981-1985, 1987-1991 and
2015-2019). Incidence raw data is allowed for all diversity class.This
data set (fish) included in the package is three list of matrices; each
matrix is a species by plots data.frame for a time period. We only use
Taxonomic diversity(TD) for demonstration below.

``` r
data("fish")
head(fish[[1]]$`1981-1985`,4)
                    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61
Agonus_cataphractus 0 0 1 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  1  0
Alosa_fallax        0 0 0 0 0 0 0 0 1  0  0  0  1  0  0  0  0  0  1  1  0  1  0  1  1  0  1  0  0  1  1  1  1  1  1  0  1  0  0  1  1  1  0  1  1  1  1  1  0  0  1  1  1  1  1  0  0  0  0  1
Ammodytes_marinus   0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
Ammodytes_tobianus  0 0 0 0 0 0 1 0 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
```

``` r
# fish.tree = fish$tree
# fish.dis = fish$dist

out.raw <- iNEXT3D(data = fish$data, 
                   diversity = "TD",
                   q = c(0, 1, 2), datatype = "incidence_raw",nboot = 10)
ggiNEXT3D(out.raw, type = 1)
```

<img src="README/README-unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.raw, type = 2)
```

<img src="README/README-unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT3D(out.raw, type = 3)
```

<img src="README/README-unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

### DATA INFORMATION FUNCTION: DataInfo3D()

We can supply the function <br><br> DataInfo3D(data, diversity = “TD”,
datatype = “abundance”, nT = NULL, PDtree, PDreftime = NULL, FDdistM,
FDtype = “AUC”, FDtau = NULL) <br><br> to compute three type
diversity(‘TD’,‘PD’,‘FD’) data information, which including sample size,
observed species richness, sample coverage estimate, and the first ten
abundance/incidence frequency counts, and so on.

``` r
DataInfo3D(dunes$data, diversity = 'TD', datatype = "abundance")
  Assemblage    n S.obs     SC f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
1         EM  373    17 0.9920  3  1  1  0  0  1  0  1  1   1
2         MO 1490    39 0.9987  2  1  1  4  2  2  3  1  1   1
3         TR 1059    42 0.9962  4  2  4  2  2  2  0  1  2   2
```

### POINT ESTIMATION FUNCTION: estimate3D()

We also supply the function <br><br> estimate3D(data, diversity = “TD”,
q = c(0, 1, 2), datatype = “abundance”, base = “coverage”, level = NULL,
nboot = 50, conf = 0.95, nT = NULL, PDtree, PDreftime = NULL, PDtype =
“meanPD”, FDdistM, FDtype = “AUC”, FDtau = NULL) <br><br> to compute
three type diversity(‘TD’,‘PD’,‘FD’) estimates with q = 0, 1, 2 for any
particular level of sample size (base=“size”) or any specified level of
sample coverage (base=“coverage”) for either abundance data
(datatype=“abundance”) or incidence data (datatype=“incidence_freq” or
“incidence_raw”). If level=NULL, this function computes the diversity
estimates for the minimum sample size/coverage among all assemblages.

For example, the following command returns the taxonomic diversity
(‘TD’) with a specified level of sample coverage of 93% for the dunes
data. For some assemblages, this coverage value corresponds to the
rarefaction part whereas the others correspond to extrapolation, as
indicated in the method of the output.

``` r
estimate3D(dunes$data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", base = "coverage",level = 0.93)
  Assemblage   SC         m      Method Order.q        qD      s.e.    qD.LCL    qD.UCL
1         EM 0.93  47.88966 Rarefaction       0 11.194508 0.5723628 10.072697 12.316318
2         EM 0.93  47.88966 Rarefaction       1  8.004329 0.3413100  7.335374  8.673284
3         EM 0.93  47.88966 Rarefaction       2  6.482469 0.2846312  5.924602  7.040336
4         MO 0.93 107.56810 Rarefaction       0 24.830302 0.8163906 23.230206 26.430398
5         MO 0.93 107.56810 Rarefaction       1 16.735683 0.5061083 15.743729 17.727637
6         MO 0.93 107.56810 Rarefaction       2 12.571520 0.4233112 11.741846 13.401195
7         TR 0.93 114.37603 Rarefaction       0 27.935556 0.8959244 26.179577 29.691536
8         TR 0.93 114.37603 Rarefaction       1 20.360646 0.5638913 19.255439 21.465853
9         TR 0.93 114.37603 Rarefaction       2 16.357398 0.4804488 15.415736 17.299060
```

### EMPIRICAL AND ASYMPTOTIC DIVERSITY FUNCTION: AO3D

<br><br> AO3D( data,diversity = “TD”,q = seq(0,2,0.2),datatype =
“abundance”, nboot = 50,conf = 0.95,nT = NULL,method = c(“Asymptotic”,
“Observed”), PDtree,PDreftime = NULL,PDtype = “meanPD”, FDdistM,FDtype =
“AUC”,FDtau = NULL ) <br><br>

The function AO3D() can compute three type
diversity(‘TD’,‘PD’,‘FD’),which including empirical diversity and
asymptotic diversity. For either abundance data (datatype=“abundance”)
or incidence data(datatype=“incidence_freq” or “incidence_raw”)with any
specified level of q can be compute.

For example, the following command returns an empirical taxonomic
diversity(‘TD’) and asymptotic taxonomic diversity(‘TD’) for dunes data,
along with its confidence interval, for a specified q level from 0 to 2.

``` r
out1 <- AO3D(dunes$data, diversity = 'TD', datatype = "abundance", method = c("Asymptotic", "Observed"),nboot = 5,conf = 0.95)

out1
   Order.q        qD       s.e.    qD.LCL    qD.UCL Assemblage     Method
1      0.0 21.487936  6.2179769  9.300925 33.674946         EM Asymptotic
2      0.2 16.936433  2.9506297 11.153305 22.719561         EM Asymptotic
3      0.4 13.885407  1.2425274 11.450098 16.320716         EM Asymptotic
4      0.6 11.862781  0.5783225 10.729290 12.996272         EM Asymptotic
5      0.8 10.496859  0.4821520  9.551859 11.441860         EM Asymptotic
6      1.0  9.541765  0.5042795  8.553395 10.530135         EM Asymptotic
7      1.2  8.847719  0.5139625  7.840371  9.855067         EM Asymptotic
8      1.4  8.325449  0.5099681  7.325930  9.324968         EM Asymptotic
9      1.6  7.920890  0.5002279  6.940461  8.901319         EM Asymptotic
10     1.8  7.600095  0.4890426  6.641589  8.558601         EM Asymptotic
11     2.0  7.340810  0.4782234  6.403510  8.278111         EM Asymptotic
12     0.0 40.998658  2.7481544 35.612374 46.384941         MO Asymptotic
13     0.2 34.209437  0.9715199 32.305293 36.113581         MO Asymptotic
14     0.4 29.060182  0.3577466 28.359011 29.761352         MO Asymptotic
15     0.6 25.175831  0.4467931 24.300133 26.051530         MO Asymptotic
16     0.8 22.235885  0.5318755 21.193429 23.278342         MO Asymptotic
17     1.0 19.987465  0.5576188 18.894553 21.080378         MO Asymptotic
18     1.2 18.241913  0.5565624 17.151071 19.332756         MO Asymptotic
19     1.4 16.863382  0.5444817 15.796217 17.930546         MO Asymptotic
20     1.6 15.755819  0.5275033 14.721931 16.789706         MO Asymptotic
21     1.8 14.851614  0.5080952 13.855766 15.847463         MO Asymptotic
22     2.0 14.102888  0.4874604 13.147484 15.058293         MO Asymptotic
23     0.0 45.996223 17.2515133 12.183878 79.808568         TR Asymptotic
24     0.2 38.731480  7.6521019 23.733636 53.729324         TR Asymptotic
25     0.4 33.399637  3.0887934 27.345714 39.453561         TR Asymptotic
26     0.6 29.505773  1.1925180 27.168481 31.843065         TR Asymptotic
27     0.8 26.635682  0.4662209 25.721906 27.549458         TR Asymptotic
28     1.0 24.481173  0.2105288 24.068544 24.893802         TR Asymptotic
29     1.2 22.827199  0.1566725 22.520126 23.134271         TR Asymptotic
30     1.4 21.527891  0.1836786 21.167888 21.887895         TR Asymptotic
31     1.6 20.484834  0.2283208 20.037333 20.932334         TR Asymptotic
32     1.8 19.631065  0.2724989 19.096977 20.165153         TR Asymptotic
33     2.0 18.920295  0.3118545 18.309071 19.531518         TR Asymptotic
34     0.0 17.000000  1.5811388 13.901025 20.098975         EM  Empirical
35     0.2 14.471457  1.0255475 12.461421 16.481493         EM  Empirical
36     0.4 12.567736  0.6403398 11.312693 13.822779         EM  Empirical
37     0.6 11.145510  0.3930937 10.375060 11.915959         EM  Empirical
38     0.8 10.079552  0.2515846  9.586456 10.572649         EM  Empirical
39     1.0  9.272102  0.1871846  8.905227  9.638977         EM  Empirical
40     1.2  8.651511  0.1684919  8.321272  8.981749         EM  Empirical
41     1.4  8.166818  0.1672663  7.838982  8.494654         EM  Empirical
42     1.6  7.782094  0.1693320  7.450209  8.113979         EM  Empirical
43     1.8  7.471913  0.1704998  7.137740  7.806087         EM  Empirical
44     2.0  7.218106  0.1702588  6.884405  7.551808         EM  Empirical
45     0.0 39.000000  1.0000000 37.040036 40.959964         MO  Empirical
46     0.2 32.970780  0.6064201 31.782218 34.159341         MO  Empirical
47     0.4 28.269000  0.5268270 27.236438 29.301562         MO  Empirical
48     0.6 24.647500  0.5252337 23.618060 25.676939         MO  Empirical
49     0.8 21.864493  0.5266566 20.832265 22.896721         MO  Empirical
50     1.0 19.712107  0.5275484 18.678131 20.746083         MO  Empirical
51     1.2 18.026804  0.5304494 16.987142 19.066465         MO  Empirical
52     1.4 16.686924  0.5348825 15.638574 17.735275         MO  Empirical
53     1.6 15.604656  0.5392157 14.547812 16.661499         MO  Empirical
54     1.8 14.717297  0.5421461 13.654711 15.779884         MO  Empirical
55     2.0 13.979950  0.5430364 12.915619 15.044282         MO  Empirical
56     0.0 42.000000  2.4899799 37.119729 46.880271         TR  Empirical
57     0.2 36.253363  1.8114448 32.702996 39.803730         TR  Empirical
58     0.4 31.844736  1.4270242 29.047820 34.641652         TR  Empirical
59     0.6 28.487492  1.2214653 26.093464 30.881520         TR  Empirical
60     0.8 25.922807  1.1066850 23.753745 28.091870         TR  Empirical
61     1.0 23.941938  1.0343156 21.914717 25.969159         TR  Empirical
62     1.2 22.387845  0.9827904 20.461611 24.314078         TR  Empirical
63     1.4 21.146894  0.9430966 19.298459 22.995330         TR  Empirical
64     1.6 20.138290  0.9111937 18.352383 21.924197         TR  Empirical
65     1.8 19.304760  0.8849557 17.570279 21.039241         TR  Empirical
66     2.0 18.605455  0.8630694 16.913870 20.297040         TR  Empirical
```

### GRAPHIC DISPLAYS FUNCTION: ggAO3D()

Plots q-profile, time-profile, and tau-profile based on the outcome of
AO3D using the ggplot2 package.

The function ggAO3D(), which extends ggplot2 with default arguments, is
described as follows: <br><br> ggAO3D(outcome, profile = “q”) <br><br>
Here outcome is the object of AO3D’s output and profile is an object
versus to diversity.Default is profile = “q”.Note that profile = “time”
is allowed for only when diversity = “PD” and profile = “tau” profile is
allowed for only when diversity = “FD” and FDtype = “tau_values”.

``` r
# q profile curve""
ggAO3D(out1,profile = "q")
```

<img src="README/README-unnamed-chunk-17-1.png" width="672" style="display: block; margin: auto;" />

The argument profile = “time” in ggAO3D function creates a separate plot
for each order.q. Therefore the different assemblages will seperated by
different colours as shown below:

``` r
# time curves, separating by "Order.q"
data(data.inc)
data <- data.inc$data
tree <- data.inc$tree
nT <- data.inc$nT
out2 <- AO3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", nT = nT, nboot = 5, method = c("Asymptotic", "Observed"), PDtree = tree, PDreftime = seq(0.1, 82.8575, length.out = 40))
ggAO3D(out2, profile = "time")
```

<img src="README/README-unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

The argument profile = “tau” in ggAO3D function creates a separate plot
for each order.q. Therefore the different assemblages will seperated by
different colours as shown below:

``` r
# tau curves, separating by "Order.q"
data <- dunes$data
distM <-  dunes$dist
out3 <- AO3D(data, diversity = 'FD', q = c(0, 1, 2), datatype = "abundance", nboot = 5, method = c("Asymptotic", "Observed"), FDtau = seq(0, 1, 0.1), FDdistM = distM, FDtype = 'tau_values')
ggAO3D(out3, profile = "tau")
```

<img src="README/README-unnamed-chunk-19-1.png" width="672" style="display: block; margin: auto;" />

### How to cite

-   Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H.,
    Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in
    alpha diversity: a framework integrating taxonomic, phylogenetic and
    functional diversity and the iNEXT.3D standardization. Methods in
    Ecology and Evolution, 12, 1926-1940.

### Referance

Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas,
M and Magurran, A. E. (2021). Measuring temporal change in alpha
diversity: a framework integrating taxonomic, phylogenetic and
functional diversity and the iNEXT.3D standardization. Methods in
Ecology and Evolution, 12, 1926-1940.
