\name{iNEXT.3D-package}
\alias{iNEXT.3D-package}
\alias{iNEXT.3D}
\docType{package}
\title{
iNterpolation and EXTrapolation in three Dimensions
}
\description{
\code{iNEXT.3D} is the extension of R package \code{iNEXT} (Hsieh et al., 2016). \code{iNEXT.3D} focuses on three measures of Hill numbers of order q: species richness (\code{q = 0}), Shannon diversity (\code{q = 1}, the exponential of Shannon entropy) and Simpson diversity (\code{q = 2}, the inverse of Simpson concentration) and extend Hill number to three dimensions: taxonomic diversity (TD), phylogenetic diversity (PD), and functional diversity (FD) under Hill-Chao family frame work (Chao et al., 2019). Besides, \code{iNEXT.3D} also provide statistic estimation for three dimensions biodiversity (Chao et al., 2021). For each diversity measure, \code{iNEXT.3D} uses the observed sample of abundance or incidence data (called the "reference sample") to compute diversity estimates and the associated confidence intervals for the following two types of rarefaction and extrapolation (R/E): 

* Sample-size-based R/E sampling curves:\code{iNEXT.3D} computes diversity estimates for rarefied and extrapolated samples up to an appropriate size. This type of sampling curve plots the diversity estimates with respect to sample size.

* Coverage-based R/E sampling curves: \code{iNEXT.3D} computes diversity estimates for rarefied and extrapolated samples with sample completeness (as measured by sample coverage) up to an appropriate coverage. This type of sampling curve plots the diversity estimates with respect to sample coverage. 

\code{iNEXT.3D} also plots the above two types of sampling curves and a sample completeness curve by `ggiNEXT3D`. The sample completeness curve provides a bridge between these two types of curves.

This package contains six main functions: \cr\cr
1. \code{iNEXT3D} Computing interpolation and extrapolation of Hill number with order q for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage or sample size. \cr\cr
2. \code{ggiNEXT3D} Visualizing the output from the function \code{iNEXT.3D}. \cr\cr
3. \code{AO3D} Computing empirical and asymptotic diversity of Hill number with order q for three dimensions: taxonomic, phylogenetic and functional diversity. \cr\cr
4. \code{ggAO3D} Visualizing the output from the function `AO3D`. \cr\cr
5. \code{estimate3D} Computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage. \cr\cr
6. \code{DataInfo3D} Exhibit basic data information \cr\cr

An Online version of \code{iNEXT.3D} is also available for users without an R background: \cr
https://chao.shinyapps.io/iNEXT_3D/. \cr
In the detailed Online \code{iNEXT.3D} User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr
}
\author{
Anne Chao, K.H. Hu \cr

Maintainer: K.H. Hu <river19951224@gmail.com>
}
\references{
Chao, A., Chiu, C.-H., Villeger, S., Sun, I.-F., Thorn, S., Lin, Y.-C., Chiang, J. M. and Sherwin, W. B. (2019). An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs, 89, e01343. 10.1002/ecm.1343. \cr\cr
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
T.C. Hsieh, K. H. Ma, and Chao, A. (2016). iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456. \cr\cr
}
\keyword{ package }
