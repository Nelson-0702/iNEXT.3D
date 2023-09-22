\name{iNEXT.3D-package}
\alias{iNEXT.3D-package}
\alias{iNEXT.3D}
\docType{package}
\title{
Interpolation and Extrapolation for three dimensions of biodiversity
}
\description{
\code{iNEXT.3D} (Interpolation and Extrapolation for three dimensions of biodiversity) is a sequel to iNEXT (Hsieh et al., 2016). 
Here the three dimensions (3D) of diversity includes taxonomic diversity (TD), phylogenetic diversity (PD) and functional diversity (FD). An online version "iNEXT.3D Online"
(https://chao.shinyapps.io/iNEXT_3D/) is also available for users without an R background.

A unified framework based on Hill numbers (for TD) and their generalizations (Hill-Chao numbers, for PD and FD) 
is adopted to quantify 3D. In this framework, TD quantifies the effective number of species, mean-PD (PD divided by 
tree depth) quantifies the effective number of lineages, and FD quantifies the effective number of virtual functional 
groups (or functional "species"). Thus, TD, mean-PD, and FD are all in the same units of species/lineage equivalents
and can be meaningfully compared; see Chao et al. (2021) for a review of the unified framework. 

\code{iNEXT.3D} also features two types of R/E sampling curves:

* Sample-size-based R/E sampling curves:\code{iNEXT.3D} computes diversity estimates for rarefied and 
extrapolated samples up to an appropriate size. This type of sampling curve plots the diversity estimates with respect to sample size.

* Coverage-based R/E sampling curves: \code{iNEXT.3D} computes diversity estimates for rarefied and 
extrapolated samples with sample completeness (as measured by sample coverage) up to an appropriate coverage. 
This type of sampling curve plots the diversity estimates with respect to sample coverage. 

\code{iNEXT.3D} also plots the above two types of sampling curves and a sample completeness curve by `ggiNEXT3D`. The sample completeness curve provides a bridge between these two types of curves.

This package contains six main functions: \cr\cr

1. \code{iNEXT3D} computing interpolation and extrapolation with order q for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage or sample size. \cr

2. \code{ggiNEXT3D} visualizing the output from the function \code{iNEXT.3D}. \cr

3. \code{AO3D} computing empirical and asymptotic diversity with order q for three dimensions: taxonomic, phylogenetic and functional diversity. \cr

4. \code{ggAO3D} visualizing the output from the function \code{AO3D}. \cr

5. \code{estimate3D} computes taxonomic, phylogenetic, and functional diversity (Hill-Chao family with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage. \cr

6. \code{DataInfo3D} provides basic data information. \cr

}
\author{
Anne Chao, K.H. Hu \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
*	Chao, A., Chiu, C.-H., Villeger, S., Sun, I.-F., Thorn, S., Lin, Y.-C., Chiang, J. M. and Sherwin, W. B. (2019). An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs, 89, e01343. 10.1002/ecm.1343. \cr\cr
*	Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
*	T.C. Hsieh, K. H. Ma, and Chao, A. (2016). iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456. \cr\cr
}
\keyword{ package }
