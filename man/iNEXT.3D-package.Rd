\name{iNEXT.3D-package}
\alias{iNEXT.3D-package}
\alias{iNEXT.3D}
\docType{package}
\title{
Interpolation and Extrapolation for three dimensions of biodiversity
}
\description{
\code{iNEXT.3D} (iNterpolation and EXTrapolation for three dimensions of biodiversity) is a sequel to iNEXT (Hsieh et al., 2016). 
Here the three dimensions (3D) of diversity includes taxonomic diversity (TD), phylogenetic diversity (PD) and functional diversity (FD). An online version "iNEXT.3D Online"
(https://chao.shinyapps.io/iNEXT_3D/) is also available for users without an R background.

A unified framework based on Hill numbers (for TD) and their generalizations (Hill-Chao numbers, for PD and FD) 
is adopted to quantify 3D. In this framework, TD quantifies the effective number of species, PD quantifies the effective total branch length, mean-PD (PD divided by 
tree depth) quantifies the effective number of lineages, and FD quantifies the effective number of virtual functional 
groups (or functional "species"). Thus, TD, mean-PD, and FD are all in the same units of species/lineage equivalents 
and can be meaningfully compared; see Chao et al. (2014) for the basic standardization theory for TD, and Chao et al. (2021) for a review of the unified theory for 3D.

For each of the three dimensions of biodiversity, \code{iNEXT.3D} features two statistical analyses (non-asymptotic and asymptotic):

(1) A non-asymptotic approach based on interpolation and extrapolation for 3D diversity (i.e., Hill-Chao numbers)

\code{iNEXT.3D} computes the estimated 3D diversity for standardized samples with a common sample size or sample completeness. This approach aims to compare diversity estimates for equally-large (with a common sample size) or equally-complete (with a common sample coverage) samples; it is based on the seamless rarefaction and extrapolation (R/E) sampling curves of Hill-Chao numbers for q = 0, 1 and 2. For each dimension of biodiversity, \code{iNEXT.3D} offers three types of R/E sampling curves:

* Sample-size-based (or size-based) R/E sampling curves: This type of sampling curve plots the diversity estimates with respect to sample size.

* Coverage-based R/E sampling curves: This type of sampling curve plots the diversity estimates with respect to sample coverage.

* Sample completeness curve: This curve depicts how sample coverage varies with sample size. The sample completeness curve provides a bridge between the size- and coverage-based R/E sampling curves.

(2) An asymptotic approach to infer asymptotic 3D diversity (i.e., Hill-Chao numbers)

\code{iNEXT.3D} computes the estimated asymptotic 3D diversity and also plots 3D diversity profiles for q between 0 and 2, in comparison with the observed diversity. Typically, the asymptotic estimates for q \eqn{\geq} 1 are reliable, but for q < 1 (especially for q = 0, species richness), the asymptotic estimates represent only lower bounds.   \cr

There are six main functions in this package: \cr

1. \code{iNEXT3D} Computes 3D diversity estimates of order q = 0, 1 and 2 for rarefied and extrapolated samples at specified sample coverage values and sample sizes. \cr

2. \code{ggiNEXT3D} Visualizes the output from the function \code{iNEXT3D}. \cr

3. \code{ObsAsy3D} Computes observed and asymptotic diversity of order q for 3D diversity. \cr

4. \code{ggObsAsy3D} Visualizes the output from the function \code{ObsAsy3D}. \cr

5. \code{estimate3D} Computes 3D diversity of order q = 0, 1 and 2 with a particular user-specified level of sample size or sample coverage. \cr

6. \code{DataInfo3D} Provides basic data information. \cr

}
\author{
Anne Chao, Kai-Hsiang Hu \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
*	Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
*	T.C. Hsieh, K. H. Ma, and Chao, A. (2016). iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods in Ecology and Evolution, 7, 1451-1456. \cr\cr
}
\keyword{package}
