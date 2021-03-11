# iNEXT3D -------------------------------------------------------------------
#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT3D}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param data a matrix, data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
#' @param class a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which 
#' integrates several threshold functional diversity to get diversity.
#' @param q a numerical vector of the order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param rowsum a logical variable to check if the input object is raw data (species by sites matrix, \code{rowsum=FALSE}) or iNEXT default input (abundance counts or incidence frequencies, \code{rowsum=TRUE}).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots} .
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If NULL, then \code{endpoint} \code{=} double reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXT3D()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXT3D()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{class = 'PD'}.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc.
#' It is necessary when \code{class = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param reftime is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}. It will be use when \code{class = 'PD'}.
#' @param PDtype desired phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). It will be use when \code{class = 'PD'}. Default is \code{"PD"}.
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{class = 'FD' or 'AUC'}.
#' @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{class = 'FD'}.
#' @importFrom reshape2 dcast
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @importFrom stats rmultinom
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats optimize
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
#' 
#' @examples
#' ## example for abundance based data (list of vector)
#' # class = 'TD'
#' data(spider)
#' out1 <- iNEXT3D(spider, class = 'TD', q = c(0,1,2), datatype = "abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- iNEXT3D(data, class = 'PD', tree = tree, datatype = "abundance", q = c(0, 1, 2), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- iNEXT3D(data[,1], class = 'FD', distM = dij, datatype = "abundance", nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- iNEXT3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' t <- round(seq(10, 500, length.out = 20))
#' out5 <- iNEXT3D(ant$h500m, class = 'TD', q = 1, datatype = "incidence_freq", size = t)
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- iNEXT3D(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- iNEXT3D(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- iNEXT3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 0)
#' out8
#' 
#' @export
#' 
iNEXT3D <- function(data, class, q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, conf = 0.95, nboot = 50, 
                  tree = NULL, nT = NULL, reftime=NULL, PDtype = 'PD', distM, threshold = NULL) {
  if ( !(class %in% c('TD', 'PD', 'FD', 'AUC')) ) 
    stop("Please select one of below class: 'TD', 'PD', 'FD', 'AUC'", call. = FALSE)
  
  if (class == 'TD') {
    out = iNEXTTD(data, q = q, datatype = datatype, size = size, endpoint = endpoint, knots = knots, conf = conf, nboot = nboot)
  } else if (class == 'PD') {
    out = iNEXTPD(data, q = q, datatype = datatype, size = size, endpoint = endpoint, knots = knots, conf = conf, nboot = nboot, tree = tree, reftime = reftime, type = PDtype, nT = nT)
  } else if (class == 'FD') {
    out = iNEXTFD(data, q = q, datatype = datatype, size = size, endpoint = endpoint, knots = knots, conf = conf, nboot = nboot, distM = distM, threshold = threshold)
  } else if (class == 'AUC') {
    out = iNEXTAUC(data, q = q, datatype = datatype, size = size, endpoint = endpoint, knots = knots, conf = conf, nboot = nboot, distM = distM)
  }
  
  return(out)
}


# estimate3D -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage 
#' 
#' \code{estimate3D}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param data a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
#' by species incidence frequencies in each column or list.
#' @param class a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which 
#' integrates several threshold functional diversity to get diversity.
#' @param q a numerical vector of the order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param nboot the number of bootstrap times to obtain confidence interval. If confidence interval is not desired, use 0 to skip this time-consuming step.
#' @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1). 
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes. 
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes. 
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{class = 'PD'}.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc.
#' It is necessary when \code{class = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param reftime is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}. It will be use when \code{class = 'PD'}.
#' @param PDtype desired phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). It will be use when \code{class = 'PD'}. Default is \code{"PD"}.
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{class = 'FD' or 'AUC'}.
#' @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{class = 'FD'}.
#' @return a \code{data.frame} of species diversity table including the sample size, sample coverage,
#' method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#' @examples
#' # class = 'TD'
#' data(spider)
#' out1 <- estimate3D(spider, class = 'TD', q = c(0,1,2), datatype = "abundance", base = "size")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- estimate3D(data, class = 'PD', tree = tree, datatype = "abundance", base = "coverage")
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- estimate3D(data, class = 'FD', distM = dij, datatype = "abundance", base = "size")
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- estimate3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", nboot = 0, base = "coverage")
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- estimate3D(ant, class = 'TD', q = c(0,1,2), datatype = "incidence_freq", base="coverage", level=0.985)
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- estimate3D(data, class = 'PD', nT = nT, tree = tree, datatype = "incidence_raw", base = "size")
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- estimate3D(data, class = 'FD', distM = dij, datatype = "incidence_freq", base = "coverage")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- estimate3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20, base = "size")
#' out8
#' 
#' @export
estimate3D <- function (data, class, q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot=50,
                       conf = 0.95, tree, nT, reftime = NULL, PDtype = 'PD', distM, threshold = NULL) 
{
  if ( !(class %in% c('TD', 'PD', 'FD', 'AUC')) ) 
    stop("Please select one of below class: 'TD', 'PD', 'FD', 'AUC'", call. = FALSE)
  
  if (class == 'TD') {
    out = estimateTD(data, q = q, datatype = datatype, base = base, nboot = nboot, conf = conf, level = level)
  } else if (class == 'PD') {
    out = estimatePD(data, q = q, datatype = datatype, base = base, nboot = nboot, conf = conf, level = level, tree = tree, reftime = reftime, type = PDtype, nT = nT)
  } else if (class == 'FD') {
    out = estimateFD(data, q = q, datatype = datatype, base = base, nboot = nboot, conf = conf, level = level, distM = distM, threshold = threshold)
  } else if (class == 'AUC') {
    out = estimateAUC(data, q = q, datatype = datatype, base = base, nboot = nboot, conf = conf, level = level, distM = distM, tau = NULL)
  }
  
  return(out)
}


# Asy3D -------------------------------------------------------------------
#' Asymptotic diversity q profile 
#' 
#' \code{Asy3D} The estimated and empirical diversity of order q 
#' 
#' @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
#' @param class a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which 
#' integrates several threshold functional diversity to get diversity.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is seq(0, 2, by = 0.2).
#' @param datatype  data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' or species-by-site incidence frequencies data (\code{datatype = "incidence_freq"}). Default is "abundance".
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{class = 'PD'}.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc.
#' It is necessary when \code{class = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param reftime is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}. It will be use when \code{class = 'PD'}.
#' @param PDtype desired phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). It will be use when \code{class = 'PD'}. Default is \code{"PD"}.
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{class = 'FD' or 'AUC'}.
#' @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{class = 'FD'}.
#' @return a table of diversity q profile by 'Estimated'.
#' @examples
#' ## example for abundance-based data
#' # class = 'TD'
#' data(spider)
#' out1 <- Asy3D(spider, class = 'TD', datatype = "abundance")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- Asy3D(data, class = 'PD', datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- Asy3D(data, class = 'FD', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- Asy3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- Asy3D(ant, class = 'TD', datatype = "incidence_freq")
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- Asy3D(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = seq(0, 2, by = 0.25))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- Asy3D(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- Asy3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20)
#' out8
#' 
#' @export
Asy3D <- function(data, class = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, 
                 tree, nT, reftime = NULL, PDtype = 'PD', distM, threshold = NULL) {
  if ( !(class %in% c('TD', 'PD', 'FD', 'AUC')) ) 
    stop("Please select one of below class: 'TD', 'PD', 'FD', 'AUC'", call. = FALSE)
  
  if (class == 'TD') {
    out = AsyTD(data, q = q, datatype = datatype, nboot = nboot, conf = conf)
  } else if (class == 'PD') {
    out = AsyPD(data, q = q, datatype = datatype, nboot = nboot, conf = conf, tree = tree, reftime = reftime, type = PDtype, nT = nT)
  } else if (class == 'FD') {
    out = AsyFD(data, q = q, datatype = datatype, nboot = nboot, conf = conf, distM = distM, threshold = threshold)
  } else if (class == 'AUC') {
    out = AsyAUC(data, q = q, datatype = datatype, nboot = nboot, conf = conf, distM = distM, tau = NULL)
  }
  
  return(out)
}


# Obs3D -------------------------------------------------------------------
#' Empirical diversity q profile 
#' 
#' \code{Obs3D} The estimated and empirical diversity of order q 
#' 
#' @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
#' @param class a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which 
#' integrates several threshold functional diversity to get diversity.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is seq(0, 2, by = 0.2).
#' @param datatype  data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' or species-by-site incidence frequencies data (\code{datatype = "incidence_freq"}). Default is "abundance".
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{class = 'PD'}.
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc.
#' It is necessary when \code{class = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param reftime is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#' the pooled assemblage. Default is \code{NULL}. It will be use when \code{class = 'PD'}.
#' @param PDtype desired phylogenetic diversity type: \code{PDtype = "PD"} for Chao et al. (2010) phylogenetic diversity
#' and \code{PDtype = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). It will be use when \code{class = 'PD'}. Default is \code{"PD"}.
#' @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{class = 'FD' or 'AUC'}.
#' @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{class = 'FD'}.
#' 
#' @return a table of diversity q profile by 'Empirical'
#' 
#' @examples
#' ## example for abundance-based data
#' # class = 'TD'
#' data(spider)
#' out1 <- Obs3D(spider, class = 'TD', datatype = "abundance")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- Obs3D(data, class = 'PD', datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- Obs3D(data, class = 'FD', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- Obs3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- Obs3D(ant, class = 'TD', datatype = "incidence_freq")
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- Obs3D(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = seq(0, 2, by = 0.25))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- Obs3D(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- Obs3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20)
#' out8
#' 
#' @export
Obs3D <- function(data, class = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95,
                 tree, nT, reftime = NULL, PDtype = 'PD', distM, threshold = NULL) {
  if ( !(class %in% c('TD', 'PD', 'FD', 'AUC')) ) 
    stop("Please select one of below class: 'TD', 'PD', 'FD', 'AUC'", call. = FALSE)
  
  if (class == 'TD') {
    out = ObsTD(data, q = q, datatype = datatype, nboot = nboot, conf = conf)
  } else if (class == 'PD') {
    out = ObsPD(data, q = q, datatype = datatype, nboot = nboot, conf = conf, tree = tree, reftime = reftime, type = PDtype, nT = nT)
  } else if (class == 'FD') {
    out = ObsFD(data, q = q, datatype = datatype, nboot = nboot, conf = conf, distM = distM, threshold = threshold)
  } else if (class == 'AUC') {
    out = ObsAUC(data, q = q, datatype = datatype, nboot = nboot, conf = conf, distM = distM, tau = NULL)
  }
  
  return(out)
}


# ggiNEXT3D -------------------------------------------------------------------
#' ggplot2 extension for an iNEXT object
#' 
#' \code{ggiNEXT3D}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT3D}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation \cr (\code{facet.var="None"}); 
#'  a separate plot for each diversity order (\code{facet.var="Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"}); 
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"}); 
#'  use different colors for diversity orders (\code{color.var="Order.q"}); 
#'  use different colors for sites (\code{color.var="Assemblage"}); 
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).  
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' 
#' @examples
#' ## example for abundance based data (list of vector)
#' # class = 'TD'
#' data(spider)
#' out1 <- iNEXT3D(spider, class = 'TD', q = c(0,1,2), datatype = "abundance")
#' ggiNEXT3D(out1, facet.var = "Assemblage")
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- iNEXT3D(data, class = 'PD', tree = tree, datatype = "abundance", q = 0, nboot = 30)
#' ggiNEXT3D(out2, type = c(1, 3))
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- iNEXT3D(data[,1], class = 'FD', distM = dij, datatype = "abundance", nboot = 0)
#' ggiNEXT3D(out3)
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- iNEXT3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", nboot = 0)
#' ggiNEXT3D(out4)
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' t <- round(seq(10, 500, length.out = 20))
#' out5 <- iNEXT3D(ant$h500m, class = 'TD', q = 1, datatype = "incidence_freq", size = t)
#' ggiNEXT3D(out5)
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- iNEXT3D(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2))
#' ggiNEXT3D(out6, facet.var = "Order.q", color.var = "Assemblage")
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- iNEXT3D(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' ggiNEXT3D(out7)
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- iNEXT3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 0)
#' ggiNEXT3D(out8)
#' 
#' @export
#' 
#' 
#' 
ggiNEXT3D = function(outcome, type = 1:3, facet.var = "Assemblage", color.var = "Order.q"){
  if (sum(names(outcome) %in% c('DataInfo', 'iNextEst', 'AsyEst')) == 3) {
    class = 'TD'
    plottable = outcome$iNextEst
  } else if (sum(names(outcome) %in% c('PDInfo', 'PDiNextEst', 'PDAsyEst')) == 3) {
    class = 'PD'
    plottable = outcome$PDiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qPD', 'qD.LCL' = 'qPD.LCL', 'qD.UCL' = 'qPD.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qPD', 'qD.LCL' = 'qPD.LCL', 'qD.UCL' = 'qPD.UCL'))
    
  } else if (sum(names(outcome) %in% c('FDInfo', 'FDiNextEst', 'FDAsyEst')) == 3) {
    class = 'FD'
    plottable = outcome$FDiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qFD', 'qD.LCL' = 'qFD.LCL', 'qD.UCL' = 'qFD.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qFD', 'qD.LCL' = 'qFD.LCL', 'qD.UCL' = 'qFD.UCL'))
    
  } else if (sum(names(outcome) %in% c('AUCiNextEst', 'AUCAsyEst')) == 2) {
    class = 'AUC'
    plottable = outcome$AUCiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qAUC', 'qD.LCL' = 'qAUC.LCL', 'qD.UCL' = 'qAUC.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qAUC', 'qD.LCL' = 'qAUC.LCL', 'qD.UCL' = 'qAUC.UCL'))
    
  } else {stop("Please use the outcome from specified function 'iNEXT3D'")}
  
  SPLIT <- c("None", "Order.q", "Assemblage", "Both")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  TYPE <-  c(1, 2, 3)
  if(sum(!(type %in% TYPE)) >= 1)
    stop("invalid plot type")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  
  if(facet.var == "Order.q") color.var <- "Assemblage"
  if(facet.var == "Assemblage") color.var <- "Order.q"
  
  if ('m' %in% colnames(plottable$size_based) & 'm' %in% colnames(plottable$coverage_based)) datatype = 'abundance'
  if ('nt' %in% colnames(plottable$size_based) & 'nt' %in% colnames(plottable$coverage_based)) datatype = 'incidence'
  
  
  out = lapply(type, function(i) type_plot(x_list = plottable, i, class, datatype, facet.var, color.var))

  return(out)
}


# type_plot -------------------------------------------------------------------
type_plot = function(x_list, type, class, datatype, facet.var, color.var) {
  x_name <- colnames(x_list$size_based)[2]
  xlab_name <- ifelse(datatype == "incidence", "sampling units", "individuals")
  
  if (class == 'TD') {
    ylab_name = "Species Diversity"
  } else if (class == 'FD') {
    ylab_name = "Functional Diversity"
  } else if (class == 'AUC') {
    ylab_name = "Functional Diversity (AUC)"
  } else if (class == 'PD' & unique(x_list$size_based$Type) == 'PD') {
    ylab_name = "Phylogenetic Diversity"
  } else if (class == 'PD' & unique(x_list$size_based$Type) == 'meanPD') {
    ylab_name = "Phylogenetic Hill Diversity"
  } 
  
  
  if (type == 1) {
    output <- x_list$size_based
    output$y.lwr <- output$qD.LCL
    output$y.upr <- output$qD.UCL
    id <- match(c(x_name, "Method", "qD", "qD.LCL", "qD.UCL", "Assemblage", "Order.q"), names(output), nomatch = 0)
    output[,1:7] <- output[, id]
    
    xlab_name <- paste0("Number of ", xlab_name)
    
  } else if (type == 2) {
    output <- x_list$size_based
    if (length(unique(output$Order.q)) > 1) output <- subset(output, Order.q == unique(output$Order.q)[1])
    output$y.lwr <- output$SC.LCL
    output$y.upr <- output$SC.UCL
    id <- match(c(x_name, "Method", "SC", "SC.LCL", "SC.UCL", "Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(output), nomatch = 0)
    output[,1:10] <- output[, id]
    
    xlab_name <- paste0("Number of ", xlab_name)
    ylab_name <- "Sample Coverage"
    
  } else if (type == 3) {
    output <- x_list$coverage_based
    output$y.lwr <- output$qD.LCL
    output$y.upr <- output$qD.UCL
    id <- match(c("goalSC", "Method", "qD", "qD.LCL", "qD.UCL", "Assemblage", "Order.q", "SC", x_name), names(output), nomatch = 0)
    output[,1:9] <- output[, id]
    
    xlab_name <- "Sample Coverage"
    
  }
  
  if (facet.var == "None" & color.var == "None" & length(unique(output$Order.q)) > 1 & length(unique(output$Assemblage)) > 1) {
    color.var <- "Order.q"
    facet.var <- "Assemblage"
    warning ("invalid color.var and facet.var setting, the iNEXT3D object consists multiple orders and assemblage, change setting as Order.q and Assemblage")
  } else if (facet.var == "None" & color.var == "None" & length(unique(output$Order.q)) > 1) {
    color.var <- "Order.q"
    warning ("invalid color.var setting, the iNEXT3D object consists multiple orders, change setting as Order.q")
  } else if (facet.var == "None" & color.var == "None" & length(unique(output$Assemblage)) > 1) { 
    color.var <- "Assemblage" 
    warning ("invalid color.var setting, the iNEXT3D object consists multiple assemblage, change setting as Assemblage")
    }
    
  
  
  title <- c("Sample-size-based sampling curve", "Sample completeness curve", "Coverage-based sampling curve")[type]
  colnames(output)[1:7] <- c("x", "Method", "y", "LCL", "UCL", "Assemblage", "Order.q")
  
  if (class == 'PD') {
    output$Reftime <- round(output$Reftime, 3)
    output$Reftime <- factor(paste0("Ref.time = ", output$Reftime), levels = paste0("Ref.time = ", unique(output$Reftime)))
  }
  if (class == 'FD') {
    output$threshold <- round(output$threshold, 3)
    output$threshold <- factor(paste0("tau = ", output$threshold), levels = paste0("tau = ", unique(output$threshold)))
  }
  
  if (color.var == "None") {
    if (levels(factor(output$Order.q)) > 1 & length(unique(output$Assemblage)) > 1) {
      warning ("invalid color.var setting, the iNEXT3D object consists multiple assemblages and orders, change setting as Both")
      color.var <- "Both"
      output$col <- output$shape <- paste(output$Assemblage, output$Order.q, sep="-")
      
    } else if (length(unique(output$Assemblage)) > 1) {
      warning ("invalid color.var setting, the iNEXT3D object consists multiple assemblages, change setting as Assemblage")
      color.var <- "Assemblage"
      output$col <- output$shape <- output$Assemblage
    } else if (levels(factor(output$Order.q)) > 1){
      warning ("invalid color.var setting, the iNEXT3D object consists multiple orders, change setting as Order.q")
      color.var <- "Order.q"
      output$col <- output$shape <- factor(output$Order.q)
    } else {
      output$col <- output$shape <- rep(1, nrow(output))
    }
  } else if (color.var == "Order.q") {     
    output$col <- output$shape <- factor(output$Order.q)
  } else if (color.var == "Assemblage") {
    if (length(unique(output$Assemblage)) == 1) {
      warning ("invalid color.var setting, the iNEXT3D object do not consist multiple assemblages, change setting as Order.q")
      output$col <- output$shape <- factor(output$Order.q)
    }
    output$col <- output$shape <- output$Assemblage
  } else if (color.var == "Both") {
    if (length(unique(output$Assemblage)) == 1) {
      warning ("invalid color.var setting, the iNEXT3D object do not consist multiple assemblages, change setting as Order.q")
      output$col <- output$shape <- factor(output$Order.q)
    }
    output$col <- output$shape <- paste(output$Assemblage, output$Order.q, sep="-")
  }
  
  if (type == 2) output$col = output$shape = output$Assemblage
  
  data.sub = output
  output$Method[output$Method == "Observed"] = "Rarefaction"
  output$lty <- factor(output$Method, levels = c("Rarefaction", "Extrapolation"))
  output$col <- factor(output$col)
  data.sub <- data.sub[which(data.sub$Method == "Observed"),]
  
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7",  "#0072B2", "#D55E00"))
  
  g <- ggplot(output, aes_string(x = "x", y = "y", colour = "col")) + 
    geom_line(aes_string(linetype = "lty"), lwd=1.5) +
    geom_point(aes_string(shape = "shape"), size=5, data = data.sub) +
    geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr", fill = "factor(col)", colour = "NULL"), alpha = 0.2) +
    scale_fill_manual(values = cbPalette) +
    scale_colour_manual(values = cbPalette) +
    guides(linetype = guide_legend(title = "Method"),
           colour = guide_legend(title = "Guides"), 
           fill = guide_legend(title = "Guides"), 
           shape = guide_legend(title = "Guides"))
  
  g = g + theme_bw() + 
    labs(x = xlab_name, y = ylab_name) + 
    ggtitle(title) + 
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
  
  
  if (facet.var == "Order.q") {
    if(length(levels(factor(output$Order.q))) == 1 & type != 2){
      warning("invalid facet.var setting, the iNEXT3D object do not consist multiple orders.")      
    } else {
      odr_grp <- labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      
      if (class == 'PD') {
        g <- g + facet_wrap(Reftime ~ Order.q, nrow = 1, labeller = odr_grp)
      } else if (class == 'FD') {
        g <- g + facet_wrap(threshold ~ Order.q, nrow = 1, labeller = odr_grp)
      } else {g <- g + facet_wrap( ~ Order.q, nrow = 1, labeller = odr_grp)}
      
      if (color.var == "Both") {
        g <- g + guides(colour = guide_legend(title = "Guides", ncol = length(levels(factor(output$Order.q))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
      if(type == 2){
        g <- g + theme(strip.background = element_blank(), strip.text.x = element_blank())
        
      }
    }
  }
  
  if(facet.var == "Assemblage"){
    if(length(unique(output$Assemblage)) == 1) {
      warning("invalid facet.var setting, the iNEXT3D object do not consist multiple assemblages")
    }else{
      if (class == 'PD') {
        g <- g + facet_wrap(Reftime ~ Assemblage, nrow = 1)
      } else if (class == 'FD') {
        g <- g + facet_wrap(threshold ~ Assemblage, nrow = 1)
      } else {g <- g + facet_wrap( ~ Assemblage, nrow = 1)}
      
      if(color.var == "Both"){
        g <- g + guides(colour = guide_legend(title = "Guides", nrow = length(levels(factor(output$Order.q)))),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  
  if(facet.var == "Both"){
    if(length(levels(factor(output$Order.q))) == 1 | length(unique(output$Assemblage)) == 1){
      warning("invalid facet.var setting, the iNEXT3D object do not consist multiple assemblages or orders.")
    }else{
      odr_grp <- labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      
      if (class == 'PD') {
        g <- g + facet_wrap(Assemblage + Reftime ~ Order.q, labeller = odr_grp)
      } else if (class == 'FD') {
        g <- g + facet_wrap(Assemblage + threshold ~ Order.q, labeller = odr_grp)
      } else {g <- g + facet_wrap(Assemblage ~ Order.q, labeller = odr_grp)}
      
      if(color.var == "both"){
        g <- g +  guides(colour = guide_legend(title = "Guides", nrow = length(levels(factor(output$Assemblage))), byrow = TRUE),
                         fill = guide_legend(title = "Guides"))
      }
    }
  }
  
  return(g)
}


# ggAsy3D -------------------------------------------------------------------
#' ggplot for Asymptotic diversity
#'
#' \code{ggAsy3D} Plots q-profile, time-profile, and tau-profile based on the outcome of \code{Asy3D} using the ggplot2 package.\cr
#' It will only show the confidence interval of 'Estimated'.
#'
#' @param outcome the outcome of the functions \code{Asy3D} .\cr
#' @return a figure of asymptotic of empirical three-divrsity\cr\cr
#'
#' @examples
#' ## example for abundance-based data
#' # class = 'TD'
#' data(spider)
#' out1 <- rbind(Asy3D(spider, class = 'TD', datatype = "abundance"), Obs3D(spider, class = 'TD', datatype = "abundance"))
#' ggAsy3D(out1)
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- Asy3D(data, class = 'PD', datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25), nboot = 30, PDtype = "meanPD")
#' ggAsy3D(out2, profile = "q")
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- Asy3D(data, class = 'FD', distM = dij, datatype = "abundance", q = c(0, 1, 2), threshold = seq(0, 0.6, 0.1), nboot = 0)
#' ggAsy3D(out3, profile = "tau")
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- Asy3D(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' ggAsy3D(out4)
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- rbind(Asy3D(ant, class = 'TD', datatype = "incidence_freq"), Obs3D(ant, class = 'TD', datatype = "incidence_freq"))
#' ggAsy3D(out5)
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- Asy3D(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2), reftime = seq(0.1, 82.8575, length.out = 40))
#' ggAsy3D(out6, profile = "time")
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- Asy3D(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' ggAsy3D(out7)
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- Asy3D(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20)
#' ggAsy3D(out8)
#'
#' @export
ggAsy3D <- function(outcome, profile = 'q'){
  if (sum(unique(outcome$Method) %in% c("Asymptotic", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'Asy3D'")
  
  if (!(profile %in% c('q', 'time', 'tau')))
    stop("Please select one of 'q', 'time', 'tau' profile.")
  
  if (sum(colnames(outcome)[1:6] == c('Order.q', 'qD', 'qD.LCL', 'qD.UCL', 'Assemblage', 'Method')) == 6) {
    class = 'TD'
  } else if (sum(colnames(outcome)[1:6] == c('Order.q', 'qPD', 'qPD.LCL', 'qPD.UCL', 'Assemblage', 'Method')) == 6) {
    class = 'PD'
  } else if (sum(colnames(outcome)[1:6] == c('Order.q', 'qFD', 'qFD.LCL', 'qFD.UCL', 'Assemblage', 'Method')) == 6) {
    class = 'FD'
  } else if (sum(colnames(outcome)[1:6] == c('Order.q', 'qAUC', 'qAUC.LCL', 'qAUC.UCL', 'Assemblage', 'Method')) == 6) {
    class = 'AUC'
  } else {stop("Please use the outcome from specified function 'Asy3D'")}
  
  ## TD & q-profile ##
  if (class == 'TD') {
    out = ggplot(outcome, aes(x = Order.q, y = qD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method == 'Asymptotic')) out = out + labs(x = 'Order q', y = 'Asymptotic Taxonomic Diversity')
      if (unique(outcome$Method == 'Empirical')) out = out + labs(x = 'Order q', y = 'Empirical Taxonomic Diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qD.LCL, ymax = qD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Taxonomic Diversity')
    }
  }
  
  ## PD & q-profile ##
  if (class == 'PD' & profile == 'q') {
    outcome$Reftime = paste('Reftime = ', round(outcome$Reftime, 3), sep = '')
    out = ggplot(outcome, aes(x = Order.q, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Asymptotic Phylogenetic Diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Empirical Phylogenetic Diversity')
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Asymptotic Phylogenetic Hill Diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Empirical Phylogenetic Hill Diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Phylogenetic Diversity')
      if (unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Phylogenetic Hill Diversity')
    }
    out = out + facet_grid(.~Reftime, scales = "free_y")
  }
  
  ## PD & time-profile ##
  if (class == 'PD' & profile == 'time') {
    outcome$Order.q = paste('q = ', outcome$Order.q, sep = '')
    out = ggplot(outcome, aes(x = Reftime, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Asymptotic Phylogenetic Diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Empirical Phylogenetic Diversity')
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Asymptotic Phylogenetic Hill Diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Empirical Phylogenetic Hill Diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Phylogenetic Diversity')
      if (unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Phylogenetic Hill Diversity')
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## FD & q-profile ##
  if (class == 'FD' & profile == 'q') {
    outcome$tau = paste('Tau = ', round(outcome$tau, 3), sep = '')
    out = ggplot(outcome, aes(x = Order.q, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic Functional Diversity')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'Order q', y = 'Empirical Functional Diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Functional Diversity')
    }
    out = out + facet_grid(.~tau, scales = "free_y")
  }
  
  ## FD & tau-profile ##
  if (class == 'FD' & profile == 'tau') {
    outcome$Order.q = paste('Order q = ', outcome$Order.q, sep = '')
    out = ggplot(outcome, aes(x = tau, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'tau', y = 'Asymptotic Functional Diversity')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'tau', y = 'Empirical Functional Diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'tau', y = 'Functional Diversity')
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## AUC & q-profile ##
  if (class == 'AUC') {
    out = ggplot(outcome, aes(x = Order.q, y = qAUC, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qAUC.LCL, ymax = qAUC.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic Functional Diversity (AUC)')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'Order q', y = 'Empirical Functional Diversity (AUC)')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qAUC.LCL, ymax = qAUC.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Functional Diversity (AUC)')
    }
  }
  
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  
  out = out +
    scale_colour_manual(values = cbPalette) + theme_bw() + 
    scale_fill_manual(values = cbPalette) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 16),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    guides(linetype = guide_legend(keywidth = 2.5))
  
  return(out)
}


#' @useDynLib iNEXT3D, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



