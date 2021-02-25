# iNEXT -------------------------------------------------------------------
#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param data a matrix, data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
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
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXT()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXT()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
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
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
#' @examples
#' ## example for abundance based data (list of vector)
#' # class = 'TD'
#' data(spider)
#' out1 <- iNEXT(spider, class = 'TD', q = c(0,1,2), datatype = "abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- iNEXT(data, class = 'PD', tree = tree, datatype = "abundance", q = c(0, 1, 2), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- iNEXT(data[,1], class = 'FD', distM = dij, datatype = "abundance", nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- iNEXT(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' t <- round(seq(10, 500, length.out = 20))
#' out5 <- iNEXT(ant$h500m, class = 'TD', q = 1, datatype = "incidence_freq", size = t)
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- iNEXT(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = c(0, 1, 2))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- iNEXT(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- iNEXT(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 0)
#' out8
#' 
#' @export
#' 
iNEXT <- function(data, class, q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, conf = 0.95, nboot = 50, 
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


# estimateD -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage 
#' 
#' \code{estimateD}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param data a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
#' by species incidence frequencies in each column or list.
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
#' out1 <- estimateD(spider, class = 'TD', q = c(0,1,2), datatype = "abundance", base = "size")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- estimateD(data, class = 'PD', tree = tree, datatype = "abundance", base = "coverage")
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- estimateD(data, class = 'FD', distM = dij, datatype = "abundance", base = "size")
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- estimateD(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", nboot = 0, base = "coverage")
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- estimateD(ant, class = 'TD', q = c(0,1,2), datatype = "incidence_freq", base="coverage", level=0.985)
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- estimateD(data, class = 'PD', nT = nT, tree = tree, datatype = "incidence_raw", base = "size")
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- estimateD(data, class = 'FD', distM = dij, datatype = "incidence_freq", base = "coverage")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- estimateD(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20, base = "size")
#' out8
#' 
#' @export
estimateD <- function (data, class, q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot=50,
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


# AsyD -------------------------------------------------------------------
#' Asymptotic diversity q profile 
#' 
#' \code{AsyD} The estimated and empirical diversity of order q 
#' 
#' @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
#' @param class a nomial selection which kind of diversity you want to calculate. You should choose one from below: 'TD', 'PD', 'FD', 'AUC'. Default is 'TD'.
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
#' out1 <- AsyD(spider, class = 'TD', datatype = "abundance")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- AsyD(data, class = 'PD', datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- AsyD(data, class = 'FD', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- AsyD(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- AsyD(ant, class = 'TD', datatype = "incidence_freq")
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- AsyD(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = seq(0, 2, by = 0.25))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- AsyD(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- AsyD(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20)
#' out8
#' 
#' @export
AsyD <- function(data, class = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, 
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


# ObsD -------------------------------------------------------------------
#' Empirical diversity q profile 
#' 
#' \code{ObsD} The estimated and empirical diversity of order q 
#' 
#' @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
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
#' out1 <- ObsD(spider, class = 'TD', datatype = "abundance")
#' out1
#' 
#' # class = 'PD'
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out2 <- ObsD(data, class = 'PD', datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25), nboot = 30)
#' out2
#' 
#' # class = 'FD'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out3 <- ObsD(data, class = 'FD', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out3
#' 
#' # class = 'AUC'
#' data(FunDdata.abu)
#' data <- FunDdata.abu$data
#' dij <-  FunDdata.abu$dij
#' out4 <- ObsD(data = data[,2], class = 'AUC', distM = dij, datatype = "abundance", q = seq(0, 2, 0.5), nboot = 0)
#' out4
#' 
#' ## example for incidence-based data
#' # class = 'TD'
#' data(ant)
#' out5 <- ObsD(ant, class = 'TD', datatype = "incidence_freq")
#' out5
#' 
#' # class = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- ObsD(data, class = 'PD', nT = nT, datatype = "incidence_raw", tree = tree, q = seq(0, 2, by = 0.25))
#' out6
#' 
#' # class = 'FD'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out7 <- ObsD(data, class = 'FD', distM = dij, datatype = "incidence_freq")
#' out7
#' 
#' # class = 'AUC'
#' data(FunDdata.inc)
#' data <- FunDdata.inc$data
#' dij <-  FunDdata.inc$dij
#' out8 <- ObsD(data, class = 'AUC', distM = dij, datatype = "incidence_freq", nboot = 20)
#' out8
#' 
#' @export
ObsD <- function(data, class = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95,
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


# ggiNEXT -------------------------------------------------------------------
#' ggplot2 extension for an iNEXT object
#' 
#' \code{ggiNEXT}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
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
#' @param grey a logical variable to display grey and white ggplot2 theme. 
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, class = 'TD', q = 0, datatype = "abundance")
#' ggiNEXT(x = out1, type = 1)
#' ggiNEXT(x = out1, type = 2)
#' ggiNEXT(x = out1, type = 3)
#' 
#' # single-assemblage incidence data with three orders q
#' data(ant)
#' size <- round(seq(10, 500, length.out=20))
#' y <- iNEXT(ant$h500m, class = 'TD', q = c(0,1,2), datatype = "incidence_freq", size = size)
#' ggiNEXT(y, se=FALSE, color.var="Order.q")
#' 
#' # multiple-assemblage abundance data with three orders q
#' data(spider)
#' z <- iNEXT(spider, class = 'TD', q = c(0,1,2), datatype = "abundance")
#' ggiNEXT(z, facet.var="Assemblage", color.var="Order.q")
#' ggiNEXT(z, facet.var="Both", color.var="Both")
#' @export
#' 
ggiNEXT <- function(x, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE){  
  UseMethod("ggiNEXT", x)
}


# ggAsyD -------------------------------------------------------------------
#' ggplot for Asymptotic diversity
#'
#' \code{ggAsyD} Plots q-profile based on the outcome of \code{AsyD} using the ggplot2 package.\cr
#' It will only show the confidence interval of 'Estimated'.
#'
#' @param outcome the outcome of the functions \code{AsyD} .\cr
#' @return a figure of estimated sample completeness with order q\cr\cr
#'
#' @examples
#' ## Type (1) example for abundance-based data
#' ## Ex.1
#' data(spider)
#' out1 <- AsyD(spider, class = 'TD', datatype = "abundance")
#' ggAsyD(out1)
#' 
#' ## Type (2) example for incidence-based data
#'
#' ## Ex.2
#' data(ant)
#' out2 <- AsyD(ant, class = 'TD', datatype = "incidence_freq", nboot = 0)
#' ggAsyD(out2)
#'
#' @export
ggAsyD <- function(outcome){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  if (sum(unique(outcome$Method) %in% c("Asymptotic", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'AsyD'")
  ggplot(outcome, aes(x=Order.q, y=qD, colour=Assemblage, lty=Method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = outcome[outcome$Method=="Asymptotic",],
                aes(ymin=qD.LCL, ymax=qD.UCL, fill=Assemblage), alpha=0.2, linetype=0) +
    # geom_ribbon(data = outcome[outcome$method=="Empirical",],
    #             aes(ymin=qD.LCL, ymax=qD.UCL, fill=Assemblage), alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    scale_linetype_manual(values = c("Asymptotic"=1, "Empirical"=2)) +
    labs(x="Order q", y="Species diversity") +
    # theme_bw(base_size = 18) +
    theme(text=element_text(size=18)) +
    theme(legend.position="bottom", legend.box = "vertical",
          legend.key.width = unit(1.2,"cm"),
          # plot.margin = unit(c(1.5,0.3,1.2,0.3), "lines"),
          legend.title=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-5,-10))
}


#' @useDynLib iNEXT3D, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


