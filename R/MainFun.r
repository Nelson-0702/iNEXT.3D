#' Data information for reference samples
#' 
#' \code{DataInfo3D} provides basic data information for 3D diversity based on a reference sample.
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is a single matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc.
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of \code{PDtree}).  
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' 
#' @return a data.frame including basic data information.\cr\cr 
#' For abundance data, basic information shared by TD, mean-PD and FD includes assemblage name (\code{Assemblage}),
#' sample size (\code{n}), observed species richness (\code{S.obs}), sample coverage estimates of the reference sample (\code{SC(n)}), 
#' sample coverage estimate for twice the reference sample size (\code{SC(2n)}). Other additional information is given below.\cr\cr
#' (1) TD: the first five species abundance counts (\code{f1}--\code{f5}).\cr\cr
#' (2) Mean-PD: the observed total branch length in the phylogenetic tree (\code{PD.obs}), 
#' the number of singletons (\code{f1*}) and doubletons (\code{f2*}) in the node/branch abundance set, as well as the total branch length 
#' of those singletons (\code{g1}) and of those doubletons (\code{g2}), and the reference time (\code{Reftime}).\cr\cr
#' (3) FD (\code{FDtype = "AUC"}): the minimum distance among all non-diagonal elements in the distance matrix (\code{dmin}), the mean distance
#' (\code{dmean}), and the maximum distance (\code{dmax}) in the distance matrix.\cr \cr
#' (4) FD (\code{FDtype = "tau_value"}): the number of singletons (\code{a1*}) and of doubletons (\code{a2*}) among the functionally indistinct
#'  set at the specified threshold level \code{'Tau'}, as well as the total contribution of singletons (\code{h1}) and of doubletons (\code{h2})
#'   at the specified threshold level \code{'Tau'}.\cr\cr
#'  
#'  For incidence data, the basic information for TD includes assemblage name (\code{Assemblage}), number of sampling units (\code{T}), 
#'  total number of incidences (\code{U}), observed species richness (\code{S.obs}), 
#'  sample coverage estimates of the reference sample (\code{SC(T)}), sample coverage estimate for twice the reference sample size
#'  (\code{SC(2T)}), as well as the first five species incidence frequency counts (\code{Q1}--\code{Q5}). For mean-PD and FD, output is similar to that
#'  for abundance data.
#'  
#' 
#' @examples
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' DataInfo3D(Brazil_rainforest_abun_data, diversity = 'TD', datatype = "abundance")
#' 
#' \donttest{
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' DataInfo3D(data, diversity = 'PD', datatype = "abundance", PDtree = tree)
#' }
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' DataInfo3D(data, diversity = 'FD', datatype = "abundance", 
#'            FDdistM = distM, FDtype = 'AUC')
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' DataInfo3D(Fish_incidence_data, diversity = 'TD', datatype = "incidence_raw")
#' 
#' \donttest{
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' DataInfo3D(data, diversity = 'PD', datatype = "incidence_raw", PDtree = tree)
#' }
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' DataInfo3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'            FDdistM = distM, FDtype = 'AUC')
#' 
#'
#' @export
DataInfo3D <- function(data, diversity = 'TD', datatype = "abundance", nT = NULL, PDtree, PDreftime = NULL, FDdistM, FDtype = "AUC", FDtau = NULL){
  
  if (diversity == 'TD') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT, to.datalist = TRUE)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    out <- TDinfo(data, datatype)
    
  } 
  
  
  if (diversity == 'PD') {
    
    # if (datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    if (datatype == "incidence_freq") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.")
    
    checkdatatype = check.datatype(data, datatype, nT = nT, raw.to.inci = F)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    nT = checkdatatype[[3]]
    
    checktree = check.tree(data, datatype, PDtree, PDreftime, nT)
    PDreftime = checktree[[1]]
    mytree = checktree[[2]]
    mydata = checktree[[3]]
    
    
    if(datatype=='abundance'){
      
      out <- lapply(mydata, function(x){
        datainf(data = x, datatype, phylotr = mytree,reft = PDreftime) %>% mutate(Reftime = PDreftime)
      }) %>% do.call(rbind,.) %>% 
        mutate(Assemblage = rep(names(mydata), each = length(PDreftime))) %>%
        select(Assemblage, n, S.obs, `SC(n)`, `SC(2n)`, PD.obs, `f1*`, `f2*`, g1, g2, Reftime)
      
    }else if (datatype=='incidence_raw'){
      
      out <- lapply(mydata, function(x){
        datainf(data = x, datatype, phylotr = mytree,reft = PDreftime) %>% mutate(Reftime = PDreftime)
      }) %>% do.call(rbind,.) %>% 
        mutate(Assemblage = rep(names(mydata), each = length(PDreftime))) %>%
        select(Assemblage,`T`, U, S.obs, `SC(T)`, `SC(2T)`, PD.obs, `Q1*`, `Q2*`, R1, R2, Reftime)
      
    }
    
  } 
  
  
  if (diversity == 'FD' & FDtype == 'tau_values') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, FDtau)
    FDtau = checkdistM[[1]]
    distM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    out = lapply(FDtau, function(tau) {
      
      out <- lapply(1:length(dat), function(i) {
        
        x = dat[[i]]
        aivi = data_transform(x, distM, tau, datatype, integer = TRUE)
        
        if (datatype == "abundance") {
          
          n = sum(x)
          f1 = sum(x == 1)
          f2 = sum(x == 2)
          f0.hat <- ifelse(f2 == 0, (n-1) / n * f1 * (f1-1) / 2, (n-1) / n * f1^2 / 2 / f2) 
          A <- ifelse(f1 > 0, n * f0.hat / (n * f0.hat + f1), 1)
          Chat <- 1 - f1/n * A
          Chat2n <- Coverage(x, "abundance", 2*sum(x))
          
          multiple = tibble('Assemblage' = names(dat)[i], 
                            'n' = n, 
                            'S.obs' = sum(x > 0), 
                            'SC(n)' = Chat, 'SC(2n)' = Chat2n, 
                            'a1*' = sum(aivi$ai == 1), 'a2*' = sum(aivi$ai == 2), 
                            'h1' = sum(aivi$vi[aivi$ai == 1,]), 'h2' = sum(aivi$vi[aivi$ai == 2,]),
                            'Tau' = tau)
          
        } 
        else if (datatype == "incidence_freq") {

          nT = x[1]
          x = x[-1]
          U <- sum(x)
          Q1 = sum(x == 1)
          Q2 = sum(x == 2)
          Q0.hat <- ifelse(Q2 == 0, (nT-1) / nT * Q1 * (Q1-1) / 2, (nT-1) / nT * Q1^2 / 2 / Q2)
          A <- ifelse(Q1 > 0, nT * Q0.hat / (nT * Q0.hat + Q1), 1)
          Chat <- 1 - Q1/U * A
          Chat2T <- Coverage(c(nT,x), "incidence_freq", 2*nT)

          multiple = tibble('Assemblage' = names(dat)[i],
                            'T' = nT,
                            'U' = U,
                            'S.obs' = sum(x > 0),
                            'SC(T)' = Chat, 'SC(2T)' = Chat2T,
                            'a1*' = sum(aivi$ai == 1), 'a2*' = sum(aivi$ai == 2),
                            'h1' = sum(aivi$vi[aivi$ai == 1,]), 'h2' = sum(aivi$vi[aivi$ai == 2,]),
                            'Tau' = tau)
        }
        
        return(multiple)
        
      }) %>% do.call(rbind,.)
      
      
    }) %>% do.call(rbind,.)
    
    rownames(out) = NULL
  } 
  
  
  if (diversity == 'FD' & FDtype == 'AUC') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, threshold = FALSE)
    distM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    Tau = t(sapply(dat, function(i){
      
      if(datatype=='abundance') {
        
        tmp <- matrix(i/sum(i),ncol =1)
        
      }
      else if(datatype=='incidence_freq'){

        tmp <- matrix(i[-1]/sum(i[-1]), ncol = 1)

      }
      
      dmean <- sum ( (tmp %*% t(tmp) ) * distM)
      distM <- distM[tmp > 0, tmp > 0]
      dmin <- min(distM[lower.tri(distM)])
      dmax <- max(distM[distM > 0])
      
      c(dmin, dmean, dmax)
    }))
    
    if (datatype == "abundance") {
      
      out <- cbind(TDinfo(dat, datatype)[,1:5], Tau)
      colnames(out)[6:8] = c("dmin", "dmean", "dmax")
      rownames(out) = NULL
      
    } else {
      
      out <- cbind(TDinfo(dat, datatype)[,1:6], Tau)
      colnames(out)[7:9] = c("dmin", "dmean", "dmax")
      rownames(out) = NULL
      
    }
    
    return(out)
  }
  
  return(out)
}


#' @useDynLib iNEXT.3D, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



#' iNterpolation and EXTrapolation of Hill-Chao number
#' 
#' \code{iNEXT3D}: Interpolation and extrapolation of Hill-Chao number with order q
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}.
#' @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
#' If NULL, then \code{endpoint} \code{=} double reference sample size.
#' @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
#' each knot represents a particular sample size for which diversity estimate will be calculated.  
#' If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXT3D()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
#' If the \code{endpoint} is larger than the reference sample size, then \code{iNEXT3D()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is a single matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param PDtree (required argument only for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).  
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of \code{PDtree}).  
#' @param PDtype (argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' @param FDcut_number (argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is \code{FDcut_number = 50}. A larger value can be set to obtain more accurate AUC value.
#' 
#' @importFrom reshape2 dcast
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import tidytree
#' @import tibble
#' @importFrom stats rmultinom
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats optimize
#' 
#' @return a list of three objects: \cr
#' \code{$TDInfo} (or \code{$PDInfo}, \code{$FDInfo}) for summarizing data information. Details please refer to the output of \code{DataInfo3D}. \cr\cr
#' \code{$TDiNextEst} (or \code{$PDiNextEst}, \code{$FDiNextEst}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics by \code{"$size_based"} and \code{"$coverage_based"}. \cr\cr
#'    In \code{"$size_based"}, the output includes:
#'    \item{Assemblage}{the name of assemblage.} 
#'    \item{Order.q}{the diversity order of q.}
#'    \item{m, mT}{the target number of sample size (or sampling units).}
#'    \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target sample size is less than, equal to, or greater than the size of the reference sample.}
#'    \item{qTD, qPD, qFD}{the estimated diversity estimate.}
#'    \item{qTD.LCL, qPD.LCL, qFD.LCL and qTD.UCL, qPD.UCL, qFD.UCL}{the bootstrap lower and upper confidence limits for the diversity.}
#'    \item{SC}{the standardized coverage value.}
#'    \item{SC.LCL, SC.UCL}{the bootstrap lower and upper confidence limits for coverage.}
#'    \item{Reftime}{reference times for PD.}
#'    \item{Type}{"PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages) for PD.}
#'    \item{Tau}{the threshold of functional distinctiveness between any two species for FD (under FDtype = tau_values).}
#'  Similar output is obtained for \code{"$coverage_based"}. \cr
#' \code{$TDAsyEst} (or \code{$PDAsyEst}, \code{$FDAsyEst}) for showing asymptotic diversity estimates along with related statistics: \cr\cr
#'    \item{Assemblage}{the name of assemblage.} 
#'    \item{qTD, qPD, qFD}{the diversity order of q.}
#'    \item{TD_obs, PD_obs, FD_obs}{the observed diversity.}
#'    \item{TD_asy, PD_asy, FD_asy}{the asymptotic diversity estimate.}
#'    \item{s.e.}{standard error of asymptotic diversity.}
#'    \item{qTD.LCL, qPD.LCL, qFD.LCL and qTD.UCL, qPD.UCL, qFD.UCL}{the bootstrap lower and upper confidence limits for asymptotic diversity.}
#'    \item{Reftime}{reference times for PD.}
#'    \item{Type}{"PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages) for PD.}
#'    \item{Tau}{the threshold of functional distinctiveness between any two species for FD (under FDtype = tau_values).}
#' 
#' 
#' @examples
#' \donttest{
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' output_TD_abun <- iNEXT3D(Brazil_rainforest_abun_data, diversity = 'TD', q = c(0,1,2), 
#'                                 datatype = "abundance")
#' output_TD_abun
#' 
#' 
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' output_PD_abun <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "abundance", 
#'                                 nboot = 20, PDtree = tree)
#' output_PD_abun
#' 
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' output_FD_abun <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 0, 
#'                                 FDdistM = distM, FDtype = 'AUC')
#' output_FD_abun
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' output_TD_inci <- iNEXT3D(Fish_incidence_data, diversity = 'TD', q = c(0, 1, 2), 
#'                                 datatype = "incidence_raw")
#' output_TD_inci
#' 
#' 
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' output_PD_inci <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), 
#'                                 datatype = "incidence_raw", nboot = 20, PDtree = tree)
#' output_PD_inci
#' 
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' output_FD_inci <- iNEXT3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, 
#'                                 FDdistM = distM, FDtype = 'AUC')
#' output_FD_inci
#' }
#' 
#' 
#' @export
iNEXT3D <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, nboot = 50, conf = 0.95, nT = NULL, 
                    PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50) {
  
  if ( !(diversity %in% c('TD', 'PD', 'FD')) ) 
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)
  
  if (diversity == 'TD') {
    
    data.original = data
    datatype.original = datatype
    
    checkdatatype = check.datatype(data, datatype, nT = nT, to.datalist = TRUE)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    size = check.size(data, datatype, size, endpoint, knots)
    
    Fun <- function(x, q, size, assem_name){
      
      x <- as.numeric(unlist(x))
      unconditional_var <- TRUE
      
      if(datatype == "abundance"){
        if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
        out <- iNEXT.Ind(Spec=x, q=q, m=size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots=knots, nboot=nboot, conf=conf,unconditional_var)
      }
      
      if(datatype == "incidence_freq"){
        t <- x[1]
        y <- x[-1]

        if(t>sum(y)){
          warning("Insufficient data to provide reliable estimators and associated s.e.")
        }

        if(sum(x)==0) stop("Zero incidence frequencies in one or more sample sites")

        out <- iNEXT.Sam(Spec=x, q=q, t=size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots=knots, nboot=nboot, conf=conf)
      }
      
      if(unconditional_var){
        
        out <- lapply(out, function(out_) cbind(Assemblage = assem_name, out_))
        
      }else{
        
        out[[1]] <- cbind(Assemblage = assem_name, out[[1]])
      }
      
      out
    }
    
    z <- qnorm(1-(1-conf)/2)
    
    if(is.null(names(data))){
      names(data) <- sapply(1:length(data), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:length(data), function(i) {
      tmp <- Fun(data[[i]],q,size[[i]],names(data)[i])
      tmp
    })
    
    out <- list(size_based = do.call(rbind,lapply(out, function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out, function(out_){out_[[2]]})))
    
    index <- rbind(asyTD(data, datatype, c(0, 1, 2), nboot, conf),
                   obsTD(data, datatype, c(0, 1, 2), nboot, conf))
    index = index[order(index$Assemblage),]
    LCL <- index$qTD.LCL[index$Method=='Asymptotic']
    UCL <- index$qTD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qTD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Observed & index$Order.q==0] <- index$Observed[index$LCL<index$Observed & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "qTD", "TD_obs", "TD_asy", "s.e.", "qTD.LCL", "qTD.UCL")
    
    
    out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
    out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
    # out$size_based <- as_tibble(out$size_based)
    # out$coverage_based <- as_tibble(out$coverage_based)
    
    info <- DataInfo3D(data.original, diversity = 'TD', datatype.original, nT)
    
    out <- list("TDInfo"=info, "TDiNextEst"=out, "TDAsyEst"=index)
  } 
  
  if (diversity == 'PD') {
    
    # if (datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    if (datatype == "incidence_freq") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.")
    
    checkdatatype = check.datatype(data, datatype, nT = nT, raw.to.inci = F)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    nT = checkdatatype[[3]]
    
    checktree = check.tree(data, datatype, PDtree, PDreftime, nT)
    PDreftime = checktree[[1]]
    mytree = checktree[[2]]
    mydata = checktree[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    size = check.size(mydata, datatype, size, endpoint, knots)
    PDtype = check.PDtype(PDtype)
    
    
    out <- inextPD(datalist = mydata, datatype = datatype, phylotr = mytree, q = q, reft = PDreftime, m=size,
                   cal = PDtype, nboot = nboot, conf = conf, unconditional_var = TRUE)
    out$size_based = out$size_based %>% select(-c('s.e.', 'SC.s.e.'))
    out$coverage_based = out$coverage_based %>% select(-('s.e.'))
    
    ## AsyEst table ##
    index <- rbind(asymPD(datalist = mydata, datatype = datatype, phylotr = mytree,q = c(0, 1, 2), 
                          reft = PDreftime, cal = PDtype, nboot, conf),
                   EmpPD(datalist = mydata, datatype = datatype, phylotr = mytree,q = c(0, 1, 2), 
                         reft = PDreftime, cal = PDtype, nboot, conf))
    index = index[order(index$Assemblage),]
    LCL <- index$qPD.LCL[index$Method=='Asymptotic']
    UCL <- index$qPD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qPD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Observed & index$Order.q==0] <- index$Observed[index$LCL<index$Observed & index$Order.q==0]
    index$Order.q <- c('q = 0 PD','q = 1 PD','q = 2 PD')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "qPD", "PD_obs", "PD_asy", "s.e.", "qPD.LCL", "qPD.UCL")
    index$Reftime = PDreftime
    index$Type = PDtype
    
    info <- DataInfo3D(data, diversity = 'PD', datatype = datatype, nT, PDtree = PDtree, PDreftime = PDreftime)
    
    
    out = list("PDInfo"=info, "PDiNextEst"=out, "PDAsyEst"=index)
    
  } 
  
  if (diversity == 'FD' & FDtype == 'tau_values') {
    
    data.original = data
    datatype.original = datatype
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, FDtau)
    FDtau = checkdistM[[1]]
    dist = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    size = check.size(dat, datatype, size, endpoint, knots)
    
    
    FUN <- function(e){
      if(inherits(dat, "list")){
        ## size-based
        temp1 = iNextFD(datalist = dat, dij = dist, q = q, datatype = datatype, tau = FDtau,
                        nboot = nboot, conf = conf, m = size)
        temp1$qFD.LCL[temp1$qFD.LCL<0] <- 0;temp1$SC.LCL[temp1$SC.LCL<0] <- 0
        temp1$SC.UCL[temp1$SC.UCL>1] <- 1
        if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nT'
        
        ## coverage-based
        temp2 <- lapply(1:length(dat), function(i) invChatFD(datalist = dat[i], dij = dist, q = q, datatype = datatype,
                                                             level = unique(Coverage(data = dat[[i]], datatype = datatype, m = size[[i]])), 
                                                             nboot = nboot, conf = conf, tau = FDtau)) %>% do.call(rbind,.)
        temp2$qFD.LCL[temp2$qFD.LCL<0] <- 0
        
        if (datatype == 'incidence_freq') colnames(temp2)[colnames(temp2) == 'm'] = 'nT'
        # temp1$Type = "FD"
        # temp2$Type = "FD"
        ans <- list(size_based = temp1, coverage_based = temp2)
        return(ans)
      }else{
        return(NULL)
      }
    }
    out <- tryCatch(FUN(e), error = function(e){return()})
    out$size_based = out$size_based %>% select(-c('s.e.', 'SC.s.e.'))
    out$coverage_based = out$coverage_based %>% select(-('s.e.'))
    
    ## AsyEst table ##
    index <- rbind(FDtable_est(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype, 
                               nboot = nboot, conf = conf, tau = FDtau),
                   FDtable_mle(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype, 
                               nboot = nboot, conf = conf, tau = FDtau))
    index <- index %>% arrange(., Assemblage)
    LCL <- index$qFD.LCL[index$Method=='Asymptotic']
    UCL <- index$qFD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qFD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Observed & index$Order.q==0] <- index$Observed[index$LCL<index$Observed & index$Order.q==0]
    index$Order.q <- c('q = 0 FD(single tau)','q = 1 FD(single tau)','q = 2 FD(single tau)')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "qFD", "FD_obs", "FD_asy", "s.e.", "qFD.LCL", "qFD.UCL")
    index$Tau = FDtau
    
    info <- DataInfo3D(data.original, diversity = 'FD', datatype = datatype.original, FDdistM = FDdistM, FDtype = 'tau_values', FDtau = FDtau, nT = nT)
    out = list("FDInfo" = info, "FDiNextEst" = out, "FDAsyEst" = index)
    
  } 
  
  if (diversity == 'FD' & FDtype == 'AUC') {
   
    data.original = data
    datatype.original = datatype
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, threshold = FALSE)
    dist = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    size = check.size(dat, datatype, size, endpoint, knots)
    FDcut_number = check.FDcut_number(FDcut_number)
    
    
    FUN <- function(e){
      if(inherits(dat, "list")){
        ## size-based
        temp1 = AUCtable_iNextFD(datalist = dat, dij = dist, q = q, datatype = datatype,
                                 tau = NULL, nboot = nboot, conf = conf, m = size, FDcut_number = FDcut_number)
        temp1$qFD.LCL[temp1$qFD.LCL<0] <- 0; temp1$SC.LCL[temp1$SC.LCL<0] <- 0
        temp1$SC.UCL[temp1$SC.UCL>1] <- 1
        if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nT'
        
        ## coverage-based
        temp2 <- lapply(1:length(dat), function(i) AUCtable_invFD(datalist = dat[i], dij = dist, q = q, datatype = datatype,
                                                                  level = unique(Coverage(data = dat[[i]], datatype = datatype, m = size[[i]])), 
                                                                  nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number)) %>% do.call(rbind,.)
        temp2$qFD.LCL[temp2$qFD.LCL<0] <- 0
        if (datatype == 'incidence_freq') colnames(temp2)[colnames(temp2) == 'm'] = 'nT'
        
        ans <- list(size_based = temp1, coverage_based = temp2)
        return(ans)
      }else{
        return(NULL)
      }
    }
    out <- tryCatch(FUN(e), error = function(e){return()})
    out$size_based = out$size_based %>% data.frame %>% select(-c('s.e.', 'SC.s.e.'))
    out$coverage_based = out$coverage_based %>% data.frame %>% select(-('s.e.'))
    
    ## AsyEst table ##
    index <- rbind(AUCtable_est(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype,
                                nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number),
                   AUCtable_mle(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype,
                                nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number))
    index = index[order(index$Assemblage),]
    LCL <- index$qFD.LCL[index$Method=='Asymptotic']
    UCL <- index$qFD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qFD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Observed & index$Order.q==0] <- index$Observed[index$LCL<index$Observed & index$Order.q==0]
    index$Order.q <- c('q = 0 FD(AUC)','q = 1 FD(AUC)','q = 2 FD(AUC)')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "qFD", "FD_obs", "FD_asy", "s.e.", "qFD.LCL", "qFD.UCL")
    
    info <- DataInfo3D(data.original, diversity = 'FD', datatype = datatype.original, FDdistM = FDdistM, FDtype = 'AUC', nT = nT)
    out = list("FDInfo" = info, "FDiNextEst" = out, "FDAsyEst" = index)
    
  }
  
  if (datatype != 'abundance'){
    out[[2]]$size_based <- rename(out[[2]]$size_based, c("mT" = "nT"))
    out[[2]]$coverage_based <- rename(out[[2]]$coverage_based, c("mT" = "nT"))
  }
  
  class(out) <- c("iNEXT3D")
  
  return(out)
}


#' ggplot2 extension for an iNEXT3D object
#' 
#' \code{ggiNEXT3D} is a \code{ggplot} extension for \code{iNEXT3D} object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve.
#' @param output an \code{iNEXT3D} object computed by \code{iNEXT3D}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation (\code{facet.var = "None"}); 
#'  a separate plot for each diversity order (\code{facet.var = "Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var = "Assemblage"}); 
#'  a separate plot for each combination of diversity order and assemblage (\code{facet.var = "Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var = "None"}); 
#'  use different colors for diversity orders (\code{color.var = "Order.q"}); 
#'  use different colors for assemblages/sites (\code{color.var = "Assemblage"}); 
#'  use different colors for combinations of diversity order and assemblage (\code{color.var = "Both"}).  
#' @return a \code{ggplot2} object for sample-size-based rarefaction/extrapolation curve (\code{type = 1}), sample completeness curve (\code{type = 2}), and coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' 
#' 
#' @examples
#' \donttest{
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' output_TD_abun <- iNEXT3D(Brazil_rainforest_abun_data, diversity = 'TD', q = c(0,1,2), 
#'                                 datatype = "abundance")
#' ggiNEXT3D(output_TD_abun, facet.var = "Assemblage")
#' 
#' 
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' output_PD_abun <- iNEXT3D(data, diversity = 'PD', q = c(0,1,2), datatype = "abundance", 
#'                                 nboot = 20, PDtree = tree)
#' ggiNEXT3D(output_PD_abun, type = c(1, 3))
#' 
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' output_FD_abun <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 0, 
#'                                 FDdistM = distM, FDtype = 'AUC')
#' ggiNEXT3D(output_FD_abun)
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' output_TD_inci <- iNEXT3D(Fish_incidence_data, diversity = 'TD', q = 1, 
#'                                 datatype = "incidence_raw")
#' ggiNEXT3D(output_TD_inci)
#' 
#' 
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' output_PD_inci <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", 
#'                                 nboot = 20, PDtree = tree)
#' ggiNEXT3D(output_PD_inci, facet.var = "Order.q", color.var = "Assemblage")
#' 
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' output_FD_inci <- iNEXT3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, 
#'                                 FDdistM = distM, FDtype = 'AUC')
#' ggiNEXT3D(output_FD_inci)
#' }
#' 
#' 
#' @export
ggiNEXT3D = function(output, type = 1:3, facet.var = "Assemblage", color.var = "Order.q"){
  
  if (sum(names(output) %in% c('TDInfo', 'TDiNextEst', 'TDAsyEst')) == 3) {
    
    class = 'TD'
    plottable = output$TDiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qTD', 'qD.LCL' = 'qTD.LCL', 'qD.UCL' = 'qTD.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qTD', 'qD.LCL' = 'qTD.LCL', 'qD.UCL' = 'qTD.UCL'))
    
  } else if (sum(names(output) %in% c('PDInfo', 'PDiNextEst', 'PDAsyEst')) == 3) {
    
    class = 'PD'
    plottable = output$PDiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qPD', 'qD.LCL' = 'qPD.LCL', 'qD.UCL' = 'qPD.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qPD', 'qD.LCL' = 'qPD.LCL', 'qD.UCL' = 'qPD.UCL'))
    
  } else if (sum(names(output) %in% c('FDInfo', 'FDiNextEst', 'FDAsyEst')) == 3) {
    
    if ("Tau" %in% colnames(output$FDiNextEst$size_based)) class = 'FD' else class = 'FD(AUC)'
    
    plottable = output$FDiNextEst
    plottable$size_based = rename(plottable$size_based, c('qD' = 'qFD', 'qD.LCL' = 'qFD.LCL', 'qD.UCL' = 'qFD.UCL'))
    plottable$coverage_based = rename(plottable$coverage_based, c('qD' = 'qFD', 'qD.LCL' = 'qFD.LCL', 'qD.UCL' = 'qFD.UCL'))
    
  } else {stop("Please use the output from the function 'iNEXT3D'")}
  
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
  if ('mT' %in% colnames(plottable$size_based) & 'mT' %in% colnames(plottable$coverage_based)) datatype = 'incidence'
  
  
  out = lapply(type, function(i) type_plot(x_list = plottable, i, class, datatype, facet.var, color.var))
  if (length(type) == 1) out = out[[1]]
  
  return(out)
}


type_plot = function(x_list, type, class, datatype, facet.var, color.var) {
  
  x_name <- colnames(x_list$size_based)[3]
  xlab_name <- ifelse(datatype == "incidence", "sampling units", "individuals")
  
  if (class == 'TD') {
    
    ylab_name = "Taxonomic diversity"
    
  } else if (class == 'FD') {
    
    ylab_name = "Functional diversity"
    
  } else if (class == 'FD(AUC)') {
    
    ylab_name = "Functional diversity (AUC)"
    
  } else if (class == 'PD' & unique(x_list$size_based$Type) == 'PD') {
    
    ylab_name = "Phylogenetic diversity"
    
  } else if (class == 'PD' & unique(x_list$size_based$Type) == 'meanPD') {
    
    ylab_name = "Mean phylogenetic diversity"
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
    ylab_name <- "Sample coverage"
    
  } else if (type == 3) {
    
    output <- x_list$coverage_based %>% tibble
    output$y.lwr <- output$qD.LCL
    output$y.upr <- output$qD.UCL
    id <- match(c("SC", "Method", "qD", "qD.LCL", "qD.UCL", "Assemblage", "Order.q", x_name), names(output), nomatch = 0)
    output[,1:8] <- output[, id]
    
    xlab_name <- "Sample coverage"
    
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
    
    output$Tau <- round(output$Tau, 3)
    output$Tau <- factor(paste0("Tau = ", output$Tau), levels = paste0("Tau = ", unique(output$Tau)))
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
  tmp = output %>% filter(Method == "Observed") %>% mutate(Method = "Extrapolation")
  output$Method[output$Method == "Observed"] = "Rarefaction"
  output = rbind(output, tmp)
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
        g <- g + facet_grid(Tau ~ Order.q, labeller = odr_grp, scales = 'free_y')
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
        g <- g + facet_grid(Tau ~ Assemblage, scales = 'free_y')
      } else {g <- g + facet_wrap(. ~ Assemblage, nrow = 1)}
      
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
        # if(length(unique(output$Reftime)) == 1) outp <- outp + theme(strip.background = element_blank(), strip.text.x = element_blank())
      } else if (class == 'FD') {
        g <- g + facet_grid(Assemblage + Tau ~ Order.q, labeller = odr_grp, scales = 'free_y')
      } else {g <- g + facet_wrap(Assemblage ~ Order.q, labeller = odr_grp)}
      
      if(color.var == "both"){
        g <- g +  guides(colour = guide_legend(title = "Guides", nrow = length(levels(factor(output$Assemblage))), byrow = TRUE),
                         fill = guide_legend(title = "Guides"))
      }
    }
  }
  
  return(g)
}


#' Compute 3D diversity estimates with a particular set of sample sizes/coverages
#' 
#' \code{estimate3D} computes 3D diversity (Hill-Chao number with q = 0, 1 and 2) with a particular set of user-specified levels of sample sizes or sample coverages. If no sample sizes or coverages are specified, this function by default computes 3D diversity estimates for the minimum sample coverage or minimum sample size among all samples extrapolated to double reference sizes.
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param base selection of sample-size-based (\code{base = "size"}) or coverage-based (\code{base = "coverage"}) rarefaction and extrapolation.
#' @param level A numerical vector specifying the particular sample sizes or sample coverages (between 0 and 1) for which 3D diversity estimates (q =0, 1 and 2) will be computed. \cr
#' If \code{base = "coverage"} (default) and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes. \cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes. 
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is a single matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of \code{PDtree}).  
#' @param PDtype (argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' @param FDcut_number (argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is \code{FDcut_number = 50}. A larger value can be set to obtain more accurate AUC value.
#' 
#' @return a data.frame of diversity table including the following arguments: (when \code{base = "coverage"})
#' \item{Assemblage}{the assemblage name}
#' \item{Order.q}{the diversity order of q.}
#' \item{SC}{the target standardized coverage value.}
#' \item{m, mT}{the corresponding sample size (or sampling units) for the standardized coverage value.}
#' \item{qTD, qPD, qFD}{the estimated diversity of order q for the target coverage value. The estimate for complete coverage (or size = infinity) represents the estimated asymptotic diversity.}
#' \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample.}
#' \item{s.e.}{standard error of diversity estimate.}
#' \item{qTD.LCL, qPD.LCL, qFD.LCL and qTD.UCL, qPD.UCL, qFD.UCL}{the bootstrap lower and upper confidence limits for the diversity of order q at the specified level. (with a default value of 0.95)}
#' \item{Reftime}{reference times for PD.}
#' \item{Type}{"PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages) for PD.}
#' \item{Tau}{the threshold of functional distinctiveness between any two species for FD (under \code{FDtype = tau_values}).}
#' Similar output is obtained for \code{base = "size"}. \cr\cr
#' 
#' 
#' @examples
#' \donttest{
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' output_est_TD_abun <- estimate3D(Brazil_rainforest_abun_data, diversity = 'TD', q = c(0,1,2), 
#'                                  datatype = "abundance", base = "size")
#' output_est_TD_abun
#' 
#' 
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' output_est_PD_abun <- estimate3D(data, diversity = 'PD', datatype = "abundance", 
#'                                  level = c(1500, 3500), base = "size", PDtree = tree)
#' output_est_PD_abun
#' 
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' output_est_FD_abun <- estimate3D(data, diversity = 'FD', datatype = "abundance", 
#'                                  level = c(0.93, 0.97), base = "coverage", nboot = 10, 
#'                                  FDdistM = distM, FDtype = 'AUC')
#' output_est_FD_abun
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' output_est_TD_inci <- estimate3D(Fish_incidence_data, diversity = 'TD', q = c(0,1,2), 
#'                                  datatype = "incidence_raw", base = "coverage")
#' output_est_TD_inci
#' 
#' 
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' output_est_PD_inci <- estimate3D(data, diversity = 'PD', datatype = "incidence_raw", 
#'                                  level = c(0.975, 0.99), base = "coverage", PDtree = tree)
#' output_est_PD_inci
#' 
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' output_est_FD_inci <- estimate3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'                                  level = c(30, 70), base = "size", nboot = 10, 
#'                                  FDdistM = distM, FDtype = 'AUC')
#' output_est_FD_inci
#' }
#' 
#' 
#' @export
estimate3D <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot = 50, conf = 0.95, nT = NULL, 
                       PDtree, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50) {
  
  if ( !(diversity %in% c('TD', 'PD', 'FD')) ) 
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)
  
  if (diversity == 'TD') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT, to.datalist = TRUE)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    base = check.base(base)
    level = check.level(data, datatype, base, level)
    
    if (base == "size") {
      out <- invSize(data, q, datatype, size = level, nboot, conf = conf)
    } else if (base == "coverage") {
      out <- invChat(data, q, datatype, C = level, nboot, conf = conf)
    }
    out$qTD.LCL[out$qTD.LCL<0] <- 0
    
  } 
  
  if (diversity == 'PD') {
    
    # if (datatype == "incidence_freq") stop ("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    if (datatype == "incidence_freq") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.")
    
    checkdatatype = check.datatype(data, datatype, nT = nT, raw.to.inci = F)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    nT = checkdatatype[[3]]
    
    checktree = check.tree(data, datatype, PDtree, PDreftime, nT)
    PDreftime = checktree[[1]]
    mytree = checktree[[2]]
    mydata = checktree[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    base = check.base(base)
    PDtype = check.PDtype(PDtype)
    level = check.level(mydata, datatype, base, level)
    
    
    if (base == "size") {
      
      out <- inextPD(datalist = mydata, datatype = datatype, phylotr = mytree,q = q, 
                     reft = PDreftime,m = lapply(mydata, function(i) level), cal = PDtype, nboot=nboot, conf = conf, unconditional_var = FALSE)$size_based %>% 
        select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
      out = out %>% .[,c(1:4, 9, 5:8, 10:11)]
      
    } else if (base == "coverage") {
      
      out <- invChatPD(datalist = mydata, datatype = datatype, phylotr = mytree, q = q,
                       reft = PDreftime, cal = PDtype, level = level, nboot, conf)
      
    }
    out$qPD.LCL[out$qPD.LCL<0] <- 0
  } 
  
  if (diversity == 'FD' & FDtype == 'tau_values') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, FDtau)
    FDtau = checkdistM[[1]]
    FDdistM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    base = check.base(base)
    level = check.level(dat, datatype, base, level)
    
    
    if (base == "size") {
      
      out = iNextFD(datalist = dat,dij = FDdistM,q = q,datatype = datatype,tau = FDtau,
                    nboot = nboot,conf = conf,m = lapply(1:length(dat), function(i) level)) %>% 
        select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
      out$qFD.LCL[out$qFD.LCL<0] <- 0
      # out$SC.LCL[out$SC.LCL<0] <- 0
      # out$SC.UCL[out$SC.UCL>1] <- 1
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nT'
      out = out %>% .[,c(1:4, 9, 5:8, 10)]
      
    } else if (base == "coverage") {
      
      out <- invChatFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                       level = level, nboot = nboot, conf = conf, tau = FDtau)
      out$qFD.LCL[out$qFD.LCL<0] <- 0
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nT'
      
    }
    
  } 
  
  if (diversity == 'FD' & FDtype == 'AUC') {
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, threshold = FALSE)
    FDdistM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    base = check.base(base)
    level = check.level(dat, datatype, base, level)
    FDcut_number = check.FDcut_number(FDcut_number)
    
    
    if (base == 'size') {
      
      out = AUCtable_iNextFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                             tau = NULL, nboot = nboot, conf = conf, m = lapply(1:length(dat), function(i) level),
                             FDcut_number = FDcut_number) %>% 
        select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
      out$qFD.LCL[out$qFD.LCL<0] <- 0
      # out$SC.LCL[out$SC.LCL<0] <- 0
      # out$SC.UCL[out$SC.UCL>1] <- 1
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nT'
      out = out %>% .[,c(1:4, 9, 5:8)]
      
    } else if (base == 'coverage') {
      
      out <- AUCtable_invFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                            level = level, nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number)
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nT'
      
    }
    out$qFD.LCL[out$qFD.LCL<0] <- 0
    
    out = data.frame(out)
  }
  
  if (datatype != 'abundance'){
    out <- rename(out, c("mT" = "nT"))
  }
  
  return(out)
}


#' Asymptotic diversity and observed diversity of order q
#' 
#' \code{ObsAsy3D}: The estimated asymptotic and observed diversity of order q 
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list of matrices/data.frames (species by sampling units); data can also be input as a single matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (\code{nT}, see below) must be specified. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence/occurrence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is a single matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified (i.e., \code{names(nT) = NULL}), then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param method computing type. Select \code{'Asymptotic'} or \code{'Observed'}.
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of \code{PDtree}).  
#' @param PDtype (argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' @param FDcut_number (argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is \code{FDcut_number = 50}. A larger value can be set to obtain more accurate AUC value.
#' 
#' @return a table of diversity table including the following arguments: 
#' \item{Assemblage}{the assemblage name.}
#' \item{Order.q}{the diversity order of q.}
#' \item{qTD, qPD, qFD}{the estimated asymptotic diversity or observed diversity of order q.} 
#' \item{s.e.}{standard error of diversity.}
#' \item{qTD.LCL, qPD.LCL, qFD.LCL and qTD.UCL, qPD.UCL, qFD.UCL}{the bootstrap lower and upper confidence limits for the diversity.}
#' \item{Method}{"Asymptotic" means asymptotic diversity and "Observed" means observed diversity.}
#' \item{Reftime}{the reference times for PD.}
#' \item{Type}{"PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages) for PD.}
#' \item{Tau}{the threshold of functional distinctiveness between any two species for FD (under FDtype = tau_values).}
#' 
#' 
#' @examples
#' \donttest{
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' output_ObsAsy_TD_abun <- ObsAsy3D(Brazil_rainforest_abun_data, diversity = 'TD', 
#'                                   datatype = "abundance")
#' output_ObsAsy_TD_abun
#' 
#' 
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' output_ObsAsy_PD_abun <- ObsAsy3D(data, diversity = 'PD', q = seq(0, 2, by = 0.25), 
#'                                   datatype = "abundance", nboot = 20, PDtree = tree)
#' output_ObsAsy_PD_abun
#' 
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' output_ObsAsy_FD_abun <- ObsAsy3D(data, diversity = 'FD', q = seq(0, 2, 0.5), 
#'                                   datatype = "abundance", nboot = 10, FDdistM = distM, 
#'                                   FDtype = 'AUC')
#' output_ObsAsy_FD_abun
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' output_ObsAsy_TD_inci <- ObsAsy3D(Fish_incidence_data, diversity = 'TD', 
#'                                   datatype = "incidence_raw")
#' output_ObsAsy_TD_inci
#' 
#' 
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' output_ObsAsy_PD_inci <- ObsAsy3D(data, diversity = 'PD', q = seq(0, 2, by = 0.25), 
#'                                   datatype = "incidence_raw", nboot = 20, PDtree = tree)
#' output_ObsAsy_PD_inci
#' 
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' output_ObsAsy_FD_inci <- ObsAsy3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'                                   nboot = 20, FDdistM = distM, FDtype = 'AUC')
#' output_ObsAsy_FD_inci
#' }
#' 
#' 
#' @export
ObsAsy3D <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT = NULL, method = c('Asymptotic', 'Observed'),
                 PDtree, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50) {
  
  if ( !(diversity %in% c('TD', 'PD', 'FD')) ) 
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)
  
  if (diversity == "TD") {
    checkdatatype = check.datatype(data, datatype, nT = nT, to.datalist = TRUE)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    
    
    if (sum(method == "Asymptotic") == length(method)) 
      
      out = asyTD(data, datatype, q, nboot, conf) else if (sum(method == "Observed") == length(method)) 
        
        out = obsTD(data, datatype, q, nboot, conf) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
          
          out = rbind(asyTD(data, datatype, q, nboot, conf), 
                      obsTD(data, datatype, q, nboot, conf))
  }
  
  if (diversity == "PD") {
    
    # if (datatype == "incidence_freq")
    #   stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    if (datatype == "incidence_freq") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.")
    
    checkdatatype = check.datatype(data, datatype, nT = nT, raw.to.inci = F)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    nT = checkdatatype[[3]]
    
    checktree = check.tree(data, datatype, PDtree, PDreftime, nT)
    PDreftime = checktree[[1]]
    mytree = checktree[[2]]
    mydata = checktree[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    PDtype = check.PDtype(PDtype)
    
    
    if (sum(method == "Asymptotic") == length(method)) 
      
      out = asymPD(datalist = mydata, datatype = datatype, phylotr = mytree, 
                   q = q, reft = PDreftime, cal = PDtype, nboot, conf) else if (sum(method == "Observed") == length(method)) 
                     
                     out = EmpPD(datalist = mydata, datatype = datatype, phylotr = mytree, 
                                 q = q, reft = PDreftime, cal = PDtype, nboot, conf) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
                                   
                                   out = rbind(asymPD(datalist = mydata, datatype = datatype, phylotr = mytree, 
                                                      q = q, reft = PDreftime, cal = PDtype, nboot, conf), 
                                               EmpPD(datalist = mydata, datatype = datatype, phylotr = mytree, 
                                                     q = q, reft = PDreftime, cal = PDtype, nboot, conf))
    
  }
  
  if (diversity == "FD" & FDtype == "tau_values") {
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, FDtau)
    FDtau = checkdistM[[1]]
    distM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    
    
    if (sum(method == "Asymptotic") == length(method)) 
      out = FDtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                        nboot = nboot, conf = conf, tau = FDtau) else if (sum(method == "Observed") == length(method)) 
                          
                          out = FDtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                            nboot = nboot, conf = conf, tau = FDtau) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
                                              
                                              out = rbind(FDtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                      nboot = nboot, conf = conf, tau = FDtau), 
                                                          FDtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                      nboot = nboot, conf = conf, tau = FDtau))
    
  }
  
  if (diversity == "FD" & FDtype == "AUC") {
    
    checkdatatype = check.datatype(data, datatype, nT = nT)
    datatype = checkdatatype[[1]]
    data = checkdatatype[[2]]
    
    checkdistM = check.dist(data, datatype, FDdistM, threshold = FALSE)
    distM = checkdistM[[2]]
    dat = checkdistM[[3]]
    
    q = check.q(q)
    conf = check.conf(conf)
    nboot = check.nboot(nboot)
    FDcut_number = check.FDcut_number(FDcut_number)
    
    if (sum(method == "Asymptotic") == length(method)) 
      out = AUCtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                         nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number) else if (sum(method == "Observed") == length(method)) 
                           
                           out = AUCtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                              nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
                                                
                                                out = rbind(AUCtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                         nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number), 
                                                            AUCtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                         nboot = nboot, conf = conf, tau = NULL, FDcut_number = FDcut_number))
    
  }
  
  return(out)
}


#' ggplot2 extension for plotting q-profile, time-profile, and tau-profile
#'
#' \code{ggObsAsy3D} is a \code{ggplot} extension for \code{ObsAsy3D} object to plot q-profile, time-profile, and tau-profile based on the output of \code{ObsAsy3D} using the ggplot2 package.\cr
#' Note, it will only show the confidence interval of \code{'Asymptotic'} when setting \code{method = "both"} in \code{ObsAsy3D}.
#' 
#' @param output the output of the function \code{ObsAsy3D}.\cr
#' @param profile a selection of profile versus to diversity. User can choose \code{'q'}, \code{'time'}, and \code{'tau'}. Default is \code{'q'} profile. \code{'time'} profile for only when \code{diversity = "PD"}. \code{'tau'} profile for only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}.\cr
#' @return a figure of asymptotic and observed diversity in q-profile, time-profile, or tau-profile.\cr\cr
#'
#' @examples
#' \donttest{
#' # diversity = 'TD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' output_ObsAsy_TD_abun <- ObsAsy3D(Brazil_rainforest_abun_data, diversity = 'TD', 
#'                                   datatype = "abundance")
#' ggObsAsy3D(output_ObsAsy_TD_abun)
#' 
#' 
#' # diversity = 'PD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_phylo_tree)
#' data <- Brazil_rainforest_abun_data
#' tree <- Brazil_rainforest_phylo_tree
#' output_ObsAsy_PD_abun <- ObsAsy3D(data, diversity = 'PD', q = seq(0, 2, by = 0.5), 
#'                                   datatype = "abundance", nboot = 20, PDtree = tree, 
#'                                   PDtype = "meanPD")
#' ggObsAsy3D(output_ObsAsy_PD_abun, profile = "q")
#' 
#' 
#' # diversity = 'FD' for abundance data
#' data(Brazil_rainforest_abun_data)
#' data(Brazil_rainforest_distance_matrix)
#' data <- Brazil_rainforest_abun_data
#' distM <- Brazil_rainforest_distance_matrix
#' output_ObsAsy_FD_abun <- ObsAsy3D(data, diversity = 'FD', q = seq(0, 2, 0.5), 
#'                                   datatype = "abundance", nboot = 10, 
#'                                   FDdistM = distM, FDtype = 'AUC')
#' ggObsAsy3D(output_ObsAsy_FD_abun)
#' 
#' 
#' # diversity = 'TD' for incidence data
#' data(Fish_incidence_data)
#' output_ObsAsy_TD_inci <- ObsAsy3D(Fish_incidence_data, diversity = 'TD', datatype = "incidence_raw")
#' ggObsAsy3D(output_ObsAsy_TD_inci)
#' 
#' 
#' # diversity = 'PD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_phylo_tree)
#' data <- Fish_incidence_data
#' tree <- Fish_phylo_tree
#' output_ObsAsy_PD_inci <- ObsAsy3D(data, diversity = 'PD', q = c(0, 1, 2), 
#'                                   datatype = "incidence_raw", nboot = 20, PDtree = tree, 
#'                                   PDreftime = seq(0.1, 1, length.out = 40), method = "Observed")
#' ggObsAsy3D(output_ObsAsy_PD_inci, profile = "time")
#' 
#' 
#' # diversity = 'FD' for incidence data
#' data(Fish_incidence_data)
#' data(Fish_distance_matrix)
#' data <- Fish_incidence_data
#' distM <- Fish_distance_matrix
#' output_ObsAsy_FD_inci <- ObsAsy3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, 
#'                                   FDdistM = distM, FDtype = 'AUC')
#' ggObsAsy3D(output_ObsAsy_FD_inci)
#' }
#' 
#' 
#' @export
ggObsAsy3D <- function(output, profile = 'q'){
  
  if (sum(unique(output$Method) %in% c("Asymptotic", "Observed")) == 0)
    stop("Please use the output from specified function 'ObsAsy3D'")
  
  if (!(profile %in% c('q', 'time', 'tau')))
    stop("Please select one of 'q', 'time', 'tau' profile.")
  
  if (sum(colnames(output)[1:7] == c('Assemblage', 'Order.q', 'qTD', 's.e.', 'qTD.LCL', 'qTD.UCL', 'Method')) == 7) {
    
    class = 'TD'
    
  } else if (sum(colnames(output)[1:7] == c('Assemblage', 'Order.q', 'qPD', 's.e.', 'qPD.LCL', 'qPD.UCL', 'Method')) == 7) {
    
    class = 'PD'
    
  } else if (sum(colnames(output)[1:7] == c('Assemblage', 'Order.q', 'qFD', 's.e.', 'qFD.LCL', 'qFD.UCL', 'Method')) == 7) {
    
    if ("Tau" %in% colnames(output)) class = 'FD' else class = 'FD(AUC)'
    
  } else {stop("Please use the output from specified function 'ObsAsy3D'")}
  
  ## TD & q-profile ##
  if (class == 'TD') {
    out = ggplot(output, aes(x = Order.q, y = qTD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qTD.LCL, ymax = qTD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method == 'Asymptotic')) out = out + labs(x = 'Order q', y = 'Asymptotic taxonomic diversity')
      if (unique(output$Method == 'Observed')) out = out + labs(x = 'Order q', y = 'Observed taxonomic diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qTD.LCL, ymax = qTD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Taxonomic diversity')
    }
  }
  
  ## PD & q-profile ##
  if (class == 'PD' & profile == 'q') {
    
    output$Reftime = paste('Reftime = ', round(output$Reftime, 3), sep = '')
    out = ggplot(output, aes(x = Order.q, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method) == 'Asymptotic' & unique(output$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Asymptotic phylogenetic diversity')
      if (unique(output$Method) == 'Observed' & unique(output$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Observed phylogenetic diversity')
      if (unique(output$Method) == 'Asymptotic' & unique(output$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Asymptotic mean phylogenetic diversity')
      if (unique(output$Method) == 'Observed' & unique(output$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Observed mean phylogenetic diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(output$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Phylogenetic diversity')
      if (unique(output$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Mean phylogenetic diversity')
      
    }
    
    out = out + facet_grid(.~Reftime, scales = "free_y")
  }
  
  ## PD & time-profile ##
  if (class == 'PD' & profile == 'time') {
    
    output$Order.q = paste('q = ', output$Order.q, sep = '')
    out = ggplot(output, aes(x = Reftime, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method) == 'Asymptotic' & unique(output$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Asymptotic phylogenetic diversity')
      if (unique(output$Method) == 'Observed' & unique(output$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Observed phylogenetic diversity')
      if (unique(output$Method) == 'Asymptotic' & unique(output$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Asymptotic mean phylogenetic diversity')
      if (unique(output$Method) == 'Observed' & unique(output$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Observed mean phylogenetic diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(output$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Phylogenetic diversity')
      if (unique(output$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Mean phylogenetic diversity')
      
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## FD & q-profile ##
  if (class == 'FD' & profile == 'q') {
    
    output$Tau = paste('Tau = ', round(output$Tau, 3), sep = '')
    out = ggplot(output, aes(x = Order.q, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic functional diversity')
      if (unique(output$Method) == 'Observed') out = out + labs(x = 'Order q', y = 'Observed functional diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Functional diversity')
    }
    out = out + facet_grid(.~Tau, scales = "free_y")
  }
  
  ## FD & tau-profile ##
  if (class == 'FD' & profile == 'tau') {
    
    output$Order.q = paste('Order q = ', output$Order.q, sep = '')
    out = ggplot(output, aes(x = Tau, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method) == 'Asymptotic') out = out + labs(x = 'Tau', y = 'Asymptotic functional diversity')
      if (unique(output$Method) == 'Observed') out = out + labs(x = 'Tau', y = 'Observed functional diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Tau', y = 'Functional diversity')
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## AUC & q-profile ##
  if (class == 'FD(AUC)') {
    
    out = ggplot(output, aes(x = Order.q, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(output$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(output$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic Functional diversity (AUC)')
      if (unique(output$Method) == 'Observed') out = out + labs(x = 'Order q', y = 'Observed Functional diversity (AUC)')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = output %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Functional diversity (AUC)')
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


## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c(".", "AUC_L", "AUC_R", "AUC_se", "Assemblage", "Inode", 
                         "Method", "Order.q", "PD.obs", "Q1*", "Q2*", "R1", "R2", 
                         "Reftime", "S.obs", "SC", "SC(2T)", "SC(2n)", "SC(T)", 
                         "SC(n)", "SC.LCL", "SC.UCL", "SC.s.e.", "SC_se", "Tau", 
                         "Type", "U", "ai", "branch.height", "branch.length", 
                         "branch.length.new", "cumsum.length", "e", "edgelengthv", 
                         "f1*", "f2*", "g1", "g2", "label", "label.new", "length.new", 
                         "m", "nT", "newlabel", "newlable", "node", "node.age", "qFD", 
                         "qFD.LCL", "qFD.UCL", "qPD", "qPD.LCL", "qPD.UCL", "qTD", 
                         "qTD.LCL", "qTD.UCL", "quantile", "refT", "s.e.", "se", 
                         "species", "tgroup", "tmp", "vi", "x", "y"
))




