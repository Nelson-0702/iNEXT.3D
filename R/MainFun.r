#' Data Information for three diversity.
#' 
#' \code{DataInfo3D} Exhibit basic data information
#' 
#' @param data a \code{matrix}, \code{data.frame} (species by sites), or \code{list} of species abundance/incidence frequencies.\cr 
#' If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
#' @param diversity selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage. \cr
#' If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. \cr
#' It is necessary when \code{diversity = 'PD'} and \code{datatype = "incidence_raw"}.
#' @param PDtree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. It is necessary when \code{diversity = 'PD'}.
#' @param PDreftime Select several reference time points for \code{diversity = 'PD'}. Default is NULL. 
#' @param FDdistM a pair wise distance matrix for all pairs of observed species in the pooled assemblage. It will be use when \code{diversity = 'FD'}.
#' @param FDtype a binary selection for FD. \code{FDtype = "tau_values"} computes diversity under certain threshold values. \code{FDtype = "AUC"} computes an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}. It will be use when \code{diversity = 'FD'} and \code{FDtype = "tau_values"}.
#' 
#' @return a data.frame of basic data information including assemblage name (Assemblage), sample size (n) or total sampling units (T), total incidences (U), observed species richness (S.obs), sample coverage estimate (SC).\cr\cr
#' Besides, show the first ten species abundance (or incidence) frequency counts in the reference sample in TD. (f1-f10 or Q1-Q10)\cr\cr
#' In PD, show the the observed total branch length in the phylogenetic tree (PD.obs), the number of singletons and doubletons in the node/branch set (f1*-f2*), the total branch length of those singletons/doubletons in the node/branch set (g1-g2), reference time (Reftime).\cr\cr
#' In FD (FDtype = "tau_values"), show the number of singletons and doubletons in the functional group (a1*-a2*), the total contribution of those singletons/doubletons in the functional group (h1-h2), the threshold of functional distinctiveness between any two species (Tau).\cr\cr
#' In FD (FDtype = "AUC"), show the the minimum distance among all non-zero elements in the distance matrix (dmin), the mean distance between any two individuals randomly selected from the pooled data (dmean), the maximum distance among all elements in the distance matrix (dmax).\cr
#' 
#' @examples
#' ## example for abundance-based data
#' # diversity = 'TD'
#' data(dunes)
#' DataInfo3D(dunes$data, diversity = 'TD', datatype = "abundance")
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' DataInfo3D(data, diversity = 'PD', datatype = "abundance", PDtree = tree)
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' DataInfo3D(data, diversity = 'FD', datatype = "abundance", FDdistM = distM, FDtype = 'tau_values')
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' DataInfo3D(data, diversity = 'FD', datatype = "abundance", FDdistM = distM, FDtype = 'AUC')
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' DataInfo3D(fish$data, diversity = 'TD', datatype = "incidence_raw")
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' DataInfo3D(data, diversity = 'PD', datatype = "incidence_raw", nT = nT, PDtree = tree)
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' DataInfo3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'            FDdistM = distM, FDtype = 'tau_values')
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' DataInfo3D(data, diversity = 'FD', datatype = "incidence_raw", FDdistM = distM, FDtype = 'AUC')
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
    
    if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    
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
        select(Assemblage, n, S.obs, SC, PD.obs, `f1*`, `f2*`, g1, g2, Reftime)
      
    }else if (datatype=='incidence_raw'){
      
      out <- lapply(mydata, function(x){
        datainf(data = x, datatype, phylotr = mytree,reft = PDreftime) %>% mutate(Reftime = PDreftime)
      }) %>% do.call(rbind,.) %>% 
        mutate(Assemblage = rep(names(mydata), each = length(PDreftime))) %>%
        select(Assemblage,`T`, U, S.obs, SC, PD.obs, `Q1*`, `Q2*`, R1, R2, Reftime)
      
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
    
    
    # out = lapply(FDtau, function(tau) {
    #   
    #   out <- TDinfo(lapply(dat, function(x) data_transform(x, distM, tau, datatype, integer = TRUE)$ai[,1]), datatype)
    #   
    #   if (datatype == "abundance") out$n = sapply(dat, function(x) sum(x)) else if (datatype == "incidence_freq") 
    #     out$U = sapply(dat, function(x) sum(x[-1]))
    #   
    #   out$SC = sapply(dat, function(x) {
    #     
    #     if (datatype == "abundance") {
    #       
    #       n = sum(x)
    #       f1 = sum(x == 1)
    #       f2 = sum(x == 2)
    #       f0.hat <- ifelse(f2 == 0, (n-1) / n * f1 * (f1-1) / 2, (n-1) / n * f1^2 / 2 / f2) 
    #       A <- ifelse(f1 > 0, n * f0.hat / (n * f0.hat + f1), 1)
    #       Chat <- 1 - f1/n * A
    #       
    #     } else if (datatype == "incidence_freq") {
    #       
    #       nT = x[1]
    #       x = x[-1]
    #       U <- sum(x)
    #       Q1 = sum(x == 1)
    #       Q2 = sum(x == 2)
    #       Q0.hat <- ifelse(Q2 == 0, (nT-1) / nT * Q1 * (Q1-1) / 2, (nT-1) / nT * Q1^2 / 2 / Q2) 
    #       A <- ifelse(Q1 > 0, nT * Q0.hat / (nT * Q0.hat + Q1), 1)
    #       Chat <- 1 - Q1/U * A
    #       
    #     }
    #     
    #     Chat
    #   })
    #   
    #   colnames(out)[colnames(out) %in% paste0("f", 1:10)] = paste0("a", 1:10, "'")
    #   
    #   out$Tau = tau
    #   
    #   return(out)
    #   
    # }) %>% do.call(rbind,.)
    
    
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
          
          multiple = tibble('Assemblage' = names(dat)[i], 
                            'n' = n, 
                            'S.obs' = sum(x > 0), 
                            'SC' = Chat, 
                            'a1*' = sum(aivi$ai == 1), 'a2*' = sum(aivi$ai == 2), 
                            'h1' = sum(aivi$vi[aivi$ai == 1,]), 'h2' = sum(aivi$vi[aivi$ai == 2,]),
                            'Tau' = tau)
          
        } else if (datatype == "incidence_freq") {
          
          nT = x[1]
          x = x[-1]
          U <- sum(x)
          Q1 = sum(x == 1)
          Q2 = sum(x == 2)
          Q0.hat <- ifelse(Q2 == 0, (nT-1) / nT * Q1 * (Q1-1) / 2, (nT-1) / nT * Q1^2 / 2 / Q2) 
          A <- ifelse(Q1 > 0, nT * Q0.hat / (nT * Q0.hat + Q1), 1)
          Chat <- 1 - Q1/U * A
          
          multiple = tibble('Assemblage' = names(dat)[i], 
                            'T' = nT, 
                            'U' = U,
                            'S.obs' = sum(x > 0), 
                            'SC' = Chat, 
                            'a1*' = sum(aivi$ai == 1), 'a2*' = sum(aivi$ai == 2), 
                            'h1' = sum(aivi$vi[aivi$ai == 1,]), 'h2' = sum(aivi$vi[aivi$ai == 2,]),
                            'Tau' = tau)
        }
        
        return(multiple)
        
      }) %>% do.call(rbind,.)
      
      
    }) %>% do.call(rbind,.)
    
    
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
        
      }else if(datatype=='incidence_freq'){
        
        tmp <- matrix(i[-1]/sum(i[-1]), ncol = 1)
        
      }
      
      dmean <- sum ( (tmp %*% t(tmp) ) * distM)
      distM <- distM[tmp > 0, tmp > 0]
      dmin <- min(distM[lower.tri(distM)])
      dmax <- max(distM[distM > 0])
      
      c(dmin, dmean, dmax)
    }))
    
    if (datatype == "abundance") {
      
      out <- cbind(TDinfo(dat, datatype)[,1:4], Tau)
      colnames(out)[5:7] = c("dmin", "dmean", "dmax")
      rownames(out) = NULL
      
    } else {
      
      out <- cbind(TDinfo(dat, datatype)[,1:5], Tau)
      colnames(out)[6:8] = c("dmin", "dmean", "dmax")
      rownames(out) = NULL
      
    }
    
    return(out)
  }
  
  return(out)
}


#' @useDynLib iNEXT.3D, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL



#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT3D}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_freq"}, data can be input as a vector of incidence frequencies (for a single assemblage), matrix/data.frame (species by assemblages), or a list of incidence frequencies; the first entry in all types of input must be the number of sampling units in each assemblage. \cr
#' (c) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}), or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
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
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified, then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (required only when \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).  
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
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
#' @return a list of three objects: \code{$DataInfo} (or \code{$PDInfo}, \code{$FDInfo}, \code{$AUCInfo}) for summarizing data information; 
#' \code{$iNextEst} (or \code{$PDiNextEst}, \code{$FDiNextEst}, \code{$AUCiNextEst}) for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} (or \code{$PDAsyEst}, \code{$FDAsyEst}, \code{$AUCAsyEst}) for showing asymptotic diversity estimates along with related statistics.  
#' 
#' @examples
#' ## example for abundance based data (list of vector)
#' # diversity = 'TD'
#' data(dunes)
#' out1 <- iNEXT3D(dunes$data, diversity = 'TD', q = c(0,1,2), datatype = "abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' out2 <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "abundance", nboot = 30, PDtree = tree)
#' out2
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out3 <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 30, FDdistM = distM, FDtype = 'tau_values')
#' out3
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out4 <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 0, FDdistM = distM)
#' out4
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' out5 <- iNEXT3D(fish$data, diversity = 'TD', q = 1, datatype = "incidence_raw")
#' out5
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", nT = nT, PDtree = tree)
#' out6
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out7 <- iNEXT3D(data, diversity = 'FD', datatype = "incidence_raw", FDdistM = distM, FDtype = 'tau_values')
#' out7
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out8 <- iNEXT3D(data, diversity = 'FD', FDdistM = distM, datatype = "incidence_raw", nboot = 0)
#' out8
#' 
#' @references
#' Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.
#' @export
#' 
iNEXT3D <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, nboot = 50, conf = 0.95, nT = NULL, 
                    PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL) {
  
  if ( !(diversity %in% c('TD', 'PD', 'FD')) ) 
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)
  
  if (diversity == 'TD') {
    
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
    
    z <- qnorm(1-(1-0.95)/2)
    
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
    LCL <- index$qD.LCL[index$Method=='Asymptotic']
    UCL <- index$qD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    
    
    out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
    out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
    out$size_based <- as_tibble(out$size_based)
    out$coverage_based <- as_tibble(out$coverage_based)
    
    info <- DataInfo3D(data, diversity = 'TD', datatype, nT)
    
    
    out <- list("DataInfo"=info, "iNextEst"=out, "AsyEst"=index)
  } 
  
  if (diversity == 'PD') {
    
    if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    
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
    if (nboot > 0) index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('q = 0 PD','q = 1 PD','q = 2 PD')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Phylogenetic Diversity", "Phylogenetic Observed", "Phylogenetic Estimator", "s.e.", "LCL", "UCL")
    index$Reftime = PDreftime
    index$Type = PDtype
    
    info <- DataInfo3D(data, diversity = 'PD', datatype = datatype, nT, PDtree = PDtree, PDreftime = PDreftime)
    
    
    out = list("PDInfo"=info, "PDiNextEst"=out, "PDAsyEst"=index)
    
  } 
  
  if (diversity == 'FD' & FDtype == 'tau_values') {
    
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
        if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nt'
        
        ## coverage-based
        temp2 <- lapply(1:length(dat), function(i) invChatFD(datalist = dat[i], dij = dist, q = q, datatype = datatype,
                                                             level = unique(Coverage(data = dat[[i]], datatype = datatype, m = size[[i]])), 
                                                             nboot = nboot, conf = conf, tau = FDtau)) %>% do.call(rbind,.)
        temp2$qFD.LCL[temp2$qFD.LCL<0] <- 0
        
        if (datatype == 'incidence_freq') colnames(temp2)[colnames(temp2) == 'm'] = 'nt'
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
    if (nboot > 0) index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('q = 0 FD(single tau)','q = 1 FD(single tau)','q = 2 FD(single tau)')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Functional Diversity", "Functional Observed", "Functional Estimator", "s.e.", "LCL", "UCL")
    index$Tau = FDtau
    
    info <- DataInfo3D(data, diversity = 'FD', datatype = datatype, FDdistM = FDdistM, FDtype = 'tau_values', FDtau = FDtau, nT = nT)
    out = list("FDInfo" = info, "FDiNextEst" = out, "FDAsyEst" = index)
    
  } 
  
  if (diversity == 'FD' & FDtype == 'AUC') {
   
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
    
    
    FUN <- function(e){
      if(inherits(dat, "list")){
        ## size-based
        temp1 = AUCtable_iNextFD(datalist = dat, dij = dist, q = q, datatype = datatype,
                                 tau = NULL, nboot = nboot, conf = conf, m = size)
        temp1$qAUC.LCL[temp1$qAUC.LCL<0] <- 0; temp1$SC.LCL[temp1$SC.LCL<0] <- 0
        temp1$SC.UCL[temp1$SC.UCL>1] <- 1
        if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nt'
        
        ## coverage-based
        temp2 <- lapply(1:length(dat), function(i) AUCtable_invFD(datalist = dat[i], dij = dist, q = q, datatype = datatype,
                                                                  level = unique(Coverage(data = dat[[i]], datatype = datatype, m = size[[i]])), 
                                                                  nboot = nboot, conf = conf, tau = NULL)) %>% do.call(rbind,.)
        temp2$qAUC.LCL[temp2$qAUC.LCL<0] <- 0
        if (datatype == 'incidence_freq') colnames(temp2)[colnames(temp2) == 'm'] = 'nt'
        
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
    index <- rbind(AUCtable_est(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype,
                                nboot = nboot, conf = conf, tau = NULL),
                   AUCtable_mle(datalist = dat, dij = dist, q = c(0, 1, 2), datatype = datatype,
                                nboot = nboot, conf = conf, tau = NULL))
    index = index[order(index$Assemblage),]
    LCL <- index$qAUC.LCL[index$Method=='Asymptotic']
    UCL <- index$qAUC.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qAUC')
    index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
    if (nboot > 0) index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('q = 0 FD(AUC)','q = 1 FD(AUC)','q = 2 FD(AUC)')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Functional Diversity", "Functional Observed", "Functional Estimator", "s.e.", "LCL", "UCL")
    
    info <- DataInfo3D(data, diversity = 'FD', datatype = datatype, FDdistM = FDdistM, FDtype = 'AUC', nT = nT)
    out = list("AUCInfo" = info, "AUCiNextEst" = out, "AUCAsyEst" = index)
    
  }
  
  class(out) <- c("iNEXT3D")
  
  return(out)
}


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
#' @return a ggplot2 object for coverage-based or size-based rarefaction and extrapolation
#' 
#' @examples
#' ## example for abundance based data (list of vector)
#' # diversity = 'TD'
#' data(dunes)
#' out1 <- iNEXT3D(dunes$data, diversity = 'TD', q = c(0,1,2), datatype = "abundance")
#' ggiNEXT3D(out1, facet.var = "Assemblage")
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' out2 <- iNEXT3D(data, diversity = 'PD', q = c(0,1,2), datatype = "abundance", nboot = 30, PDtree = tree)
#' ggiNEXT3D(out2, type = c(1, 3))
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out3 <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 0, FDdistM = distM, FDtype = 'tau_values')
#' ggiNEXT3D(out3)
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out4 <- iNEXT3D(data, diversity = 'FD', datatype = "abundance", nboot = 0, FDdistM = distM)
#' ggiNEXT3D(out4)
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' out5 <- iNEXT3D(fish$data, diversity = 'TD', q = 1, datatype = "incidence_raw")
#' ggiNEXT3D(out5)
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- iNEXT3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", nT = nT, PDtree = tree)
#' ggiNEXT3D(out6, facet.var = "Order.q", color.var = "Assemblage")
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out7 <- iNEXT3D(data, diversity = 'FD', datatype = "incidence_raw", FDdistM = distM, FDtype = 'tau_values')
#' ggiNEXT3D(out7)
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out8 <- iNEXT3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, FDdistM = distM)
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
    
  } else if (sum(names(outcome) %in% c("AUCInfo", 'AUCiNextEst', 'AUCAsyEst')) == 3) {
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
  if (length(type) == 1) out = out[[1]]
  
  return(out)
}


type_plot = function(x_list, type, class, datatype, facet.var, color.var) {
  x_name <- colnames(x_list$size_based)[2]
  xlab_name <- ifelse(datatype == "incidence", "sampling units", "individuals")
  
  if (class == 'TD') {
    ylab_name = "Taxonomic diversity"
  } else if (class == 'FD') {
    ylab_name = "Functional diversity"
  } else if (class == 'AUC') {
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
    output <- x_list$coverage_based %>% data.frame
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


#' Compute species diversity with a particular of sample size/coverage 
#' 
#' \code{estimate3D}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_freq"}, data can be input as a vector of incidence frequencies (for a single assemblage), matrix/data.frame (species by assemblages), or a list of incidence frequencies; the first entry in all types of input must be the number of sampling units in each assemblage. \cr
#' (c) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}), or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection)
#' @param base selection of sample-size-based (\code{base = "size"}) or coverage-based (\code{base = "coverage"}) rarefaction and extrapolation.
#' @param level A numerical vector specifying the particular sample sizes or sample coverages (between 0 and 1). \cr
#' If \code{base = "coverage"} (default) and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes. \cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes. 
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified, then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (required only when \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).  
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' 
#' @return a \code{data.frame} of diversity table including the following arguments:
#' 'Assemblage' = the assemblage name.\cr\cr
#' 'm' (or 'nt') = the corresponding sample size (or sampling units) for the standardized coverage value. \cr\cr
#' 'Method' = Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample.\cr\cr
#' 'Order.q' = the diversity order of q.\cr\cr
#' 'SC' = the target standardized coverage value. \cr\cr
#' 'qD' (or 'qPD', 'qFD', 'qAUC') = the estimated diversity of order q for the target coverage value. The estimate for complete coverage (or size = infinity) represents the estimated asymptotic diversity. \cr\cr
#' 's.e.' = standard error of diversity estimate.\cr\cr
#' 'qD.LCL' (or 'qPD.LCL', 'qFD.LCL', 'qAUC.LCL'), 'qD.UCL' (or 'qPD.UCL', 'qFD.UCL', 'qAUC.UCL') = the bootstrap lower and upper confidence limits for the diversity of order q at the specified level (with a default value of 0.95).\cr\cr
#' 'Reftime' = reference times for PD.\cr\cr
#' 'Type' = "PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages).\cr\cr
#' 'Tau' = the threshold of functional distinctiveness between any two species.\cr
#' 
#' 
#' @examples
#' # diversity = 'TD'
#' data(dunes)
#' out1 <- estimate3D(dunes$data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", base = "size")
#' out1
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' out2 <- estimate3D(data, diversity = 'PD', datatype = "abundance", base = "coverage", PDtree = tree)
#' out2
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out3 <- estimate3D(data, diversity = 'FD', datatype = "abundance", base = "size", FDdistM = distM, FDtype = 'tau_values')
#' out3
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out4 <- estimate3D(data, diversity = 'FD', datatype = "abundance", base = "coverage", nboot = 0, FDdistM = distM)
#' out4
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' out5 <- estimate3D(fish$data, diversity = 'TD', q = c(0,1,2), datatype = "incidence_raw", base = "coverage", level=0.985)
#' out5
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- estimate3D(data, diversity = 'PD', datatype = "incidence_raw", base = "size", nT = nT, PDtree = tree)
#' out6
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out7 <- estimate3D(data, diversity = 'FD', datatype = "incidence_raw", base = "coverage", FDdistM = distM, FDtype = 'tau_values')
#' out7
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out8 <- estimate3D(data, diversity = 'FD', datatype = "incidence_raw", base = "size", nboot = 20, FDdistM = distM)
#' out8
#' 
#' @references
#' Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.
#' 
#' @export
estimate3D <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot = 50, conf = 0.95, nT = NULL, 
                       PDtree, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL) {
  
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
    out$qD.LCL[out$qD.LCL<0] <- 0
    
  } 
  
  if (diversity == 'PD') {
    
    if(datatype == "incidence_freq") stop ("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    
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
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
      out = out %>% .[,c(1:4, 9, 5:8, 10)]
      
    } else if (base == "coverage") {
      
      out <- invChatFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                       level = level, nboot = nboot, conf = conf, tau = FDtau)
      out$qFD.LCL[out$qFD.LCL<0] <- 0
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
      
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
    
    
    
    if (base == 'size') {
      
      out = AUCtable_iNextFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                             tau = NULL, nboot = nboot, conf = conf, m = lapply(1:length(dat), function(i) level)) %>% 
        select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
      out$qAUC.LCL[out$qAUC.LCL<0] <- 0
      # out$SC.LCL[out$SC.LCL<0] <- 0
      # out$SC.UCL[out$SC.UCL>1] <- 1
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
      out = out %>% .[,c(1:4, 9, 5:8)]
      
    } else if (base == 'coverage') {
      
      out <- AUCtable_invFD(datalist = dat, dij = FDdistM, q = q, datatype = datatype,
                            level = level, nboot = nboot, conf = conf, tau = NULL)
      if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
      
    }
    
  }
  
  return(out)
}


#' Asymptotic diversity and observed diversity of q profile 
#' 
#' \code{AO3D} The estimated diversity of order q 
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a vector of species abundances (for a single assemblage), matrix/data.frame (species by assemblages), or a list of species abundance vectors. \cr
#' (b) For \code{datatype = "incidence_freq"}, data can be input as a vector of incidence frequencies (for a single assemblage), matrix/data.frame (species by assemblages), or a list of incidence frequencies; the first entry in all types of input must be the number of sampling units in each assemblage. \cr
#' (c) For \code{datatype = "incidence_raw"}, data can be input as a list of matrix/data.frame (species by sampling units); data can also be input as a matrix/data.frame by merging all sampling units across assemblages based on species identity; in this case, the number of sampling units (nT, see below) must be input. 
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is seq(0, 2, by = 0.2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}), sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}), or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection)
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param nT (required only when \code{datatype = "incidence_raw"} and input data is matrix/data.frame) a vector of nonnegative integers specifying the number of sampling units in each assemblage. If assemblage names are not specified, then assemblages are automatically named as "assemblage1", "assemblage2",..., etc. 
#' @param method computing type. Select 'Asymptotic' or 'Observed'.
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (required only when \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).  
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' 
#' @return a table of diversity table including the following arguments.
#' 'Order.q' = the diversity order of q.\cr\cr
#' 'qD' (or 'qPD', 'qFD', 'qAUC') = the estimated asymptotic diversity or empirical (observed) diversity of order q. \cr\cr
#' 's.e.' = standard error of diversity. \cr\cr
#' 'qD.LCL' (or 'qPD.LCL', 'qFD.LCL', 'qAUC.LCL'), 'qD.UCL' (or 'qPD.UCL', 'qFD.UCL', 'qAUC.UCL') = the bootstrap lower and upper confidence limits for the diversity of order q at the specified level (with a default value of 0.95).\cr\cr
#' 'Assemblage' = the assemblage name.\cr\cr
#' 'Method' = "Asymptotic" or "Empirical".\cr\cr
#' 'Reftime' = reference times for PD.\cr\cr
#' 'Type' = "PD" (effective total branch length) or "meanPD" (effective number of equally divergent lineages).\cr\cr
#' 'Tau' = the threshold of functional distinctiveness between any two species.\cr
#' 
#' 
#' @examples
#' ## example for abundance-based data
#' # diversity = 'TD'
#' data(dunes)
#' out1 <- AO3D(dunes$data, diversity = 'TD', datatype = "abundance")
#' out1
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' out2 <- AO3D(data, diversity = 'PD', q = seq(0, 2, by = 0.25), 
#'              datatype = "abundance", nboot = 30, PDtree = tree)
#' out2
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out3 <- AO3D(data, diversity = 'FD', datatype = "abundance", 
#'              nboot = 50, FDdistM = distM, FDtype = 'tau_values')
#' out3
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out4 <- AO3D(data[,1:2], diversity = 'FD', q = seq(0, 2, 0.5), 
#'              datatype = "abundance", nboot = 20, FDdistM = distM)
#' out4
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' out5 <- AO3D(fish$data, diversity = 'TD', datatype = "incidence_raw")
#' out5
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- AO3D(data, diversity = 'PD', q = seq(0, 2, by = 0.25), 
#'              datatype = "incidence_raw", nT = nT, PDtree = tree)
#' out6
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out7 <- AO3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'              FDdistM = distM, FDtype = 'tau_values')
#' out7
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out8 <- AO3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, FDdistM = distM)
#' out8
#' 
#' @references
#' Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.
#' 
#' @export
AO3D <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT = NULL, method = c('Asymptotic', 'Observed'),
                 PDtree, PDreftime = NULL, PDtype = 'meanPD', FDdistM, FDtype = 'AUC', FDtau = NULL) {
  
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
    
    if (datatype == "incidence_freq") 
      stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
    
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
    
    if (sum(method == "Asymptotic") == length(method)) 
      out = AUCtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                         nboot = nboot, conf = conf,  tau = NULL) else if (sum(method == "Observed") == length(method)) 
                           
                           out = AUCtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                              nboot = nboot, conf = conf, tau = NULL) else if (sum(method == c("Asymptotic", "Observed")) == length(method)) 
                                                
                                                out = rbind(AUCtable_est(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                         nboot = nboot, conf = conf, tau = NULL), 
                                                            AUCtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype, 
                                                                         nboot = nboot, conf = conf, tau = NULL))
    
  }
  
  return(out)
}


#' ggplot for Asymptotic diversity
#'
#' \code{ggAO3D} Plots q-profile, time-profile, and tau-profile based on the outcome of \code{AO3D} using the ggplot2 package.\cr
#' It will only show the confidence interval of 'Estimated'.
#' 
#' @param outcome the outcome of the functions \code{AO3D}.\cr
#' @param profile a selection of profile versus to diversity. User can choose \code{'q'}, \code{'time'}, and \code{'tau'}. Default is \code{'q'} profile. \code{'time'} profile for only when \code{diversity = "PD"}. \code{'tau'} profile for only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}.\cr
#' @return a figure of asymptotic or empirical (observed) diversity in q-profile, time-profile, or tau-profile.\cr\cr
#'
#' @examples
#' ## example for abundance-based data
#' # diversity = 'TD'
#' data(dunes)
#' out1 <- AO3D(dunes$data, diversity = 'TD', datatype = "abundance")
#' ggAO3D(out1)
#' 
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' tree <- dunes$tree
#' out2 <- AO3D(data, diversity = 'PD', q = seq(0, 2, by = 0.25), datatype = "abundance", 
#'              nboot = 30, PDtree = tree, PDtype = "meanPD")
#' ggAO3D(out2, profile = "q")
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out3 <- AO3D(data, diversity = 'FD', q = c(0, 1, 2), datatype = "abundance", nboot = 0, 
#'              FDtau = seq(0, 0.6, 0.1), FDdistM = distM, FDtype = 'tau_values')
#' ggAO3D(out3, profile = "tau")
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' distM <- dunes$dist
#' out4 <- AO3D(data, diversity = 'FD', q = seq(0, 2, 0.5), datatype = "abundance", nboot = 0, FDdistM = distM)
#' ggAO3D(out4)
#' 
#' ## example for incidence-based data
#' # diversity = 'TD'
#' data(fish)
#' out5 <- AO3D(fish$data, diversity = 'TD', datatype = "incidence_raw")
#' ggAO3D(out5)
#' 
#' # diversity = 'PD'
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' nT <- data.inc$nT
#' out6 <- AO3D(data, diversity = 'PD', q = c(0, 1, 2), datatype = "incidence_raw", 
#'              nT = nT, PDtree = tree, PDreftime = seq(0.1, 82.8575, length.out = 40))
#' ggAO3D(out6, profile = "time")
#' 
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out7 <- AO3D(data, diversity = 'FD', datatype = "incidence_raw", 
#'              FDdistM = distM, FDtype = 'tau_values')
#' ggAO3D(out7)
#' 
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(fish)
#' data <- fish$data
#' distM <- fish$dist
#' out8 <- AO3D(data, diversity = 'FD', datatype = "incidence_raw", nboot = 20, FDdistM = distM)
#' ggAO3D(out8)
#'
#' @export
ggAO3D <- function(outcome, profile = 'q'){
  if (sum(unique(outcome$Method) %in% c("Asymptotic", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'AO3D'")
  
  if (!(profile %in% c('q', 'time', 'tau')))
    stop("Please select one of 'q', 'time', 'tau' profile.")
  
  if (sum(colnames(outcome)[1:7] == c('Order.q', 'qD', 's.e.', 'qD.LCL', 'qD.UCL', 'Assemblage', 'Method')) == 7) {
    class = 'TD'
  } else if (sum(colnames(outcome)[1:7] == c('Order.q', 'qPD', 's.e.', 'qPD.LCL', 'qPD.UCL', 'Assemblage', 'Method')) == 7) {
    class = 'PD'
  } else if (sum(colnames(outcome)[1:7] == c('Order.q', 'qFD', 's.e.', 'qFD.LCL', 'qFD.UCL', 'Assemblage', 'Method')) == 7) {
    class = 'FD'
  } else if (sum(colnames(outcome)[1:7] == c('Order.q', 'qAUC', 's.e.', 'qAUC.LCL', 'qAUC.UCL', 'Assemblage', 'Method')) == 7) {
    class = 'AUC'
  } else {stop("Please use the outcome from specified function 'AO3D'")}
  
  ## TD & q-profile ##
  if (class == 'TD') {
    out = ggplot(outcome, aes(x = Order.q, y = qD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method == 'Asymptotic')) out = out + labs(x = 'Order q', y = 'Asymptotic taxonomic diversity')
      if (unique(outcome$Method == 'Empirical')) out = out + labs(x = 'Order q', y = 'Empirical taxonomic diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qD.LCL, ymax = qD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Taxonomic diversity')
    }
  }
  
  ## PD & q-profile ##
  if (class == 'PD' & profile == 'q') {
    
    outcome$Reftime = paste('Reftime = ', round(outcome$Reftime, 3), sep = '')
    out = ggplot(outcome, aes(x = Order.q, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Asymptotic phylogenetic diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Empirical phylogenetic diversity')
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Asymptotic mean phylogenetic diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Empirical mean phylogenetic diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Type) == 'PD') out = out + labs(x = 'Order q', y = 'Phylogenetic diversity')
      if (unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Order q', y = 'Mean phylogenetic diversity')
      
    }
    
    out = out + facet_grid(.~Reftime, scales = "free_y")
  }
  
  ## PD & time-profile ##
  if (class == 'PD' & profile == 'time') {
    
    outcome$Order.q = paste('q = ', outcome$Order.q, sep = '')
    out = ggplot(outcome, aes(x = Reftime, y = qPD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qPD.LCL, ymax = qPD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Asymptotic phylogenetic diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Empirical phylogenetic diversity')
      if (unique(outcome$Method) == 'Asymptotic' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Asymptotic mean phylogenetic diversity')
      if (unique(outcome$Method) == 'Empirical' & unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Empirical mean phylogenetic diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qPD.LCL, ymax = qPD.UCL), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Type) == 'PD') out = out + labs(x = 'Reference time', y = 'Phylogenetic diversity')
      if (unique(outcome$Type) == 'meanPD') out = out + labs(x = 'Reference time', y = 'Mean phylogenetic diversity')
      
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## FD & q-profile ##
  if (class == 'FD' & profile == 'q') {
    
    outcome$Tau = paste('Tau = ', round(outcome$Tau, 3), sep = '')
    out = ggplot(outcome, aes(x = Order.q, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic functional diversity')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'Order q', y = 'Empirical functional diversity')
      
    } else {
      
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Order q', y = 'Functional diversity')
    }
    out = out + facet_grid(.~Tau, scales = "free_y")
  }
  
  ## FD & tau-profile ##
  if (class == 'FD' & profile == 'tau') {
    
    outcome$Order.q = paste('Order q = ', outcome$Order.q, sep = '')
    out = ggplot(outcome, aes(x = Tau, y = qFD, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qFD.LCL, ymax = qFD.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'Tau', y = 'Asymptotic functional diversity')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'Tau', y = 'Empirical functional diversity')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qFD.LCL, ymax = qFD.UCL), linetype = 0, alpha = 0.2)
      
      out = out + labs(x = 'Tau', y = 'Functional diversity')
    }
    out = out + facet_grid(.~Order.q, scales = "free_y")
  }
  
  ## AUC & q-profile ##
  if (class == 'AUC') {
    out = ggplot(outcome, aes(x = Order.q, y = qAUC, colour = Assemblage, fill = Assemblage))
    
    if (length(unique(outcome$Method)) == 1) {
      out = out + geom_line(size = 1.5) + geom_ribbon(aes(ymin = qAUC.LCL, ymax = qAUC.UCL, fill = Assemblage), linetype = 0, alpha = 0.2)
      
      if (unique(outcome$Method) == 'Asymptotic') out = out + labs(x = 'Order q', y = 'Asymptotic Functional diversity (AUC)')
      if (unique(outcome$Method) == 'Empirical') out = out + labs(x = 'Order q', y = 'Empirical Functional diversity (AUC)')
    } else {
      out = out + geom_line(aes(lty = Method), size = 1.5) + 
        geom_ribbon(data = outcome %>% filter(Method=="Asymptotic"), aes(ymin = qAUC.LCL, ymax = qAUC.UCL), linetype = 0, alpha = 0.2)
      
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


