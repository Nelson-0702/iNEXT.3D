# PDInfo -------------------------------------------------------------------
# Summarizes phylogenetic data information.
#
# Function \code{PDInfo} summarizes phylogenetic data statistics for specified/default reference time.
# @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
# See the function \code{\link{iNEXTPD}} for details.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
# @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
# @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
# then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
# the pooled assemblage. Default is \code{NULL}.
# @return Returns a table of phylogenetic data information, including reference sample size (\code{n}), number of sampling units (\code{nT}),
# number of observed species (\code{S.obs}), observed total branch length, i.e., Faith’s PD (\code{PD.obs}), the first two rare species
# frequency counts and their branch length sums, at specified/default reference time specified in the argument \code{reftime}.
# See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for formulas and interpretations.
# @examples
# # Datatype: abundance data
# data(data.abu)
# data <- data.abu$data
# tree <- data.abu$tree
# out <- PDInfo(data = data, datatype = "abundance", tree = tree)
# out
#
# # Datatype: incidence_raw data
# data(data.inc)
# data <- data.inc$data
# tree <- data.inc$tree
# nT <- data.inc$nT
# out <- PDInfo(data = data, nT = nT, datatype = "incidence_raw", tree = tree)
# out
# 
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
PDInfo <- function(data, datatype = "abundance", tree, reftime=NULL, nT){
  
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  
  if(class(data)[1] == "list" & datatype == "incidence_raw"){
    if(class(nT) == 'data.frame') nT = unlist(nT)
    
    data = lapply(data, function(i) data.frame(i))
    data2 = lapply(data, function(i) {
      i$species = rownames(i)
      return(i)
    })
    nT = as.vector(sapply(data, ncol))
    names(nT) = if(is.null(data)) paste0("assemblage", 1:length(data)) else names(data)
    
    data = data2[[1]]
    for(i in 2:length(data2)){
      data = full_join(data, data2[[i]], by = "species")
    }
    data[is.na(data)] = 0
    rownames(data) = data$species
    data = data[, colnames(data)!="species"]
    
  }
  
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(class(data) == "list") {
      mydata = data
      if(is.null(names(mydata))) names(mydata) <- paste0("assemblage",1:length(mydata))
    } else {
      if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      ntmp <- 0
      for(i in 1:length(nT)){
        mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
        ntmp <- ntmp+nT[i]
      }
      if(is.null(names(nT))) {
        names(mydata) <- paste0("assemblage",1:length(nT))
      }else{
        names(mydata) = names(nT)
      }
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  if(datatype=='abundance'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reftime = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,n,S.obs,PD.obs,`f1*`,`f2*`,g1,g2,Reftime)
  }else if (datatype=='incidence_raw'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reftime = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,`nT`,S.obs,PD.obs,`Q1*`,`Q2*`,R1,R2,Reftime)
  }
  
  return(data.frame(infos))
  
}


# iNEXTPD -------------------------------------------------------------------
# Interpolation (rarefaction) and extrapolation of Chao et al.’s (2010) phylogenetic diversity and mean phylogenetic diversity (phylogenetic Hill numbers)
#
# Function \code{iNEXTPD} computes phylogenetic diversity estimates for rarefied samples and extrapolated samples
# along with confidence intervals and related coverage estimates based on Chao et al.’s (2010) phylogenetic
# diversity (PD) and mean phylogenetic diversity (meanPD). See Chao et al. (2010, 2015, 2016) and
# Hsieh and Chao (2017) for pertinent background and methodologies.
# @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).\cr
# Abundance data: a species-by-site matrix/data.frame of species abundances. The row (species) names of
# data must match the species names in the phylogenetic tree and thus cannot be missing.\cr
# Incidence raw data: species-by-site raw incidence matrix/data.frame. When there are N assemblages
# and thus N matrices, users must first merge the N matrices by species identity to obtain a large
# merged incidence matrix, where the rows of the matrix refer to all species presented in the merged
# data. The row (species) names of data must match the species names in the phylogenetic tree and
# thus cannot be missing.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
# @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
# @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
# @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
# then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
# the pooled assemblage. Default is \code{NULL}.
# @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
# and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
# @param endpoint a positive integer specifying the endpoint for the rarefaction and extrapolation range.
# If \code{NULL}, then \code{endpoint} = double of the reference sample size in each assemblage. It is ignored if \code{size} is given.
# @param knots a positive integer specifying the number of equally-spaced knots between 1 and the \code{endpoint}. Default is 40.
# @param size a sequence of positive integers specifying the sample sizes for which PD or meanPD estimates will be calculated.
# If \code{NULL}, then estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
# and \code{knots}.
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
# Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.

# @return a list of three objects: 
# \code{$DataInfo} for summarizing data information; 
# \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics:
#  \item{\code{$size_based}: size-based PD or meanPD estimates along with their confidence intervals
#  (if \code{nboot > 0}) and relevant statistics information.} \cr
#  \item{\code{$coverage_based}: coverage-based diversity estimates along with confidence intervals
#  (if \code{nboot > 0}) and relevant statistics information.}
# and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
# @examples
# \donttest{
# # Datatype: abundance data
# data(data.abu)
# data <- data.abu$data
# tree <- data.abu$tree
# out <- iNEXTPD(data = data, tree = tree, datatype = "abundance", q = c(0, 1, 2), nboot = 30)
# out
#
# # Datatype: incidence_raw data
# data(data.inc)
# data <- data.inc$data
# tree <- data.inc$tree
# nT <- data.inc$nT
# out <- iNEXTPD(data = data, nT = nT, datatype = "incidence_raw", tree = tree, 
# q = c(0, 1, 2))
# out
# }
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
# Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
# \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
iNEXTPD <- function(data, nT, datatype = "abundance", tree, q = c(0,1,2), reftime=NULL, type = 'PD', endpoint = NULL, knots = 40, size = NULL, nboot = 50, conf = 0.95) {
  
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  # qq <- 0:2
  # if(is.na(pmatch(q, qq) == T) == T) stop("invalid order of q, we only compute q = 0, 1 or 2", call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)[1] == "list" & datatype == "incidence_raw"){
    
    data = lapply(data, function(i) data.frame(i))
    data2 = lapply(data, function(i) {
      i$species = rownames(i)
      return(i)
    })
    nT = as.vector(sapply(data, ncol))
    names(nT) = if(is.null(data)) paste0("assemblage", 1:length(data)) else names(data)
    
    data = data2[[1]]
    for(i in 2:length(data2)){
      data = full_join(data, data2[[i]], by = "species")
    }
    data[is.na(data)] = 0
    rownames(data) = data$species
    data = data[, colnames(data)!="species"]
    
  }
  
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  
  if(datatype=="incidence_raw"){
    if(class(data) == "list") {
      mydata = data
      if(is.null(names(mydata))) names(mydata) <- paste0("assemblage",1:length(mydata))
    } else {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      
      if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      ntmp <- 0
      for(i in 1:length(nT)){
        mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
        ntmp <- ntmp+nT[i]
      }
      if(is.null(names(nT))) {
        names(mydata) <- paste0("assemblage",1:length(nT))
      }else{
        names(mydata) = names(nT)
      }
    }
  }else{
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  if(length(knots)!=length(mydata)) knots <- rep(knots,length(mydata))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(mydata, function(x) 2*sum(x))
      }else if(datatype == "incidence_raw"){
        endpoint <- sapply(mydata, function(x) 2*ncol(x))
      }
    }else{
      if(length(endpoint)!=length(mydata)){
        endpoint <- rep(endpoint,length(mydata))
      }
    }
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }
      
      if(endpoint[i] <= ni){
        mi <- floor(seq(1,endpoint[i],length.out = knots[i]))
      }else{
        mi <- floor(c(seq(1,ni,length.out = floor(knots[i]/2)),
                      seq(ni+1,endpoint[i],length.out = knots[i]-floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else{
    if(class(size)=="numeric"|class(size)=="integer"|class(size)=="double"){
      size <- list(size = size)
    }
    if(length(size)!=length(mydata)) size <- lapply(1:length(mydata), function(x) size[[1]])
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }
  
  FUN <- function(e){
    if(class(mydata)=="list"){
      inextPD(datalist = mydata, datatype = datatype, phylotr = mytree, q = q, reft = reftime, m=size,
              cal = type, nboot=nboot, conf = conf, unconditional_var = TRUE)
    }else{
      return(NULL)
    }
  }
  out <- tryCatch(FUN(e), error = function(e){return()})
  
  ## AsyEst table ##
  index <- rbind(AsyPD(data = data, nT = nT, tree = tree, q = c(0,1,2), datatype = ifelse(datatype=='abundance','abundance','incidence_raw'), nboot = 30,conf = 0.95),
                 ObsPD(data = data, nT = nT, tree = tree, q = c(0,1,2), datatype = ifelse(datatype=='abundance','abundance','incidence_raw'), nboot = 30,conf = 0.95))
  index = index[order(index$Assemblage),]
  LCL <- index$qPD.LCL[index$Method=='Asymptotic']
  UCL <- index$qPD.UCL[index$Method=='Asymptotic']
  index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qPD')
  index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
  index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
  index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
  index[,3:4] = index[,4:3]
  colnames(index) <- c("Assemblage", "Phylogenetic Diversity", "Phylogenetic Observed", "Phylogenetic Estimator", "s.e.", "LCL", "UCL")
  info <- DataInfo3D(data, diversity = 'PD', datatype = datatype, nT, PDtree = tree, PDreftime = reftime)
  return( list("PDInfo"=info, "PDiNextEst"=out, "PDAsyEst"=index) )
}

# estimatePD -------------------------------------------------------------------
# Computes phylogenetic diversity for specified values of sample coverage
#
# Function \code{estimatePD} computes Chao et al.’s (2010, 2016) phylogenetic diversity (PD, effective total branch lengths,
# for diversity order q = 0, 1 and 2) and mean phylogenetic diversity (meanPD, phylogenetic Hill
# numbers or the effective number of lineages, q = 0, 1 and 2) at specified values of sample coverage. See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for formulas and interpretations.
# Use the function \code{iNEXTPD} to compute PD or meanPD for specified sample sizes.
# @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data). See the function \code{\link{iNEXTPD}} for details.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
# @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
# @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
# then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
# the pooled assemblage. Default is \code{NULL}.
# @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
# and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
# @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
# @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1). 
# If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes. 
# If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes. 
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
# Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
# @return Returns a table of the computed phylogenetic diversity (PD or meanPD) for specified/default diversity orders \code{q} and reference times
# for the user-specified values of sample coverage. The corresponding sample sizes and sample coverage values are also provided.
# @examples
# # Datatype: abundance data
# data(data.abu)
# data <- data.abu$data
# tree <- data.abu$tree
# out1 <- estimatePD(data = data, tree = tree, datatype = "abundance", base = "size")
# out1
#
# out2 <- estimatePD(data = data, tree = tree, datatype = "abundance", base = "coverage")
# out2
# 
# # Datatype: incidence_raw data
# data(data.inc)
# data <- data.inc$data
# tree <- data.inc$tree
# nT <- data.inc$nT
# out <- estimatePD(data = data, nT = nT, tree = tree, datatype = "incidence_raw", base = "coverage")
# out
# 
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
# Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
# \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
estimatePD <- function(data, nT, tree, datatype = "abundance", q = c(0,1,2), reftime=NULL, type = 'PD', base = "coverage", level = NULL, nboot = 50, conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if(datatype == "incidence_freq") stop("The PD can only accept 'datatype = incidence_raw'.")
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  #if (length(level)>1) stop('Currently, we only accept one fixed level of coverage.')
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(class(data) == "list") {
      mydata = data
      if(is.null(names(mydata))) names(mydata) <- paste0("assemblage",1:length(mydata))
    } else {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      
      if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      ntmp <- 0
      for(i in 1:length(nT)){
        mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
        ntmp <- ntmp+nT[i]
      }
      if(is.null(names(nT))) {
        names(mydata) <- paste0("assemblage",1:length(nT))
      }else{
        names(mydata) = names(nT)
      }
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[pool.name,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  if (is.null(level) & base == 'size') {
    if(datatype == "abundance") {
      level <- sapply(mydata, function(x) 2*sum(x))
    }else if(datatype == "incidence_raw"){
      level <- sapply(mydata, function(x) 2*ncol(x))
    }
    level <- lapply(1:length(mydata), function(i) min(level))
  } else if (is.null(level) & base == 'coverage') {
    if(datatype=='abundance'){
      level <- sapply(mydata,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })
      
    }else if(datatype=='incidence_raw'){
      level <- sapply(mydata,function(x){
        ni <- ncol(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })
    }
    level <- min(level)
  }
  
  if (base == "size") {
    tmp <- inextPD(datalist = mydata, datatype = datatype, phylotr = mytree,q = q, 
                   reft = reftime,m = level, cal = type, nboot=nboot, conf = conf, unconditional_var = FALSE)$size_based %>% 
      select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
  } else if (base == "coverage") {
    tmp <- invChatPD(datalist = mydata, datatype = datatype, phylotr = mytree, q = q,
                     reft = reftime, cal = type, level = level, nboot, conf)
  }
  
  return(tmp)
}


# AsyPD -------------------------------------------------------------------
# Computes asymptotic estimates for phylogenetic diversity and mean phylogenetic diversity (phylogenetic Hill numbers)
#
# Function \code{AsyPD} computes asymptotic phylogenetic diversity estimates with respect to specified/default
# diversity order q and reference time to infer true phylogenetic diversity (PD) or phylogenetic Hill numbers (meanPD). See Chao et al. (2015) and Hsieh and Chao (2017) for the statistical estimation detail.
# @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
# See the function \code{\link{iNEXTPD}} for details.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
# @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
# @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
# @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
# then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
# the pooled assemblage. Default is \code{NULL}.
# @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
# and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
# Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
# @return Returns a table of estimated asymptotic phylogenetic diversity estimates (\code{type = "PD"}) or
# phylogenetic Hill numbers (\code{type = "meanPD"}) with respect to specified/default order \code{q} and
# reference time specified in the argument \code{reftime}.
# @examples
# # Datatype: abundance data
# data(data.abu)
# data <- data.abu$data
# tree <- data.abu$tree
# out <- AsyPD(data = data, datatype = "abundance", tree = tree,
# q = seq(0, 2, by = 0.25), nboot = 30)
# out
#
# # Datatype: incidence_raw data
# data(data.inc)
# data <- data.inc$data
# tree <- data.inc$tree
# nT <- data.inc$nT
# out <- AsyPD(data = data, nT = nT, datatype = "incidence_raw",
# tree = tree, q = seq(0, 2, by = 0.25))
# out
# 
# @references
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
AsyPD <- function(data,nT,datatype = "abundance",tree,q = seq(0,2,by = 0.25),reftime = NULL,type = 'PD',nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
  #if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((conf < 0) | (conf < 0) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(class(data) == "list") {
      mydata = data
      if(is.null(names(mydata))) names(mydata) <- paste0("assemblage",1:length(mydata))
    } else {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      
      if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      ntmp <- 0
      for(i in 1:length(nT)){
        mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
        ntmp <- ntmp+nT[i]
      }
      if(is.null(names(nT))) {
        names(mydata) <- paste0("assemblage",1:length(nT))
      }else{
        names(mydata) = names(nT)
      }
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  FUN = function(e){
    asymPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)# mytree is pooled tree of class phylo
  }
  #out <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}


# ObsPD -------------------------------------------------------------------
# Computes observed phylogenetic diversity and phylogenetic Hill numbers
#
# Function \code{ObsPD} computes empirical or observed phylogenetic diversity (PD) and phylogenetic Hill
# numbers (meanPD, mean phylogenetic diversity) for specified/default order \code{q} and reference
# time specified in the argument \code{reftime}. See Chao et al. (2010) for details of PD and meanPD.
# @param data a matrix/data.frame of species abundances (for abundance data) or species-by-site incidence raw matrix/data.frame (for incidence data).
# See the function \code{\link{iNEXTPD}} for details.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species-by-site raw incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}.
# @param tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage.
# @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
# @param reftime a positive value or sequence specifying the reference times for diversity computation. If \code{NULL},
# then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
# the pooled assemblage. Default is \code{NULL}.
# @param type desired diversity type: \code{type = "PD"} for Chao et al. (2010) phylogenetic diversity
# and \code{type = "meanPD"} for mean phylogenetic diversity (phylogenetic Hill number). Default is \code{"PD"}.
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals.
# Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
# @return Returns a table of empirical (observed) phylogenetic diversity (\code{type = "PD"}) or phylogenetic Hill number (\code{type= "meanPD"})
# for specified/default order q and reference time.
# @examples
# # Datatype: abundance data
# data(data.abu)
# data <- data.abu$data
# tree <- data.abu$tree
# out <- ObsPD(data = data, datatype = "abundance", tree = tree,
# q = seq(0, 2, by = 0.25))
# out
#
# # Datatype: incidence_raw data
# data(data.inc)
# data <- data.inc$data
# tree <- data.inc$tree
# nT <- data.inc$nT
# out <- ObsPD(data = data, nT = nT, datatype = "incidence_raw",
# tree = tree, q = seq(0, 2, by = 0.25))
# out
# 
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
ObsPD <- function(data,nT,datatype = "abundance",tree,q = seq(0, 2, by = 0.25),reftime = NULL,type = "PD",
                   nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if(datatype == "incidence_freq") stop("The diversity = 'PD' can only accept 'datatype = incidence_raw'.")
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  if (sum(q<0)>0) stop("q must be a positive number", call. = FALSE)
  # if ((profile != "q") & (profile != "time")) stop("invalid profile", call. = FALSE)
  # if (length(tprofile_times) == 1 & is.null(tprofile_times)==F) stop("length of time should be greater than one", call. = FALSE)
  # if (sum(tprofile_times<0)>=1 & is.null(tprofile_times)==F) stop("time must be a positive number", call. = FALSE)
  # if (is.null(knots) ==F) {
  #   if ((knots < 0) | (is.numeric(knots)==F) | (knots%%1>0)) {
  #     stop('knot must be a nonnegative integer, We use "knots" = 50 to calculate!', call. = FALSE)
  #   }
  # }
  if(c("numeric") %in% class(data) | c("integer") %in% class(data) | c("double") %in% class(data) ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  
  if(datatype=="incidence_raw"){
    if(class(data) == "list") {
      mydata = data
      if(is.null(names(mydata))) names(mydata) <- paste0("assemblage",1:length(mydata))
    } else {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      
      if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
      ntmp <- 0
      for(i in 1:length(nT)){
        mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
        ntmp <- ntmp+nT[i]
      }
      if(is.null(names(nT))) {
        names(mydata) <- paste0("assemblage",1:length(nT))
      }else{
        names(mydata) = names(nT)
      }
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)
  
  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  
  FUN <- function(e){
    EmpPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)
  }
  #temp <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}


