# FDInfo -------------------------------------------------------------------
# Exhibit basic data information
# 
# \code{FDInfo}: exhibits basic data information
# 
# @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"}, 
# then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <- FunDdata.abu$dij
# FDInfo(data = data, distM = dij, datatype = "abundance")
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <- FunDdata.inc$dij
# FDInfo(data = data, distM = dij, datatype = "incidence_freq")
# }
FDInfo <- function(data, datatype, distM, threshold = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  
  distM <- distM[order_sp,order_sp]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }else if(datatype=='incidence_freq'){
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  } else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",call. = FALSE)
  }
  
  info <- TDInfo(lapply(dat, function(x) data_transform(x, distM, threshold, datatype)$ai %>% round), datatype)
  info$n =lapply(dat, function(x) sum(x))
  info$SC = lapply(dat, function(x) {
    n = sum(x)
    f1 = sum(x==1)
    f2 = sum(x==2)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2) 
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    Chat <- round(1 - f1/n*A, 4)
  })
  colnames(info)[colnames(info) %in% paste0("f",1:10)] = paste0("a",1:10,"'")
  info$threshold = threshold
  return(info)
}



# iNEXTFD -------------------------------------------------------------------
# Interpolation and extrapolation of functional diversity
#
# \code{iNEXTFD}: the seamless rarefaction and extrapolation sampling curves of functional diversity(FD) for q = 0, 1 and 2.
# See Chao et al. (2019) for pertinent background and methods.
# @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: a S by N matrix/data.frame
# where N is the number of assemblages. The element in i-th row and k-th is the abundance of species i in assemblage k. Please note
# that the rownames of data must be the species names matching the species names in distance matrix and thus can't be empty.\cr
# Type (2) incidence frequency data: the sampling unit is quadrat or transect, the observed species was only recorded as presence(detection)/absence(non-detection)
# data in each sampling unit. Likewise, the rownames of data must be the species names matching the species names in phylogeny tree and thus can't be empty. \cr
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# @param endpoint an positive interger specifying the endpoint for rarefaction and
# extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum reference sample size. It will be ignored if \code{size} is given. \cr
# @param knots a positive integer specifying the number of knots between 1 and the \code{endpoint}. Default is 40.\cr
# @param size a sequence of positive integers specifying the sample sizes for which FD estimates will be calculated. If \code{NULL}, then FD estimates will be
# calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}. \cr
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold} = dmean/2. Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table of FD estimates and sample completeness for interpolated or extrapolated sample sizes along with their confidence intervals (if \code{nboot > 0}). \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- iNEXTFD(data = data, distM = dij, datatype = "abundance", nboot = 0)
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- iNEXTFD(data = data, distM = dij, datatype = "incidence_freq")
# out
# }
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
iNEXTFD <- function(data, distM, datatype = "abundance", q = c(0,1,2), endpoint = NULL, 
                    knots = 40, size = NULL, conf = 0.95, nboot = 50, threshold = NULL, nT = NULL) {
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  
  distM <- distM[order_sp,order_sp]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }else if(datatype=='incidence_freq'){
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  } else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",call. = FALSE)
  }
  
  if(length(knots)!=length(dat)) knots <- rep(knots,length(dat))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(dat, function(x) 2*sum(x))
      }else if(datatype == "incidence_freq"){
        endpoint <- sapply(dat, function(x) 2*x[1])
      }
    }else{
      if(length(endpoint)!=length(dat)){
        endpoint <- rep(endpoint,length(dat))
      }
    }
    size <- lapply(1:length(dat),function(i){
      if(datatype == "abundance") {
        ni <- sum(dat[[i]])
      }else if(datatype == "incidence_freq"){
        ni <- dat[[i]][1]
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
    if(length(size)!=length(dat)) size = lapply(1:length(dat), function(x) size[[1]])
    size <- lapply(1:length(dat),function(i){
      if(datatype == "abundance") ni <- sum(dat[[i]]) else ni <- (dat[[i]])[1]
      
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]])) else mi <- size[[i]]
      unique(mi)
    })
  }
  
  FUN <- function(e){
    if(class(dat)=="list"){
      ## size-based
      temp1 = iNextFD(datalist = dat, dij = distM, q = q, datatype = datatype, tau = threshold,
                      nboot = nboot, conf = conf, m = size)
      temp1$qFD.LCL[temp1$qFD.LCL<0] <- 0;temp1$SC.LCL[temp1$SC.LCL<0] <- 0
      temp1$SC.UCL[temp1$SC.UCL>1] <- 1
      obs = filter(temp1, Method == "Observed")
      obs$goalSC = rep(sapply(1:length(dat), 
                              function(i)Coverage(data = dat[[i]], datatype = datatype,
                                                    m = ifelse(datatype == "incidence_freq", dat[[i]][1], sum(dat[[i]])))),
                       each = length(q))
      obs$m = rep(sapply(1:length(dat),
                         function(i) ifelse(datatype == "incidence_freq", dat[[i]][1], sum(dat[[i]]))),
                  each = length(q))
      obs = obs[!colnames(obs) %in% c("SC.s.e.", "SC.LCL", "SC.UCL")]
      if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nt'
      
      ## coverage-based
      temp2 <- lapply(1:length(dat), function(i) invChatFD(datalist = dat[i], dij = distM, q = q, datatype = datatype,
                                                           level = Coverage(data = dat[[i]], datatype = datatype, m = size[[i]]), 
                                                           nboot = nboot, conf = conf, tau = threshold)) %>% do.call(rbind,.)
      # obs_cov = Coverage(data = dat[[i]], datatype = datatype, 
      #                  m = ifelse(datatype == "incidence_freq", dat[[i]][1], sum(dat[[i]])))
      #temp2$Method[temp2$goalSC == obs_cov] = "Observed"
      temp2$qFD.LCL[temp2$qFD.LCL<0] <- 0
      temp2 = bind_rows(temp2, obs)
      
      if (datatype == 'incidence_freq') colnames(temp2)[colnames(temp2) == 'm'] = 'nt'
      temp1$Type = "FD"
      temp2$Type = "FD"
      ans <- list(size_based = temp1, coverage_based = temp2)
      return(ans)
    }else{
      return(NULL)
    }
  }
  out <- tryCatch(FUN(e), error = function(e){return()})
  
  ## AsyEst table ##
  index <- rbind(AsyFD(data = data, distM = distM, q = c(0,1,2), datatype = datatype, nboot = 30, conf = 0.95, threshold = NULL),
                 ObsFD(data = data, distM = distM, q = c(0,1,2), datatype = datatype, nboot = 30, conf = 0.95, threshold = NULL))
  index <- index %>% arrange(., Assemblage)
  LCL <- index$qFD.LCL[index$Method=='Asymptotic']
  UCL <- index$qFD.UCL[index$Method=='Asymptotic']
  index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qFD')
  index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
  index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
  index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
  index[,3:4] = index[,4:3]
  colnames(index) <- c("Assemblage", "Functional Diversity", "Functional Observed", "Functional Estimator", "s.e.", "LCL", "UCL")
  
  info <- DataInfo3D(data, diversity = 'FD', datatype = datatype, FDdistM = distM, FDtype = 'tau_values', FDtau = threshold, nT = nT)
  info$n = lapply(dat, function(x) sum(x))
  info$SC = lapply(dat, function(x) {
    n = sum(x)
    f1 = sum(x==1)
    f2 = sum(x==2)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2) 
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    Chat <- round(1 - f1/n*A, 4)
  })
  colnames(info)[colnames(info) %in% paste0("f",1:10)] = paste0("f",1:10,"'")
  return(list("FDInfo" = info, "FDiNextEst" = out, "FDAsyEst" = index))
}


# estimateFD -------------------------------------------------------------------
# Compute functional diversity with particular sample coverages
#
# \code{estimateFD}: computes functional diversity(FD) with particular user-specified levels of sample coverages.
# See Chao et al. (2019) for pertinent background and methods.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# See \code{\link{iNEXTFD}} for data details.
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}.
# @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
# @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1). 
# If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes. 
# If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes. 
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out1 <- estimateFD(data = data, distM = dij, datatype = "abundance", base = "size")
# out1
# 
# out2 <- estimateFD(data = data, distM = dij, datatype = "abundance", base = "coverage")
# out2
# 
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- estimateFD(data = data, distM = dij, datatype = "incidence_freq", base = "coverage")
# out
# }
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
estimateFD <- function(data, distM, datatype = "abundance", q = c(0,1,2), base = "coverage", threshold = NULL, level = NULL, nboot = 50, conf = 0.95, nT = NULL) {
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop = FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  
  distM <- distM[order_sp,order_sp]
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }else if(datatype=='incidence_freq'){
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  } else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",call. = FALSE)
  }
  
  if (is.null(level) & base == "size") {
    if(datatype == "abundance") {
      level <- sapply(dat, function(x) 2*sum(x))
    }else if(datatype == "incidence_freq"){
      level <- sapply(dat, function(x) 2*x[1])
    }
    level <- lapply(dat, function(i) min(level))
  } else if (is.null(level) & base == "coverage") {
    if(datatype=='abundance'){
      level <- sapply(dat,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
      
    }else if(datatype == 'incidence_freq'){
      level <- sapply(dat,function(x){
        ni <- x[1]
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
    }
    level <- min(level)
  } 
  
  if (base == "size") {
    out = iNextFD(datalist = dat,dij = distM,q = q,datatype = datatype,tau = threshold,
                   nboot = nboot,conf = conf,m = level) %>% 
      select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
    out$qFD.LCL[out$qFD.LCL<0] <- 0
    # out$SC.LCL[out$SC.LCL<0] <- 0
    # out$SC.UCL[out$SC.UCL>1] <- 1
    if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
  } else if (base == "coverage") {
    out <- invChatFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                     level = level, nboot = nboot, conf = conf, tau = threshold)
    out$qFD.LCL[out$qFD.LCL<0] <- 0
    if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
  }
  return(out)
}


# AsyFD -------------------------------------------------------------------
# Asymptotic functional diversity q profile 
#
# \code{AsyFD}: computes asymptotic functional diversity(FD) under certain threshold.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# See \code{\link{AsyFD}} for data details.
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- AsyFD(data = data, distM = dij, datatype = "abundance", q=seq(0,2,0.5), nboot=10)
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- AsyFD(data = data, distM = dij, datatype = "incidence_freq")
# out
# }
AsyFD <- function(data, distM, datatype = "abundance", q = seq(0, 2, by = 0.25), nboot = 50, conf = 0.95, threshold = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }else if(datatype=='incidence_freq'){
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  } else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",call. = FALSE)
  }
  
  out <- FDtable_est(datalist = dat, dij = distM, q = q, datatype = datatype,
                     nboot = nboot, conf = conf, tau = threshold)
  out
}


# ObsFD -------------------------------------------------------------------
# Empirical functional diversity q profile 
#
#\code{ObsFD}: computes Empirical functional diversity(FD) under certain threshold.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# See \code{\link{ObsFD}} for data details.
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param threshold a sequence between 0 and 1 specifying tau. If \code{NULL}, \code{threshold = } dmean. Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates by area under curve with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- ObsFD(data = data, distM = dij, datatype = "abundance")
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- ObsFD(data = data, distM = dij, datatype = "incidence_freq")
# out
# }
ObsFD <- function(data, distM, datatype = "abundance", q = seq(0, 2, by = 0.25), nboot = 50, conf = 0.95, threshold = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(is.null(threshold)) {
    if(datatype=='abundance') {
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }else if(datatype=='incidence_freq'){
      tmp = sapply(dat, function(x) x) %>% apply(., 1, sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    dmin <- min(distM[distM>0])
    #dmax <- max(distM)
    #threshold <- (dmean+dmin)/2
    threshold <- dmean
  } else if(sum(threshold<0)>0|sum(threshold>1)>0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",call. = FALSE)
  }
  
  out <- FDtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype,
                     nboot = nboot, conf = conf, tau = threshold)
  out
}

# AUCInfo -------------------------------------------------------------------
# Exhibit basic data information
# 
# \code{FDInfo}: exhibits basic data information
# 
# @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"}, 
# then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <- FunDdata.abu$dij
# AUCInfo(data = data, distM = dij, datatype = "abundance")
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <- FunDdata.inc$dij
# AUCInfo(data = data, distM = dij, datatype = "incidence_freq")
# }
AUCInfo <- function(data, datatype, distM, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  
  distM <- distM[order_sp,order_sp]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  threshold = t(sapply(dat, function(i){
    
    #index = i>0
    
    
    if(datatype=='abundance') {
      tmp <- matrix(i/sum(i),ncol =1)
    }else if(datatype=='incidence_freq'){
      tmp <- matrix(i[-1]/sum(i[-1]), ncol = 1)
    }
    dmean <- sum ( (tmp %*% t(tmp) ) * distM)
    distM <- distM[tmp>0,tmp>0]
    dmin <- min(distM[distM>0])
    dmax <- max(distM[distM>0])
    c(dmin, dmean, dmax)
  }))
  
  
  info <- cbind(TDInfo(dat, datatype)[,1:4], threshold)
  colnames(info)[5:7] = c("dmin", "dmean", "dmax")
  return(info)
}

# iNEXTAUC -------------------------------------------------------------------
# Interpolation and extrapolation of functional diversity by area under curve
#
# \code{iNEXTAUC}: the seamless rarefaction and extrapolation sampling curves of functional diversity(FD) by area under curve thorough several thresholds for q = 0, 1 and 2.
# @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: a S by N matrix/data.frame
# where N is the number of assemblages. The element in i-th row and k-th is the abundance of species i in assemblage k. Please note
# that the rownames of data must be the species names matching the species names in distance matrix and thus can't be empty.\cr
# Type (2) incidence frequency data: the sampling unit is quadrat or transect, the observed species was only recorded as presence(detection)/absence(non-detection)
# data in each sampling unit. Likewise, the rownames of data must be the species names matching the species names in phylogeny tree and thus can't be empty. \cr
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# @param endpoint an positive interger specifying the endpoint for rarefaction and
# extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum reference sample size. It will be ignored if \code{size} is given. \cr
# @param knots a positive integer specifying the number of knots between 1 and the \code{endpoint}. Default is 40.\cr
# @param size a sequence of positive integers specifying the sample sizes for which FD estimates will be calculated. If \code{NULL}, then FD estimates will be
# calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}. \cr
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param tau a sequence between 0 and 1 specifying tau for integrating area under curve. If \code{NULL}, \code{tau} = (0, 0.01, 0.02,..., 0.99, 1). Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table of AUC estimates and sample completeness for interpolated or extrapolated sample sizes along with their confidence intervals (if \code{nboot > 0}). \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- iNEXTAUC(data = data[,1], distM = dij, datatype = "abundance", nboot = 0)
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- iNEXTAUC(data = data, distM = dij, datatype = "incidence_freq", nboot = 0)
# out
# }
iNEXTAUC <- function(data, distM, datatype = "abundance", q = c(0,1,2), endpoint = NULL, 
                     knots = 20, size = NULL, conf = 0.95, nboot = 50, nT = NULL) {
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if(length(knots)!=length(dat)) knots <- rep(knots,length(dat))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(dat, function(x) 2*sum(x))
      }else if(datatype == "incidence_freq"){
        endpoint <- sapply(dat, function(x) 2*x[1])
      }
    }else{
      if(length(endpoint)!=length(dat)){
        endpoint <- rep(endpoint,length(dat))
      }
    }
    size <- lapply(1:length(dat),function(i){
      if(datatype == "abundance") {
        ni <- sum(dat[[i]])
      }else if(datatype == "incidence_freq"){
        ni <- dat[[i]][1]
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
    if(length(size)!=length(dat)) lapply(1:length(dat), function(x) size[[1]])
    size <- lapply(1:length(dat),function(i){
      ni <- sum(dat[[i]])
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }
  
  FUN <- function(e){
    if(class(dat) == "list"){
      ## size-based
      temp1 = AUCtable_iNextFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                              tau = NULL, nboot = nboot, conf = conf, m = size)
      temp1$qAUC.LCL[temp1$qAUC.LCL<0] <- 0; temp1$SC.LCL[temp1$SC.LCL<0] <- 0
      temp1$SC.UCL[temp1$SC.UCL>1] <- 1
      if (datatype == 'incidence_freq') colnames(temp1)[colnames(temp1) == 'm'] = 'nt'
      
      ## coverage-based
      temp2 <- lapply(1:length(dat), function(i) AUCtable_invFD(datalist = dat[i], dij = distM, q = q, datatype = datatype,
                                                                level = Coverage(data = dat[[i]], datatype = datatype, m = size[[i]]), 
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
  
  ## AsyEst table ##
  index <- rbind(AsyAUC(data = data, distM = distM, q = c(0,1,2), datatype = datatype, nboot = 20, conf = 0.95),
                 ObsAUC(data = data, distM = distM, q = c(0,1,2), datatype = datatype, nboot = 20, conf = 0.95))
  index = index[order(index$Assemblage),]
  LCL <- index$qAUC.LCL[index$Method=='Asymptotic']
  UCL <- index$qAUC.UCL[index$Method=='Asymptotic']
  index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qAUC')
  index <- cbind(index,se = (UCL - index$Asymptotic)/qnorm(1-(1-conf)/2),LCL,UCL)
  index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
  index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
  index[,3:4] = index[,4:3]
  colnames(index) <- c("Assemblage", "Functional Diversity", "Functional Observed", "Functional Estimator", "s.e.", "LCL", "UCL")
  
  info <- DataInfo3D(data, diversity = 'FD', datatype = datatype, FDdistM = distM, FDtype = 'AUC', nT = nT)
  return( list("AUCInfo" = info, "AUCiNextEst" = out, "AUCAsyEst" = index) )
}


# estimateAUC -------------------------------------------------------------------
# Compute functional diversity by area under curve with particular sample coverages
#
#\code{estimateAUC}: computes functional diversity(FD) by area under curve thorough several thresholds with particular user-specified levels of sample coverages.
# See Chao et al. (2019) for pertinent background and methods.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. Default is \code{c(0,1,2)}. \cr
# @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
# @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1). 
# If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes. 
# If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes. 
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param tau a sequence between 0 and 1 specifying tau for integrating area under curve. If \code{NULL}, \code{tau} = (0, 0.01, 0.02,..., 0.99, 1). Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out1 <- estimateAUC(data = data, distM = dij, datatype = "abundance", nboot = 0, base = "size")
# out1
# 
# out2 <- estimateAUC(data = data, distM = dij, datatype = "abundance", nboot = 0, base = "coverage")
# out2
# 
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- estimateAUC(data = data, distM = dij, datatype = "incidence_freq", nboot = 20, base = "coverage")
# out
# }
# @references
# Chao, A., Chiu C.-H. and Jost, L. (2010). functional diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
# Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of functional diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
# Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive functional diversity among multiple assemblages. Systematic Biology 66, 100-111.
estimateAUC <- function(data, distM, datatype = "abundance", q = c(0,1,2), base = "coverage", level = NULL, nboot = 50, conf = 0.95, tau = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  if (is.null(level) & base == "size") {
    if(datatype == "abundance") {
      level <- sapply(dat, function(x) 2*sum(x))
    }else if(datatype == "incidence_freq"){
      level <- sapply(dat, function(x) 2*x[1])
    }
    level <- lapply(dat, function(i) min(level))
  } else if (is.null(level) & base == "coverage") {
    if(datatype=='abundance'){
      level <- sapply(dat,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
      
    }else if(datatype=='incidence_freq'){
      level <- sapply(dat,function(x){
        ni <- x[1]
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
    }
    level <- min(level)
  } 
  
  if (base == 'size') {
    out = AUCtable_iNextFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                            tau = tau, nboot = nboot, conf = conf, m = level) %>% 
      select(-c('SC.s.e.', 'SC.LCL', 'SC.UCL'))
    out$qAUC.LCL[out$qAUC.LCL<0] <- 0
    # out$SC.LCL[out$SC.LCL<0] <- 0
    # out$SC.UCL[out$SC.UCL>1] <- 1
  } else if (base == 'coverage') {
    out <- AUCtable_invFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                          level = level, nboot = nboot, conf = conf, tau = tau)
    if (datatype == 'incidence_freq') colnames(out)[colnames(out) == 'm'] = 'nt'
  }
  return(out)
}


# AsyAUC -------------------------------------------------------------------
# Asymptotic functional diversity q profile by area under curve
#
# \code{AsyAUC}: computes asymptotic functional diversity(FD) by area under curve thorough several thresholds.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# See \code{\link{AsyAUC}} for data details.
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. \cr
# If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param tau a sequence between 0 and 1 specifying tau for integrating area under curve. If \code{NULL}, \code{tau} = (0, 0.01, 0.02,..., 0.99, 1). Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- AsyAUC(data = data[,2], distM = dij, datatype = "abundance", nboot=0)
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- AsyAUC(data = data, distM = dij, datatype = "incidence_freq", nboot=20)
# out
# }
AsyAUC <- function(data, distM, datatype = "abundance", q = seq(0, 2, by = 0.25), nboot = 50, conf = 0.95, tau = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i){
        i$species = rownames(i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  out <- AUCtable_est(datalist = dat, dij = distM, q = q, datatype = datatype,
                     nboot = nboot, conf = conf, tau = tau)
  out
  
}


# ObsAUC -------------------------------------------------------------------
# Empirical functional diversity q profile 
#
# \code{ObsAUC}: computes Empirical functional diversity(FD) by area under curve thorough several thresholds.
# @param data a matrix/data.frame of species abundances/incidences data.\cr
# See \code{\link{ObsAUC}} for data details.
# @param distM a pair wise distance matrix for all pairs of observed species in the pooled assemblage.\cr
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
# or species by sampling-units incidence frequencies (\code{datatype = "incidence_freq"}), default is \code{"abundance"}. \cr
# @param q a sequence of nonnegative integers specifying the diversity orders of FD. \cr
# If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each site to its double sample size. Default is \code{NULL}.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
# @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
# in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
# @param tau a sequence between 0 and 1 specifying tau for integrating area under curve. If \code{NULL}, \code{tau} = (0, 0.01, 0.02,..., 0.99, 1). Default is \code{NULL}.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table including the sample size, sample coverage,
# method (Interpolated or Extrapolated), and diversity estimates by area under curve with each \code{q} for the user-specified sample coverages. \cr\cr
# @examples
# \donttest{
# # Type (1) abundance data (treat incidence frequencies as abundances to save computation time.)
# data(FunDdata.abu)
# data <- FunDdata.abu$data
# dij <-  FunDdata.abu$dij
# out <- ObsAUC(data = data, distM = dij, datatype = "abundance", nboot=20)
# out
# # Type (2) incidence frequency data 
# data(FunDdata.inc)
# data <- FunDdata.inc$data
# dij <-  FunDdata.inc$dij
# out <- ObsAUC(data = data, distM = dij, datatype = "incidence_freq", nboot=20)
# out
# }
ObsAUC <- function(data, distM, datatype = "abundance", q = seq(0, 2, by = 0.25), nboot = 50, conf = 0.95, tau = NULL, nT = NULL){
  distM = as.matrix(distM)
  
  if (datatype == "incidence_raw") {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if(class(data) == "list"){
    if(length(data) == 1){
      data = data[[1]]
    }else{
      region_names = if(is.null(names(data))) paste0("region_", 1:length(data)) else names(data)
      
      data2 = lapply(data, function(i) {
        i = as.matrix(i)
        i = data.frame('species' = rownames(i), i)
        return(i)
      })
      data = data2[[1]]
      for(i in 2:length(data2)){
        data = data.frame(full_join(data, data2[[i]], by = "species"))
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data[!colnames(data) == "species"]
      names(data) = region_names
    }
  }
  
  DATATYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- matrix(data, ncol = 1)
  
  if(datatype=='incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  distM <- distM[rowSums(data)>0,rowSums(data)>0]
  data <- data[rowSums(data)>0,,drop=FALSE]
  if(nrow(data)!=nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  if(is.null(rownames(data))|is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <-  paste0('sp',1:nrow(data))
  }else{
    if(sum(rownames(data) %in% rownames(distM))!=nrow(distM))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  order_sp <- match(rownames(data),rownames(distM))
  distM <- distM[order_sp,order_sp]
  
  if(datatype=='incidence_freq'){
    data <- rbind(nT,data)
  }
  name_sp <- rownames(data)
  dat <- lapply(1:ncol(data), function(k)  {x <- data[,k];names(x) <- name_sp;x})
  if(is.null(colnames(data))) {
    names(dat) <- paste0("site",1:length(dat))
  }else{
    names(dat) = colnames(data)
  }
  
  out <- AUCtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype,
                     nboot = nboot, conf = conf, tau = tau)
  out
  
}


