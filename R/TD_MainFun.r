# DataInfo -------------------------------------------------------------------
# Exhibit basic data information
# 
# \code{DataInfo}: exhibits basic data information
# 
# @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"}, 
# then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
# @examples 
# data(spider)
# DataInfo(spider, datatype="abundance")
DataInfo <- function(x, datatype="abundance", nT){
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  # if(datatype=="incidence_raw"){
  #   if(class(x)=="list"){
  #     x <- lapply(x, as.incfreq)
  #   }else{
  #     x <- as.incfreq(x)
  #   }
  #   datatype <- "incidence"
  # }
  
  if (datatype == "incidence_raw") {
    if (class(data) == "data.frame" | class(data) == "matrix") {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      mydata = list()
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
      data = lapply(mydata, function(i){
        out = as.incfreq(i)
        return(out)
      })
    } else if (class(data) == "list") 
      data <- lapply(data, as.incfreq)
    datatype <- "incidence"
  }
  
  Fun.abun <- function(x){
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x==k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x>0)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)  #estimation of unseen species via Chao1
    A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
    Chat <- round(1 - f1/n*A, 4)
    c(n, Sobs, Chat, fk)
  }
  
  Fun.ince <- function(x){
    nT <- x[1]
    x <- x[-1]
    U <- sum(x)
    Qk <- sapply(1:10, function(k) sum(x==k))
    Q1 <- Qk[1]
    Q2 <- Qk[2]
    Sobs <- sum(x>0)
    Q0.hat <- ifelse(Q2==0, (nT-1)/nT*Q1*(Q1-1)/2, (nT-1)/nT*Q1^2/2/Q2)  #estimation of unseen species via Chao2
    A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
    Chat <- round(1 - Q1/U*A,4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  
  if(datatype == "abundance"){
    if(class(x) == "numeric" | class(x) == "integer"){
      out <- matrix(Fun.abun(x), nrow=1)
    }else if(class(x) == "list"){
      out <- do.call("rbind", lapply(x, Fun.abun))
    } else if(class(x)[1] == "matrix" | class(x) == "data.frame"){
      out <- t(apply(as.matrix(x), 2, Fun.abun))  
    }
    if(nrow(out) > 1){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
    }
    as.data.frame(out)
  }else if(datatype == "incidence"){
    if(class(x) == "numeric" | class(x) == "integer"){
      out <- matrix(Fun.ince(x), nrow=1)
    }else if(class(x) == "list"){
      out <- do.call("rbind", lapply(x, Fun.ince))
    } else if(class(x)[1] == "matrix" | class(x) == "data.frame"){
      out <- t(apply(as.matrix(x), 2, Fun.ince))  
    }
    if(nrow(out) > 1){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
    }
    as.data.frame(out)
  }
}


# iNEXTTD -------------------------------------------------------------------
# iNterpolation and EXTrapolation of Hill number
# 
# \code{iNEXTTD}: Interpolation and extrapolation of Hill number with order q
# 
# @param data a matrix, data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
# @param q a numerical vector of the order of Hill number.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param rowsum a logical variable to check if the input object is raw data (species by sites matrix, \code{rowsum=FALSE}) or iNEXTTD default input (abundance counts or incidence frequencies, \code{rowsum=TRUE}).
# @param size an integer vector of sample sizes (number of individuals or sampling units) for which diversity estimates will be computed. 
# If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{endpoint} and \code{knots} .
# @param endpoint an integer specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation. 
# If NULL, then \code{endpoint} \code{=} double reference sample size.
# @param knots an integer specifying the number of equally-spaced \code{knots} (say K, default is 40) between size 1 and the \code{endpoint};
# each knot represents a particular sample size for which diversity estimate will be calculated.  
# If the \code{endpoint} is smaller than the reference sample size, then \code{iNEXTTD()} computes only the rarefaction esimates for approximately K evenly spaced \code{knots}. 
# If the \code{endpoint} is larger than the reference sample size, then \code{iNEXTTD()} computes rarefaction estimates for approximately K/2 evenly spaced \code{knots} between sample size 1 and the reference sample size, and computes extrapolation estimates for approximately K/2 evenly spaced \code{knots} between the reference sample size and the \code{endpoint}.
# @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
# @param nboot an integer specifying the number of replications.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @importFrom reshape2 dcast
# @return a list of three objects: 
# \code{$DataInfo} for summarizing data information; 
# \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics:
# \item{\code{$size_based}: size-based PD or meanPD estimates along with their confidence intervals
# (if \code{nboot > 0}) and relevant statistics information.} \cr
# \item{\code{$coverage_based}: coverage-based diversity estimates along with confidence intervals
# (if \code{nboot > 0}) and relevant statistics information.}
# and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
# @examples
# ## example for abundance based data (list of vector)
# data(spider)
# out1 <- iNEXTTD(spider, q=c(0,1,2), datatype="abundance")
# out1$DataInfo # showing basic data information.
# out1$AsyEst # showing asymptotic diversity estimates.
# out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
# ## example for abundance based data (data.frame)
# data(bird)
# out2 <- iNEXTTD(bird, q=0, datatype="abundance")
# ggiNEXT(out2)
# \dontrun{
# ## example for incidence frequencies based data (list of data.frame)
# data(ant)
# t <- round(seq(10, 500, length.out=20))
# out3 <- iNEXTTD(ant$h500m, q=1, datatype="incidence_freq", size=t, se=FALSE)
# out3$iNextEst
# }
# 
iNEXTTD <- function(data, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50, nT)
{
  if(datatype == "incidence") stop('Please try datatype = "incidence_freq" or datatype = "incidence_raw".')  
  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]
  
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  # if (datatype == "incidence_raw") {
  #   if (class(data) == "data.frame" | class(data) == "matrix") 
  #     data <- as.incfreq(data)
  #   else if (class(data) == "list") 
  #     data <- lapply(data, as.incfreq)
  #   datatype <- "incidence"
  # }
  
  if (datatype == "incidence_raw") {
    if (class(data) == "data.frame" | class(data) == "matrix") {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      mydata = list()
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
      data = lapply(mydata, function(i){
        out = as.incfreq(i)
        return(out)
      })
    } else if (class(data) == "list") 
      data <- lapply(data, as.incfreq)
    datatype <- "incidence"
  }
  
  Fun <- function(x, q, assem_name){
    x <- as.numeric(unlist(x))
    unconditional_var <- TRUE
    if(datatype == "abundance"){
      if(sum(x)==0) stop("Zero abundance counts in one or more sample sites")
      out <- iNEXT.Ind(Spec=x, q=q, m=size, endpoint=ifelse(is.null(endpoint), 2*sum(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf,unconditional_var)
    }
    if(datatype == "incidence"){
      t <- x[1]
      y <- x[-1]
      if(t>sum(y)){
        warning("Insufficient data to provide reliable estimators and associated s.e.") 
      }
      if(sum(x)==0) stop("Zero incidence frequencies in one or more sample sites")
      
      out <- iNEXT.Sam(Spec=x, q=q, t=size, endpoint=ifelse(is.null(endpoint), 2*max(x), endpoint), knots=knots, se=se, nboot=nboot, conf=conf)  
    }
    if(unconditional_var){
      out <- lapply(out, function(out_) cbind(Assemblage = assem_name, out_))
    }else{
      out[[1]] <- cbind(Assemblage = assem_name, out[[1]])
    }
    out
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  
  z <- qnorm(1-(1-0.95)/2)
  if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    out <- Fun(data,q,'Assemblage1')
    
    out <- list(size_based = out[[1]],
                coverage_based = out[[2]])
    
    index <- rbind(AsyTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95),
                   ObsTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95))
    index = index[order(index$Assemblage),]
    LCL <- index$qD.LCL[index$Method=='Asymptotic']
    UCL <- index$qD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
  
    
  }else if(class_x=="matrix" | class_x=="data.frame"){
    if(is.null(colnames(data))){
      colnames(data) <- sapply(1:ncol(data), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:ncol(data), function(i) {
      tmp <- Fun(data[,i],q,colnames(data)[i])
      tmp
    })
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- rbind(AsyTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95),
                   ObsTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95))
    index = index[order(index$Assemblage),]
    LCL <- index$qD.LCL[index$Method=='Asymptotic']
    UCL <- index$qD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    
    
  }else if(class_x=="list"){
    if(is.null(names(data))){
      names(data) <- sapply(1:length(data), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:length(data), function(i) {
      tmp <- Fun(data[[i]],q,names(data)[i])
      tmp
    })
    
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- rbind(AsyTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95),
                   ObsTD(data = data,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95))
    index = index[order(index$Assemblage),]
    LCL <- index$qD.LCL[index$Method=='Asymptotic']
    UCL <- index$qD.UCL[index$Method=='Asymptotic']
    index <- dcast(index,formula = Assemblage+Order.q~Method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Asymptotic)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    index[,3:4] = index[,4:3]
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
  }else{
    stop("invalid class of data, data should be a object of numeric, matrix, data.frame, or list")
  }
  out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
  out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
  info <- DataInfo(data, datatype)
  
  
  z <- list("DataInfo"=info, "iNextEst"=out, "AsyEst"=index)
  class(z) <- c("iNEXT")
  return(z)
  }


# estimateTD -------------------------------------------------------------------
# Compute species diversity with a particular of sample size/coverage 
# 
# \code{estimateTD}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
# @param data a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
# If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, followed 
# by species incidence frequencies in each column or list.
# @param q a numerical vector of the order of Hill number.
# @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
# @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
# @param nboot the number of bootstrap times to obtain confidence interval. If confidence interval is not desired, use 0 to skip this time-consuming step.
# @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1). 
# If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes. 
# If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes. 
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a \code{data.frame} of species diversity table including the sample size, sample coverage,
# method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
# @examples
# \dontrun{
# data(spider)
# out1 <- estimateTD(spider, q = c(0,1,2), datatype = "abundance", base="size")
# out1
# out2 <- estimateTD(spider, q = c(0,1,2), datatype = "abundance", base="coverage")
# out2
# }
# data(ant)
# out <- estimateTD(ant, q = c(0,1,2), "incidence_freq", base="coverage", level=0.985, conf=0.95)
# out
estimateTD <- function (data, q = c(0,1,2), datatype = "abundance", base = "coverage", level = NULL, nboot=50,
                       conf = 0.95, nT) 
{
  if(datatype == "incidence") stop('Please try datatype = "incidence_freq" or datatype = "incidence_raw".')  
  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]
  
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  # if (datatype == "incidence_raw") {
  #   if (class(data) == "data.frame" | class(data) == "matrix") 
  #     data <- as.incfreq(data)
  #   else if (class(data) == "list") 
  #     data <- lapply(data, as.incfreq)
  #   datatype <- "incidence"
  # }
  
  if (datatype == "incidence_raw") {
    if (class(data) == "data.frame" | class(data) == "matrix") {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      mydata = list()
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
      data = lapply(mydata, function(i){
        out = as.incfreq(i)
        return(out)
      })
    } else if (class(data) == "list") 
      data <- lapply(data, as.incfreq)
    datatype <- "incidence"
  }
  
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  if (base == "size") {
    tmp <- invSize(data, q, datatype, size = level, nboot, conf = conf)
  }
  else if (base == "coverage") {
    tmp <- invChat(data, q, datatype, C = level, nboot, conf = conf)
  }
  tmp$qD.LCL[tmp$qD.LCL<0] <- 0
  tmp
}


# AsyTD -------------------------------------------------------------------
# Asymptotic diversity q profile 
# 
# \code{AsyTD} The estimated and empirical diversity of order q 
# 
# @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
# @param q a nonnegative value or sequence specifying the diversity order. Default is seq(0, 2, by = 0.2).
# @param datatype  data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# or species-by-site incidence frequencies data (\code{datatype = "incidence_freq"}). Default is "abundance".
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table of diversity q profile by 'Estimated' and 'Empirical'
# 
# @examples
# ## example for abundance based data (list of vector)
# # abundance data
# data(spider)
# out1 <- AsyTD(spider, datatype = "abundance")
# out1
# 
# ## example for incidence frequencies based data (list of data.frame)
# data(ant)
# out2 <- AsyTD(ant, datatype = "incidence_freq")
# out2
# 
# 
# @references
# Chao,A. and Jost,L.(2015).Estimating diversity and entropy profiles via discovery rates of new species.
AsyTD <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT){
  if(datatype == "incidence") stop('Please try datatype = "incidence_freq" or datatype = "incidence_raw".')  
  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]
  
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(data) == "data.frame" | class(data) == "matrix") {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      mydata = list()
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
      data = lapply(mydata, function(i){
        out = as.incfreq(i)
        return(out)
      })
    } else if (class(data) == "list") 
      data <- lapply(data, as.incfreq)
    datatype <- "incidence"
  }
  
  if (class(data) == "data.frame" | class(data) ==  "matrix"){
    datalist <- lapply(1:ncol(data), function(i) data[,i])
    if(is.null(colnames(data))) names(datalist) <-  paste0("data",1:ncol(data)) else names(datalist) <- colnames(data)
    data <- datalist
  } else if (class(data) == "numeric" | class(data) == "integer" | class(data) == "double") {
    data <- list(data = data)
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  if(nboot < 0 | round(nboot) - nboot != 0)
    stop("Please enter non-negative integer for nboot.")
  if(conf < 0 | conf > 1)
    stop("Please enter value between zero and one for confident interval.")
  
  
  if(datatype=="abundance"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile(data[[i]],q)
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(data[[i]])
        Abun.Mat <- rmultinom(nboot, sum(data[[i]]), Prob.hat)
        
        error <- qnorm(1-(1-conf)/2) * 
          apply(apply(Abun.Mat, 2, function(xb) Diversity_profile(xb, q)), 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i], "Method" = "Asymptotic")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile.inc(data[[i]],q)
      if(nboot > 1){
        nT <- data[[i]][1]
        Prob.hat <- EstiBootComm.Sam(data[[i]])
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          error <- qnorm(1-(1-conf)/2) * 
            apply(apply(Abun.Mat, 2, function(yb) Diversity_profile.inc(yb, q)), 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i],"Method" = "Asymptotic")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }
  return(out)
}


# ObsTD -------------------------------------------------------------------
# Empirical diversity q profile 
# 
# \code{ObsTD} The estimated and empirical diversity of order q 
# 
# @param data a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
# @param q a nonnegative value or sequence specifying the diversity order. Default is seq(0, 2, by = 0.2).
# @param datatype  data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
# or species-by-site incidence frequencies data (\code{datatype = "incidence_freq"}). Default is "abundance".
# @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
# @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
# @param nT needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the number of sampling units in each assemblage.
# If \code{names(nT) = NULL}, then assemblage are automatically named as "assemblage1", "assemblage2",..., etc. Ignored if \code{datatype = "abundance"}.
# @return a table of diversity q profile by 'Estimated' and 'Empirical'
# 
# @examples
# ## example for abundance based data (list of vector)
# # abundance data
# data(spider)
# out1 <- ObsTD(spider, datatype = "abundance")
# out1
# 
# ## example for incidence frequencies based data (list of data.frame)
# data(ant)
# out2 <- ObsTD(ant, datatype = "incidence_freq")
# out2
ObsTD <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT){
  if(datatype == "incidence") stop('Please try datatype = "incidence_freq" or datatype = "incidence_raw".')  
  TYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]
  
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  # if (datatype == "incidence_raw") {
  #   if (class(data) == "data.frame" | class(data) == "matrix") 
  #     data <- as.incfreq(data)
  #   else if (class(data) == "list") 
  #     data <- lapply(data, as.incfreq)
  #   datatype <- "incidence"
  # }
  
  if (datatype == "incidence_raw") {
    if (class(data) == "data.frame" | class(data) == "matrix") {
      if(class(nT) == 'data.frame') nT = unlist(nT)
      mydata = list()
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
      data = lapply(mydata, function(i){
        out = as.incfreq(i)
        return(out)
      })
    } else if (class(data) == "list") 
      data <- lapply(data, as.incfreq)
    datatype <- "incidence"
  }
  
  if (class(data) == "data.frame" | class(data) ==  "matrix"){
    datalist <- lapply(1:ncol(data), function(i) data[,i])
    if(is.null(colnames(data))) names(datalist) <-  paste0("data",1:ncol(data)) else names(datalist) <- colnames(data)
    data <- datalist
  } else if (class(data) == "numeric" | class(data) == "integer" | class(data) == "double") {
    data <- list(data = data)
  }
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  if(nboot < 0 | round(nboot) - nboot != 0)
    stop("Please enter non-negative integer for nboot.")
  if(conf < 0 | conf > 1)
    stop("Please enter value between zero and one for confident interval.")
  
  if(datatype=="abundance"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile_MLE(data[[i]],q)
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(data[[i]])
        Abun.Mat <- rmultinom(nboot, sum(data[[i]]), Prob.hat)
        
        error <- qnorm(1-(1-conf)/2) * 
          apply(apply(Abun.Mat, 2, function(xb) Diversity_profile_MLE(xb,q)), 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i], "Method" = "Empirical")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile_MLE.inc(data[[i]],q)
      if(nboot > 1){
        nT <- data[[i]][1]
        Prob.hat <- EstiBootComm.Sam(data[[i]])
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          error <- qnorm(1-(1-conf)/2) * 
            apply(apply(Abun.Mat, 2, function(yb) Diversity_profile_MLE.inc(yb,q)), 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i],"Method" = "Empirical")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }
  return(out)
}


