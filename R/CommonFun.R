# as.incfreq -------------------------------------------------------------------
# Transform incidence raw data to incidence frequencies (iNEXT input format) 
# 
# \code{as.incfreq}: transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units).
# @param x a \code{data.frame} or \code{matirx} of species by sites presence-absence matrix.
# @return a \code{vector} of species incidence frequencies, the first entry of the input data must be total number of sampling units.
# @examples
# data(ciliates)
# lapply(ciliates, as.incfreq)
as.incfreq <- function(x){
  class_x <- class(x)[1]
  if(class_x == "data.frame" | class_x == "matrix"){
    a <- sort(as.numeric(unique(c(unlist(x)))))
    if(!identical(a, c(0,1))){
      warning("Invalid data type, the element of species by sites presence-absence matrix should be 0 or 1. Set nonzero elements as 1.")
      x <- (x > 0)
    }
    nT <- ncol(x)
    y <- rowSums(x)
    y <- c(nT, y)
    # names(y) <- c("nT", rownames(x))
    y
  }else if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    c(1, x) 
  }else{
    stop("Invalid data type, it should be a data.frame or matrix.")
  }
}


# as.abucount -------------------------------------------------------------------
# Transform abundance raw data to abundance row-sum counts (iNEXT input format) 
# 
# \code{as.abucount}: transform abundance raw data (a species by sites matrix) to abundance rwo-sum counts data (iNEXT input format).
# @param x a \code{data.frame} or \code{matirx} of species by sites matrix.
# @return a \code{vector} of species abundance row-sum counts.
# @examples
# data(ciliates)
# lapply(ciliates, as.abucount)
as.abucount <- function(x){
  class_x <- class(x)[1]
  if(class_x == "data.frame" | class_x == "matrix"){
    y <- rowSums(x)
    y
  }else if(class_x=="numeric" | class_x=="integer" | class_x=="double"){
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    x 
  }else{
    stop("invalid data type, it should be a data.frame or matrix.")
  }
}


# EstiBootComm.Ind -------------------------------------------------------------------
# Estimation of species relative abundance distribution
# 
# \code{EstiBootComm.Ind} Estimation of species reletive abundance distribution to obtain bootstrap s.e.
# 
# @param Spec a vector of species abundances
# @return a vector of reltavie abundance
# @seealso \code{\link{EstiBootComm.Sam}}
# @examples 
# data(spider)
# EstiBootComm.Ind(spider$Girdled)
EstiBootComm.Ind <- function(Spec){
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)        #sample size
  f1 <- sum(Spec == 1)   #singleton 
  f2 <- sum(Spec == 2)   #doubleton
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  a <- f1/n*A
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  if(f0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))  	#estimation of relative abundance of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))		  							#Output: a vector of estimated relative abundance
}


# EstiBootComm.Sam -------------------------------------------------------------------
# Estimation of species detection distribution
# 
# \code{EstiBootComm.Sam} Estimation of species detection distribution to obtain bootstrap s.e.
# 
# @param Spec a vector of species incidence, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @return a vector of estimated detection probability
# @seealso \code{\link{EstiBootComm.Sam}}
# @examples 
# data(ant)
# EstiBootComm.Sam(ant$h50m)
EstiBootComm.Sam <- function(Spec){
  nT <- Spec[1]
  Spec <- Spec[-1]
  Sobs <- sum(Spec > 0)   #observed species
  Q1 <- sum(Spec == 1) 	#singleton 
  Q2 <- sum(Spec == 2) 	#doubleton
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  a <- Q1/nT*A
  b <- sum(Spec / nT * (1 - Spec / nT) ^ nT)
  
  if(Q0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  
  Prob.hat <- Spec / nT * (1 - w * (1 - Spec / nT) ^ nT)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))  	#estimation of detection probability of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))									#Output: a vector of estimated detection probability
}


# Chat.Ind -------------------------------------------------------------------
# Abundance-based sample coverage
# 
# \code{Chat.Ind} Estimation of abundance-based sample coverage function
# 
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
# @export
Chat.Ind <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}


# Chat.Sam -------------------------------------------------------------------
# Incidence-based sample coverage
# 
# \code{Chat.Sam} Estimation of incidence-based sample coverage function
# 
# @param x a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param t a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
# @export
Chat.Sam <- function(x, t){
  nT <- x[1]
  y <- x[-1]
  y <- y[y>0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)  #estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  Sub <- function(t){
    if(t < nT) {
      yy <- y[(nT-y)>=t]
      out <- 1 - sum(yy / U * exp(lgamma(nT-yy+1)-lgamma(nT-yy-t+1)-lgamma(nT)+lgamma(nT-t)))     
    }
    #if(t < nT) out <- 1 - sum(y / U * exp(lchoose(nT - y, t) - lchoose(nT - 1, t)))
    if(t == nT) out <- 1 - Q1 / U * A
    if(t > nT) out <- 1 - Q1 / U * A^(t - nT + 1)
    out
  }
  sapply(t, Sub)		
}


# invChat -------------------------------------------------------------------
# Compute species diversity with fixed sample coverage
# 
# \code{invChat} compute species diversity with fixed sample coverage
# @param x a \code{data.frame} or \code{list} for species abundance/incidence frequencies.
# @param q a numerical vector of the order of Hill number.
# @param datatype the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
# @param C a specific sample coverage to compare, which is between 0 to 1. Default is the minimum of double sample size for all sites.
# @param nboot the number of bootstrap times to obtain confidence interval. If confidence interval is not desired, use 0 to skip this time-consuming step.
# @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
# @return a \code{data.frame} with fixed sample coverage to compare species diversity.
# @examples
# data(spider)
# incChat(spider, "abundance")
# incChat(spider, "abundance", 0.85)
# 
# @export

invChat <- function (x, q, datatype = "abundance", C = NULL,nboot=0, conf = NULL) {
  qtile <- qnorm(1-(1-conf)/2)
  TYPE <- c("abundance", "incidence")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (class(x) == "numeric" | class(x) == "integer"){
    x <- list(data = x)
  }
  if (class(x) == "data.frame" | class(x) ==  "matrix"){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if (datatype == "abundance") {
    if (class(x) == "list") {
      if (is.null(C)) {
        C <- min(unlist(lapply(x, function(x) Chat.Ind(x,2*sum(x)))))
      }
      Community = rep(names(x),each = length(q)*length(C))
      out <- lapply(x, function(x_){
        est <- invChat.Ind(x_, q, C)
        if (sum(round(est$m) > 2 * sum(x_))>0) 
          warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
        
        if(nboot>1){
          Prob.hat <- EstiBootComm.Ind(x_)
          Abun.Mat <- rmultinom(nboot, sum(x_), Prob.hat)
          ses <- apply(matrix(apply(Abun.Mat,2 ,function(a) invChat.Ind(a, q,C)$qD),
                              nrow = length(q) * length(C)),1,sd)
        }else{
          ses <- rep(0,nrow(est))
        }
        est <- cbind(est,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
        est
      })
      out <- do.call(rbind,out)
      #out <- do.call(rbind, lapply(x, function(x) invChat.Ind(x, q, C)))
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-4)),(ncol(out)-2),(ncol(out)-1),(ncol(out)-3))]
      rownames(out) <- NULL
      out = out %>% select(c('Assemblage', 'goalSC', 'SC', 'm', 'Method', 'Order.q', 'qD', 'qD.LCL', 'qD.UCL'))
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence") {
    if (class(x) == "list") {
      if (is.null(C)) {
        C <- min(unlist(lapply(x, function(x) Chat.Sam(x,2*max(x)))))
      }
      Community = rep(names(x),each = length(q)*length(C))
      out <- lapply(x, function(x_){
        est <- invChat.Sam(x_, q, C)
        if (sum(round(est$t) > 2 * max(x_))>0) 
          warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
        
        if(nboot>1){
          Prob.hat <- EstiBootComm.Sam(x_)
          Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, x_[1], p)))
          Abun.Mat <- matrix(c(rbind(x_[1], Abun.Mat)),ncol=nboot)
          tmp <- which(colSums(Abun.Mat)==x_[1])
          if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
          if(ncol(Abun.Mat)==0){
            warning("Insufficient data to compute bootstrap s.e.")
          }
          ses <- apply(matrix(apply(Abun.Mat,2 ,function(a) invChat.Sam(a, q,C)$qD),nrow = length(q)* length(C)),1,sd)
        }else{
          ses <- rep(0,nrow(est))
        }
        est <- cbind(est,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
      })
      out <- do.call(rbind,out)
      #out <- do.call(rbind, lapply(x, function(x) invChat.Sam(x,q,C,nboot, conf)))
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-4)),(ncol(out)-2),(ncol(out)-1),(ncol(out)-3))]
      rownames(out) <- NULL
      out = out %>% select(c('Assemblage', 'goalSC', 'SC', 't', 'Method', 'Order.q', 'qD', 'qD.LCL', 'qD.UCL'))
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}


# invSize -------------------------------------------------------------------
invSize <- function(x, q, datatype="abundance", size=NULL, nboot=0, conf=NULL){
  qtile <- qnorm(1-(1-conf)/2)
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  if (class(x) == "numeric" | class(x) == "integer"){
    x <- list(data = x)
  }
  if (class(x) == "data.frame" | class(x) ==  "matrix"){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if(datatype=="abundance"){
    if (class(x) == "list") {
      if (is.null(size)) {
        size <- min(unlist(lapply(x, function(x) 2*sum(x))))
      } 
      Community = rep(names(x),each = length(q)*length(size))
      out <- lapply(x, function(x_){
        est <- invSize.Ind(x_, q, size)
        if(nboot>1){
          Prob.hat <- EstiBootComm.Ind(x_)
          Abun.Mat <- rmultinom(nboot, sum(x_), Prob.hat)
          ses <- apply(matrix(apply(Abun.Mat,2 ,function(a) invSize.Ind(a, q,size)$qD),
                              nrow = length(q)* length(size)),1,sd)
        }else{
          ses <- rep(0,nrow(est))
        }
        est <- cbind(est,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
        est
      })
      out <- do.call(rbind,out)
      #out <- do.call(rbind, lapply(x, function(x) invChat.Ind(x, q, C)))
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence") {
    if (class(x) == "list") {
      if (is.null(size)) {
        size <- min(unlist(lapply(x, function(x) 2*max(x))))
      }
      Community = rep(names(x),each = length(q)*length(size))
      out <- lapply(x, function(x_){
        est <- invSize.Sam(x_, q, size)
        if(nboot>1){
          Prob.hat <- EstiBootComm.Sam(x_)
          Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, x_[1], p)))
          Abun.Mat <- matrix(c(rbind(x_[1], Abun.Mat)),ncol=nboot)
          tmp <- which(colSums(Abun.Mat)==x_[1])
          if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
          if(ncol(Abun.Mat)==0){
            warning("Insufficient data to compute bootstrap s.e.")
          }
          ses <- apply(matrix(apply(Abun.Mat,2 ,function(a) invSize.Sam(a, q,size)$qD),
                              nrow = length(q)* length(size)),1,sd)
        }else{
          ses <- rep(0,nrow(est))
        }
        est <- cbind(est,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
      })
      out <- do.call(rbind,out)
      #out <- do.call(rbind, lapply(x, function(x) invChat.Sam(x,q,C,nboot, conf)))
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}


# invChat.Ind -------------------------------------------------------------------
invChat.Ind <- function (x, q, C) {
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- sum(x)
  refC <- Chat.Ind(x, n)
  f <- function(m, C) abs(Chat.Ind(x, m) - C)
  mm <- sapply(C, function(cvrg){
    if (refC == cvrg) {
      mm <- n
    }else if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = sum(x))
      mm <- opt$minimum
      # mm <- round(mm)
    }else if (refC < cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }else if (f1 == 0 & f2 > 0) {
        A <- 0
      }else if(f1 == 1 & f2 == 0) {
        A <- 0
      }else if(f1 == 0 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A==0,0,(log(n/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
      # mm <- round(mm)
    }
    mm
  })
  mm[mm < 1] <- 1
  SC <- Chat.Ind(x,mm)
  # if (sum(round(mm) > 2 * n)>0) 
  #   warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  
  out <- TD.m.est(x = x,m = mm,qs = q)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(mm,length(q))
  order <- rep(q,each = length(mm))
  SC <- rep(SC,length(q))
  data.frame(m = m,Method = method,Order.q = order,
             SC=SC,qD = out,goalSC = rep(C,length(q)))
  # if (nboot==0|is.null(conf)) {
  #   out <- TD.m.est(x = x,m = mm,qs = q)
  #   out <- subset(iNEXT.Ind(x,q,m = c(1,mm),se = FALSE),m==mm)
  #   out <- out[,c(1,2,3,5,4)]
  # }else {
  #   out <- subset(iNEXT.Ind(x,q,m = c(1,mm),se = TRUE,conf = conf,nboot = nboot), m==mm)
  #   out <- out[,c(1, 2, 3, 7, 4, 5, 6)]
  # }
  # out <- out[!duplicated(out), ]
  # out
}


# invChat.Sam -------------------------------------------------------------------
invChat.Sam <- function (x, q, C) {
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- max(x)
  refC <- Chat.Sam(x, n)
  f <- function(m, C) abs(Chat.Sam(x, m) - C)
  mm <- sapply(C, function(cvrg){
    if (refC == cvrg) {
      mm <- n
    }else if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = max(x))
      mm <- opt$minimum
      # mm <- round(mm)
    }else if (refC < cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if(f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }else if(f1 == 0) {
        A <- 0
      }else if(f1 == 1 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A==0,0,(log(U/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
      # mm <- round(mm)
    }
    mm
  })
  mm[mm < 1] <- 1
  SC <- Chat.Sam(x,mm)
  # if (sum(round(mm) > 2 * n)>0) 
  #   warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  out <- TD.m.est_inc(y = x,t_ = mm,qs = q)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(mm,length(q))
  order <- rep(q,each = length(mm))
  SC <- rep(SC,length(q))
  data.frame(t = m,Method = method,Order.q = order,
             SC=SC,qD = out,goalSC = rep(C,length(q)))
  
  # if (nboot==0|is.null(conf)) {
  #   out <- subset(iNEXT.Sam(x,q,t = c(1,mm),se = FALSE), t==mm)
  #   out <- out[,c(1,2,3,5,4)]
  # }else {
  #   out <- subset(iNEXT.Sam(x,q,t = c(1,mm),se = TRUE,conf = conf,nboot = nboot), t==mm)
  #   out <- out[, c(1, 2, 3, 7, 4, 5, 6)]
  # }
  # out <- out[!duplicated(out), ]
  # out
}


# invSize.Ind -------------------------------------------------------------------
invSize.Ind <- function(x, q, size){
  m <- NULL # no visible binding for global variable 'm'
  
  n <- sum(x)
  if(is.null(size)){
    size <- sum(x)
  }
  out <- TD.m.est(x = x,m = size,qs = q)
  SC <- Chat.Ind(x,size)
  method <- ifelse(size>n,'Extrapolation',ifelse(size<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(size,length(q))
  order <- rep(q,each = length(size))
  SC <- rep(SC,length(q))
  data.frame(m = m,Method = method,Order.q = order,SC=SC,qD = out)
  # if(nboot==0|is.null(conf)){
  #   method <- ifelse(size<sum(x), "interpolated", ifelse(size==sum(x), "observed", "extrapolated"))
  #   out <- subset(iNEXT.Ind(x,q,m = c(1,size),se = FALSE), m==size)
  #   out <- out[,c(1,2,3,5,4)]
  #   # out <- data.frame(m=size, method=method, 
  #   #                   SamCov=round(Chat.Ind(x,size),3),
  #   #                   SpeRic=round(Dqhat.Ind(x,0,size),3),
  #   #                   ShaDiv=round(Dqhat.Ind(x,1,size),3),
  #   #                   SimDiv=round(Dqhat.Ind(x,2,size),3))
  #   # colnames(out) <- c("m", "method", "SC", "q = 0", "q = 1", "q = 2")
  # }else{
  #   out <- subset(iNEXT.Ind(x,q,m = c(1,size),se = TRUE,conf = conf,nboot = nboot), m==size)
  #   out <- out[,c(1, 2, 3, 7, 4, 5, 6)]
  # }
  # out <- out[!duplicated(out), ]
  # out
}


# invSize.Sam -------------------------------------------------------------------
invSize.Sam <- function(x, q, size){
  m <- NULL # no visible binding for global variable 'm'
  
  n <- max(x)
  if(is.null(size)){
    size <- max(x)
  }
  out <- TD.m.est_inc(y = x,t_ = size,qs = q)
  SC <- Chat.Sam(x,size)
  method <- ifelse(size>n,'Extrapolation',ifelse(size<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(size,length(q))
  order <- rep(q,each = length(size))
  SC <- rep(SC,length(q))
  data.frame(t = m,Method = method,Order.q = order,SC=SC,qD = out)
  # if(nboot==0|is.null(conf)){
  #   method <- ifelse(size<max(x), "interpolated", ifelse(size==max(x), "observed", "extrapolated"))
  #   out <- subset(iNEXT.Sam(x,q,t = c(1,size),se = FALSE), t==size)
  #   out <- out[,c(1,2,3,5,4)]
  # }else{
  #   out <- subset(iNEXT.Sam(x,q,t = c(1,size),se = TRUE,conf = conf,nboot = nboot), t==size)
  #   out <- out[, c(1, 2, 3, 7, 4, 5, 6)]
  # }
  # out <- out[!duplicated(out), ]
  # out
}



