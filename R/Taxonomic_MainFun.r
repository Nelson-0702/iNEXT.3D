# DataInfo -------------------------------------------------------------------
#' Exhibit basic data information
#' 
#' \code{DataInfo}: exhibits basic data information
#' 
#' @param x a vector/matrix/list of species abundances or incidence frequencies.\cr If \code{datatype = "incidence"}, 
#' then the first entry of the input data must be total number of sampling units, followed by species incidence frequencies.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' sampling-unit-based incidence frequencies data (\code{datatype = "incidence_freq"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @return a data.frame of basic data information including sample size, observed species richness, sample coverage estimate, and the first ten abundance/incidence frequency counts.
#' @examples 
#' data(spider)
#' DataInfo(spider, datatype="abundance")
#' @export
DataInfo <- function(x, datatype="abundance"){
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(class(x)=="list"){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
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


# iNEXT -------------------------------------------------------------------
#' iNterpolation and EXTrapolation of Hill number
#' 
#' \code{iNEXT}: Interpolation and extrapolation of Hill number with order q
#' 
#' @param x a matrix, data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list. 
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
#' @param se a logical variable to calculate the bootstrap standard error and \code{conf} confidence interval.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param nboot an integer specifying the number of replications.
#' @importFrom reshape2 dcast
#' @return a list of three objects: \code{$DataInfo} for summarizing data information; 
#' \code{$iNextEst} for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#' and \code{$AsyEst} for showing asymptotic diversity estimates along with related statistics.  
#' @examples
#' ## example for abundance based data (list of vector)
#' data(spider)
#' out1 <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' out1$DataInfo # showing basic data information.
#' out1$AsyEst # showing asymptotic diversity estimates.
#' out1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' ## example for abundance based data (data.frame)
#' data(bird)
#' out2 <- iNEXT(bird, q=0, datatype="abundance")
#' ggiNEXT(out2)
#' \dontrun{
#' ## example for incidence frequencies based data (list of data.frame)
#' data(ant)
#' t <- round(seq(10, 500, length.out=20))
#' out3 <- iNEXT(ant$h500m, q=1, datatype="incidence_freq", size=t, se=FALSE)
#' out3$iNextEst
#' }
#' @export
#' 
iNEXT <- function(x, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(x)[1]
  
  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported after v2.0.8, 
         please try datatype="incidence_freq".')  
  }
  if(datatype=="incidence_freq") datatype <- "incidence"
  
  if(datatype=="incidence_raw"){
    if(class_x=="list"){
      x <- lapply(x, as.incfreq)
    }else{
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  #   if(rowsum==FALSE){
  #     if(datatype=="abundance"){
  #       if(class_x =="list"){
  #         x <- lapply(x, as.abucount)
  #       }else{
  #         x <- as.abucount(x)
  #       }
  #     }else if(datatype=="incidence"){
  #       if(class_x =="list"){
  #         x <- lapply(x, as.incfreq)
  #       }else{
  #         x <- as.incfreq(x)
  #       }
  #     }
  #   }
  
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
    out <- Fun(x,q,'Assemblage1')
    
    out <- list(size_based = out[[1]],
                coverage_based = out[[2]])
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq')
                     ,nboot = 100,conf = 0.95)
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Order.q~method,value.var = 'qD')
    index <- cbind(index[,-1],se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    colnames(index) <- c("Observed","Estimator","Est_s.e.","95% Lower","95% Upper")
    # index <- rbind(as.matrix(ChaoSpecies(x, datatype, conf)), 
    #                as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)),
    #                as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))
    rownames(index) <- c("Species Richness", "Shannon diversity", "Simpson diversity")
  
    
  }else if(class_x=="matrix" | class_x=="data.frame"){
    if(is.null(colnames(x))){
      colnames(x) <- sapply(1:ncol(x), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:ncol(x), function(i) {
      tmp <- Fun(x[,i],q,colnames(x)[i])
      tmp
    })
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95)
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Assemblage+Order.q~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, ncol(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    
    
  }else if(class_x=="list"){
    if(is.null(names(x))){
      names(x) <- sapply(1:length(x), function(i) paste0('assemblage',i))
    }
    out <- lapply(1:length(x), function(i) {
      tmp <- Fun(x[[i]],q,names(x)[i])
      tmp
    })
    
    out <- list(size_based = do.call(rbind,lapply(out,  function(out_){out_[[1]]})),
                coverage_based = do.call(rbind,lapply(out,  function(out_){out_[[2]]})))
    
    index <- AsyD(x = x,q = c(0,1,2),datatype = ifelse(datatype=='abundance','abundance','incidence_freq'),nboot = 100,conf = 0.95)
    LCL <- index$qD.LCL[index$method=='Estimated']
    UCL <- index$qD.UCL[index$method=='Estimated']
    index <- dcast(index,formula = Assemblage+Order.q~method,value.var = 'qD')
    index <- cbind(index,se = (UCL - index$Estimated)/z,LCL,UCL)
    index$LCL[index$LCL<index$Empirical & index$Order.q==0] <- index$Empirical[index$LCL<index$Empirical & index$Order.q==0]
    index$Order.q <- c('Species richness','Shannon diversity','Simpson diversity')
    # arr <- array(0, dim = c(3, 5, length(x)))
    # arr[1,,] <- t(as.matrix(ChaoSpecies(x, datatype, conf)))
    # arr[2,,] <- t(as.matrix(ChaoEntropy(x, datatype, transform=TRUE, conf)))
    # arr[3,,] <- t(as.matrix(EstSimpson(x, datatype, transform=TRUE, conf)))  
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity", "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.", "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3,1))
    # index <- dcast(as.data.frame(index), formula = Var1+Var2~Var3, value.var = "Freq")
    colnames(index) <- c("Assemblage", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
  }else{
    stop("invalid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  out$size_based$Assemblage <- as.character(out$size_based$Assemblage)
  out$coverage_based$Assemblage <- as.character(out$coverage_based$Assemblage)
  info <- DataInfo(x, datatype)
  
  
  z <- list("DataInfo"=info, "iNextEst"=out, "AsyEst"=index)
  class(z) <- c("iNEXT")
  return(z)
  }


# estimateD -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage 
#' 
#' \code{estimateD}: computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#' @param x a \code{data.frame} or \code{list} of species abundances or incidence frequencies.\cr 
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
#' @return a \code{data.frame} of species diversity table including the sample size, sample coverage,
#' method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#' @examples
#' \dontrun{
#' data(spider)
#' out1 <- estimateD(spider, q = c(0,1,2), datatype = "abundance", base="size")
#' out1
#' out2 <- estimateD(spider, q = c(0,1,2), datatype = "abundance", base="coverage")
#' out2
#' }
#' data(ant)
#' out <- estimateD(ant, q = c(0,1,2), "incidence_freq", base="coverage", level=0.985, conf=0.95)
#' out
#' @export
estimateD <- function (x, q = c(0,1,2), datatype = "abundance", base = "size", level = NULL, nboot=50,
                       conf = 0.95) 
{
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence") {
    stop("datatype=\"incidence\" was no longer supported after v2.0.8, \n         please try datatype=\"incidence_freq\".")
  }
  if (datatype == "incidence_freq") 
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "data.frame" | class(x) == "matrix") 
      x <- as.incfreq(x)
    else if (class(x) == "list") 
      x <- lapply(x, as.incfreq)
    datatype <- "incidence"
  }
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) 
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1) 
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  if (base == "size") {
    tmp <- invSize(x, q, datatype, size = level, nboot, conf = conf)
  }
  else if (base == "coverage") {
    tmp <- invChat(x, q, datatype, C = level, nboot, conf = conf)
  }
  tmp$qD.LCL[tmp$qD.LCL<0] <- 0
  tmp
}


# AsyD -------------------------------------------------------------------
#' Asymptotic diversity q profile 
#' 
#' \code{AsyD} The estimated and empirical diversity of order q 
#' 
#' @param x a matrix/data.frame (species by sites), or list of species abundances or incidence frequencies. If \code{datatype = "incidence_freq"}, then the first entry of the input data must be total number of sampling units in each column or list.
#' @param q a nonnegative value or sequence specifying the diversity order. Default is seq(0, 2, by = 0.2).
#' @param datatype  data type of input data: individual-based abundance data (\code{datatype = "abundance"}),  
#' or species-by-site incidence frequencies data (\code{datatype = "incidence_freq"}). Default is "abundance".
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @return a table of diversity q profile by 'Estimated' and 'Empirical'
#' 
#' @examples
#' ## example for abundance based data (list of vector)
#' # abundance data
#' data(spider)
#' out1 <- AsyD(spider, datatype = "abundance")
#' out1
#' 
#' ## example for incidence frequencies based data (list of data.frame)
#' data(ant)
#' out2 <- AsyD(ant, datatype = "incidence_freq", nboot = 0)
#' out2
#' 
#' 
#' @references
#' Chao,A. and Jost,L.(2015).Estimating diversity and entropy profiles via discovery rates of new species.
#' @export
AsyD <- function(x, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95){
  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported, 
         please try datatype = "incidence_freq".')  
  }
  if(datatype == "incidence_raw"){
    stop('datatype="incidence_raw" was no longer supported, 
         please try datatype = "incidence_freq".')  
  }
  TYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
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
  
  if (class(x) == "data.frame" | class(x) ==  "matrix"){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  } else if (class(x) == "numeric" | class(x) == "integer" | class(x) == "double") {
    x <- list(data = x)
  }
  
  if(datatype=="abundance"){
    out <- lapply(1:length(x),function(i){
      dq <- c(Diversity_profile(x[[i]],q),Diversity_profile_MLE(x[[i]],q))
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(x[[i]])
        Abun.Mat <- rmultinom(nboot, sum(x[[i]]), Prob.hat)
        
        error <- qnorm(1-(1-conf)/2) * 
          apply(apply(Abun.Mat, 2, function(xb) c(Diversity_profile(xb, q),Diversity_profile_MLE(xb,q))), 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = rep(q,2), "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(x)[i], "method" = rep(c("Estimated","Empirical"),each = length(q)))
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence_freq"){
    out <- lapply(1:length(x),function(i){
      dq <- c(Diversity_profile.inc(x[[i]],q),Diversity_profile_MLE.inc(x[[i]],q))
      if(nboot > 1){
        nT <- x[[i]][1]
        Prob.hat <- EstiBootComm.Sam(x[[i]])
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          error = 0
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          error <- qnorm(1-(1-conf)/2) * 
            apply(apply(Abun.Mat, 2, function(yb) c(Diversity_profile.inc(yb, q),Diversity_profile_MLE.inc(yb,q))), 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = rep(q,2), "qD" = dq,"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(x)[i],"method" = rep(c("Estimated","Empirical"),each = length(q)))
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
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
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggiNEXT(x=out1, type=1)
#' ggiNEXT(x=out1, type=2)
#' ggiNEXT(x=out1, type=3)
#' 
#'\dontrun{
#' # single-assemblage incidence data with three orders q
#' data(ant)
#' size <- round(seq(10, 500, length.out=20))
#' y <- iNEXT(ant$h500m, q=c(0,1,2), datatype="incidence_freq", size=size, se=FALSE)
#' ggiNEXT(y, se=FALSE, color.var="Order.q")
#' 
#' # multiple-assemblage abundance data with three orders q
#' z <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' ggiNEXT(z, facet.var="Assemblage", color.var="Order.q")
#' ggiNEXT(z, facet.var="Both", color.var="Both")
#'}
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
#' out1 <- AsyD(spider, datatype = "abundance")
#' ggAsyD(out1)
#' 
#' ## Type (2) example for incidence-based data
#'
#' ## Ex.2
#' data(ant)
#' out2 <- AsyD(ant, datatype = "incidence_freq", nboot = 0)
#' ggAsyD(out2)
#'
#' @export
ggAsyD <- function(outcome){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  if (sum(unique(outcome$method) %in% c("Estimated", "Empirical")) == 0)
    stop("Please use the outcome from specified function 'AsyD'")
  ggplot(outcome, aes(x=Order.q, y=qD, colour=Assemblage, lty=method)) +
    geom_line(size=1.2) +
    scale_colour_manual(values = cbPalette) +
    geom_ribbon(data = outcome[outcome$method=="Estimated",],
                aes(ymin=qD.LCL, ymax=qD.UCL, fill=Assemblage), alpha=0.2, linetype=0) +
    # geom_ribbon(data = outcome[outcome$method=="Empirical",],
    #             aes(ymin=qD.LCL, ymax=qD.UCL, fill=Assemblage), alpha=0.2, linetype=0) +
    scale_fill_manual(values = cbPalette) +
    scale_linetype_manual(values = c("Estimated"=1, "Empirical"=2)) +
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


