# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{TD.m.est} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @param qs a numerical vector of the order of Hill number
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
TD.m.est = function(x, m, qs){ ## here q is allowed to be a vector containing non-integer values.
  n <- sum(x)
  #xv_matrix = as.matrix(xv)
  ifi <- table(x);ifi <- cbind(i = as.numeric(names(ifi)),fi = ifi)
  obs <- Diversity_profile_MLE(x,qs)
  RFD_m <- RTD(ifi, n, n-1, qs)
  #asymptotic value
  asy <- Diversity_profile(x,qs)
  #beta
  beta <- rep(0,length(qs))
  beta0plus <- which(asy != RFD_m)
  beta[beta0plus] <- (obs[beta0plus]-RFD_m[beta0plus])/(asy[beta0plus]-RFD_m[beta0plus])
  #Extrapolation, 
  ETD = function(m,qs){
    m = m-n
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        1/ ((1/(n+m))+(1-1/(n+m))*sum(ifi[,2]*ifi[,1]/n*(ifi[,1]-1)/(n-1)) )
      } 
    })
    return(out)
  }
  Sub = function(m){
    if(m<n){
      if(m == round(m)) { mRTD[-1,mRTD[1,] == m] 
      } else { (ceiling(m) - m)*mRTD[-1, mRTD[1,] == floor(m)] + (m - floor(m))*mRTD[-1, mRTD[1,] == ceiling(m)] }
    }else if(m==n){
      obs
    }else{
      ETD(m,qs)
    }
  }
  
  if (sum(m < n) != 0) {
    int.m = sort(unique(c(floor(m[m<n]), ceiling(m[m<n]))))
    mRTD = rbind(int.m, sapply(int.m, function(k) RTD(ifi,n,k,qs)))
  }
  as.vector(t(sapply(m, Sub))) 
}


# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{TD.m.est_inc} Estimation of interpolation and extrapolation of incidence-based Hill number
# 
# @param y a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param t_ a integer vector of rarefaction/extrapolation sample size
# @param qs a numerical vector of the order of Hill number
# @return a vector of estimated interpolation and extrapolation function of Hill number with order q
# @export
TD.m.est_inc <- function(y, t_, qs){
  nT <- y[1]
  y <- y[-1]
  y <- y[y > 0]
  U <- sum(y)
  #xv_matrix = as.matrix(xv)
  iQi <- table(y);iQi <- cbind(i = as.numeric(names(iQi)),Qi = iQi)
  obs <- Diversity_profile_MLE.inc(c(nT,y),qs)
  RFD_m <- RTD_inc(iQi, nT, nT-1, qs)
  # RFD_m2 <- RTD_inc(iQi, nT, nT-2, qs)
  # whincr <- which(RFD_m != RFD_m2)
  # Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  asy <- Diversity_profile.inc(c(nT,y),qs)
  beta <- rep(0,length(qs))
  # beta0plus <- which(asy != obs)
  # beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
  beta0plus <- which(asy != RFD_m)
  beta[beta0plus] <- (obs[beta0plus]-RFD_m[beta0plus])/(asy[beta0plus]-RFD_m[beta0plus])
  ETD = function(m,qs){
    m = m-nT
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 ){
        1/ ((1/(nT+m))*(nT/U)+(1-1/(nT+m))*sum(iQi[,2]*iQi[,1]/(U^2)*(iQi[,1]-1)/(1-1/nT)) )
      } 
    })
    return(out)
  }
  Sub = function(m){
    if(m<nT){
      if(m == round(m)) { mRTD_inc[-1,mRTD_inc[1,]==m] 
      } else { (ceiling(m)-m)*mRTD_inc[-1,mRTD_inc[1,]==floor(m)]+(m-floor(m))*mRTD_inc[-1,mRTD_inc[1,]==ceiling(m)] }
    }else if(m==nT){
      obs
    }else{
      ETD(m,qs)
    }
  }
  
  if (sum(t_ < nT) != 0) {
    int.t_ = sort(unique(c(floor(t_[t_<nT]), ceiling(t_[t_<nT]))))
    mRTD_inc = rbind(int.t_, sapply(int.t_, function(k) RTD_inc(iQi,nT,k,qs)))
  }
  as.vector(t(sapply(t_, Sub)))
}


# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{iNEXT.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param Spec a vector of species abundances
# @param q a numerical vector of the order of Hill number
# @param m a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Sam}}
# @examples
# data(spider)
# # q = 0 with specific endpoint
# iNEXT.Ind(spider$Girdled, q=0, endpoint=500)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Ind(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600))
iNEXT.Ind <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), knots=40, nboot=200, conf=0.95, unconditional_var = TRUE)
{
  qtile <- qnorm(1-(1-conf)/2)
  n <- sum(Spec)
  if(is.null(m)) {
    if(endpoint <= n) {
      m <- floor(seq(1, endpoint, length=floor(knots)))
    } else {
      m <- c(floor(seq(1, sum(Spec)-1, length.out=floor(knots/2)-1)), sum(Spec), floor(seq(sum(Spec)+1, to=endpoint, length.out=floor(knots/2))))
    }
    m <- c(1, m[-1])
  } else if(is.null(m)==FALSE) {	
    if(max(m)>n | length(m[m==n])==0)  m <- c(m, n-1, n, n+1)
    m <- sort(m)
  }
  m <- unique(m)
  #====conditional on m====
  Dq.hat <- TD.m.est(Spec,m,q)
  C.hat <- Coverage(Spec, 'abundance', m)
  #====unconditional====
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Ind(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[round(Dq.hat_unc$m) == n] = "Observed"
  }
  
  if(nboot > 1 & length(Spec) > 1) {
    Prob.hat <- EstiBootComm.Ind(Spec)
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(x) TD.m.est(x, m, q)),
                          nrow = length(Dq.hat)),1,sd, na.rm=TRUE)
    
    ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(x) Coverage(x, 'abundance', m)),nrow = length(m)),
                        1, sd, na.rm=TRUE)
    if(unconditional_var){
      ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(x) invChat.Ind(x, q,unique(Dq.hat_unc$SC))$qD),
                            nrow = nrow(Dq.hat_unc)),1,sd, na.rm=TRUE)
    }
  } else {
    ses_m <- rep(NA,length(Dq.hat))
    ses_C_on_m <- rep(NA,length(m))
    if(unconditional_var){
      ses_C <- rep(NA,nrow(Dq.hat_unc))
    }
  }
  out_m <- data.frame("m"=rep(m,length(q)), "qD"=Dq.hat, 
                      "qD.LCL"=Dq.hat-qtile*ses_m,
                      "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), 
                      "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  # out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$m<n, "Rarefaction", ifelse(out_m$m==n, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(m))
  id_m <- match(c("m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- data.frame(Dq.hat_unc, 'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                       'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("SC","m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}


# iNterpolation and EXTrapolation of incidence-based Hill number
# \code{iNEXT.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
# 
# @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numerical vector of the order of Hill number
# @param t a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Ind}}
# @examples
# data(ant)
# # q = 0 with specific endpoint
# iNEXT.Sam(ant$h50m, q=0, endpoint=100)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)))
iNEXT.Sam <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), knots=40, nboot=200, conf=0.95, unconditional_var = TRUE)
{
  qtile <- qnorm(1-(1-conf)/2)
  if(which.max(Spec)!=1) 
    stop("invalid data structure!, first element should be number of sampling units")
  
  nT <- Spec[1]
  if(is.null(t)) {
    if(endpoint <= nT) {
      t <- floor(seq(1, endpoint, length.out=floor(knots)))
    } else {
      t <- c(floor(seq(1, nT-1, length.out=floor(knots/2)-1)), nT, floor(seq(nT+1, to=endpoint, length.out=floor(knots/2))))
    }
    t <- c(1, t[-1])
  } else if(is.null(t)==FALSE) {	
    if(max(t)>nT | length(t[t==nT])==0)  t <- c(t, nT-1, nT, nT+1)
    t <- sort(t)
  }
  t <- unique(t)
  #====conditional on m====
  Dq.hat <- TD.m.est_inc(Spec,t,q)
  C.hat <- Coverage(Spec, "incidence_freq", t)
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[round(Dq.hat_unc$nt) == nT] = "Observed"
  }
  
  if(nboot > 1 & length(Spec) > 2){
    Prob.hat <- EstiBootComm.Sam(Spec)
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
    Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
    tmp <- which(colSums(Abun.Mat)==nT)
    if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
    if(ncol(Abun.Mat)==0){
      out <- cbind("t"=t, "qD"=Dq.hat, "SC"=C.hat)
      warning("Insufficient data to compute bootstrap s.e.")
    }else{		
      ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(y) TD.m.est_inc(y, t, q)),
                            nrow = length(Dq.hat)),1,sd, na.rm=TRUE)
      
      ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(y) Coverage(y, "incidence_freq", t)),nrow = length(t)),
                          1, sd, na.rm=TRUE)
      if(unconditional_var){
        ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(y) invChat.Sam(y, q,unique(Dq.hat_unc$SC))$qD),
                              nrow = nrow(Dq.hat_unc)),1,sd, na.rm=TRUE)
      }
    }
  }else {
    ses_m <- rep(NA,length(Dq.hat))
    ses_C_on_m <- rep(NA,length(t))
    if(unconditional_var){
      ses_C <- rep(NA,nrow(Dq.hat_unc))
    }
  }
  
  out_m <- data.frame("nt"=rep(t,length(q)), "qD"=Dq.hat, 
                      "qD.LCL"=Dq.hat-qtile*ses_m,
                      "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  # out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$nt<nT, "Rarefaction", ifelse(out_m$nt==nT, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(t))
  id_m <- match(c("nt", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- data.frame(Dq.hat_unc,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                        'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("SC","nt", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}


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
  TYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, TYPE))) 
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1) 
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (inherits(x, c("numeric", "integer"))){
    x <- list(data = x)
  }
  if (inherits(x, c("data.frame", "matrix"))){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if (datatype == "abundance") {
    if (inherits(x, "list")) {
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
        est <- cbind(est,s.e.=ses,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
        est
      })
      out <- do.call(rbind,out)
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-4)),(ncol(out)-2),(ncol(out)-1),(ncol(out)-3))]
      rownames(out) <- NULL
      out = out %>% select(c('Assemblage', 'SC', 'm', 'Method', 'Order.q', 'qD', 's.e.', 'qD.LCL', 'qD.UCL'))
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence_freq") {
    if (inherits(x, "list")) {
      Community = rep(names(x),each = length(q)*length(C))
      out <- lapply(x, function(x_){
        est <- invChat.Sam(x_, q, C)
        if (sum(round(est$nt) > 2 * max(x_))>0) 
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
        est <- cbind(est,s.e.=ses,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
      })
      out <- do.call(rbind,out)
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-4)),(ncol(out)-2),(ncol(out)-1),(ncol(out)-3))]
      rownames(out) <- NULL
      out = out %>% select(c('Assemblage', 'SC', 'nt', 'Method', 'Order.q', 'qD', 's.e.', 'qD.LCL', 'qD.UCL'))
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}


invChat.Ind <- function (x, q, C) {
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- sum(x)
  refC <- Coverage(x, 'abundance', n)
  f <- function(m, C) abs(Coverage(x, 'abundance', m) - C)
  mm <- sapply(C, function(cvrg){
    if (refC == cvrg) {
      mm <- n
    }else if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = sum(x))
      mm <- opt$minimum
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
    }
    mm
  })
  mm[mm < 1] <- 1
  SC <- C
  # if (sum(round(mm) > 2 * n)>0) 
  #   warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  
  out <- TD.m.est(x = x,m = mm,qs = q)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(mm,length(q))
  order <- rep(q,each = length(mm))
  SC <- rep(SC,length(q))
  data.frame(m = m,Method = method,Order.q = order,
             SC=SC,qD = out)
}


invChat.Sam <- function (x, q, C) {
  x <- x[x>0] ####added by yhc
  m <- NULL
  n <- max(x)
  refC <- Coverage(x, "incidence_freq", n)
  f <- function(m, C) abs(Coverage(x, "incidence_freq", m) - C)
  mm <- sapply(C, function(cvrg){
    if (refC == cvrg) {
      mm <- n
    }else if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = max(x))
      mm <- opt$minimum
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
    }
    mm
  })
  mm[mm < 1] <- 1
  SC <- C
  # if (sum(round(mm) > 2 * n)>0) 
  #   warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  out <- TD.m.est_inc(y = x,t_ = mm,qs = q)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(mm,length(q))
  order <- rep(q,each = length(mm))
  SC <- rep(SC,length(q))
  data.frame(nt = m,Method = method,Order.q = order,
             SC=SC,qD = out)
  
}


invSize <- function(x, q, datatype="abundance", size=NULL, nboot=0, conf=NULL){
  qtile <- qnorm(1-(1-conf)/2)
  TYPE <- c("abundance", "incidence_freq")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if (inherits(x, c("numeric", "integer"))){
    x <- list(data = x)
  }
  if (inherits(x, c("data.frame", "matrix"))){
    datalist <- lapply(1:ncol(x), function(i) x[,i])
    if(is.null(colnames(x))) names(datalist) <-  paste0("data",1:ncol(x)) else names(datalist) <- colnames(x)
    x <- datalist
  }
  if(datatype=="abundance"){
    if (inherits(x, "list")) {
      Community = rep(names(x),each = length(q)*length(size))
      out <- lapply(x, function(x_){
        est <- invSize.Ind(x_, q, size)
        if(nboot>1){
          Prob.hat <- EstiBootComm.Ind(x_)
          Abun.Mat <- rmultinom(nboot, sum(x_), Prob.hat)
          ses <- apply(matrix(apply(Abun.Mat,2 ,function(a) invSize.Ind(a, q,size)$qD),
                              nrow = length(q) * length(size)),1,sd)
        }else{
          ses <- rep(0,nrow(est))
        }
        est <- cbind(est,s.e.=ses,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
        est
      })
      out <- do.call(rbind,out)
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }else if (datatype == "incidence_freq") {
    if (inherits(x, "list")) {
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
        est <- cbind(est,s.e.=ses,qD.LCL=est$qD-qtile*ses,qD.UCL=est$qD+qtile*ses)
      })
      out <- do.call(rbind,out)
      out$Assemblage <- Community
      out <- out[,c(ncol(out),seq(1,(ncol(out)-1)))]
      rownames(out) <- NULL
    }else {
      stop("Wrong data format, dataframe/matrix or list would be accepted")
    }
  }
  out
}


invSize.Ind <- function(x, q, size){
  m <- NULL # no visible binding for global variable 'm'
  
  n <- sum(x)
  if(is.null(size)){
    size <- sum(x)
  }
  out <- TD.m.est(x = x,m = size,qs = q)
  SC <- Coverage(x,'abundance',size)
  method <- ifelse(size>n,'Extrapolation',ifelse(size<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(size,length(q))
  order <- rep(q,each = length(size))
  SC <- rep(SC,length(q))
  data.frame(m = m,Method = method,Order.q = order,SC=SC,qD = out)
}


invSize.Sam <- function(x, q, size){
  m <- NULL # no visible binding for global variable 'm'
  
  n <- max(x)
  if(is.null(size)){
    size <- max(x)
  }
  out <- TD.m.est_inc(y = x,t_ = size,qs = q)
  SC <- Coverage(x,"incidence_freq",size)
  method <- ifelse(size>n,'Extrapolation',ifelse(size<n,'Rarefaction','Observed'))
  method <- rep(method,length(q))
  m <- rep(size,length(q))
  order <- rep(q,each = length(size))
  SC <- rep(SC,length(q))
  data.frame(nt = m,Method = method,Order.q = order,SC=SC,qD = out)
}



Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  sortx = sort(unique(x))
  tab = table(x)
  Sub_q012 <- function(q){
    if(q==0){
      length(x) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }else if(q==1){
      A <- sum(tab*sortx/n*(digamma(n)-digamma(sortx)))
      B <- TD1_2nd(n,f1,f2)
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(tab[sortx>=q]*exp(lchoose(sortx[sortx>=q],q)-lchoose(n,q)))
      A^(1/(1-q))
    }
  }
  ans <- rep(0,length(q))
  q_part1 = which(abs(q-round(q))==0)
  if(length(q_part1)>0){
    ans[q_part1] <- sapply(q[q_part1], Sub_q012)
  }
  q_part2 <- which(!abs(q-round(q))==0)
  if(length(q_part2)>0){
    ans[q_part2] <- TDq(ifi = cbind(i = sortx, fi = tab),n = n,qs = q[q_part2],f1 = f1,A = p1)
  }
  ans
}


Diversity_profile.inc <- function(data,q){
  nT = data[1]
  Yi = data[-1]
  Yi <- Yi[Yi!=0]
  U <- sum(Yi)
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  Sobs <- length(Yi)
  
  if(Q2>0 & Q1>0){
    A <- 2*Q2/((nT-1)*Q1+2*Q2)
  }
  else if(Q2==0 & Q1>1){
    A <- 2/((nT-1)*(Q1-1)+2)
  }
  else{
    A <- 1
  }
  
  Q0hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)
  B <- sapply(q,function(q) ifelse(A==1,0,(Q1/nT)*(1-A)^(-nT+1)*round((A^(q-1)-sum(sapply(c(0:(nT-1)),function(r) choose(q-1,r)*(A-1)^r))), 12)))
  qD <- (U/nT)^(q/(q-1))*(qTDFUN(q,Yi,nT) + B)^(1/(1-q))
  qD[which(q==0)] = Sobs+Q0hat
  yi <- Yi[Yi>=1 & Yi<=(nT-1)]
  delta <- function(i){
    (yi[i]/nT)*sum(1/c(yi[i]:(nT-1)))
  }
  if(sum(q %in% 1)>0){
    C_ <- ifelse(A==1,0,(Q1/nT)*(1-A)^(-nT+1)*(-log(A)-sum(sapply(c(1:(nT-1)),function(r) (1-A)^r/r))))
    qD[which(q==1)] <- exp((nT/U)*( sum(sapply(c(1:length(yi)),function(i) delta(i))) + C_)+log(U/nT))
  }
  return(qD)
}


Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}


Diversity_profile_MLE.inc <- function(data,q){
  Yi = data[-1]
  U = sum(Yi)
  Yi <- Yi[Yi!=0]
  ai <- Yi/U
  qD = qTD_MLE(q,ai)
  qD[which(q==1)] <- exp(-sum(ai*log(ai)))
  return(qD)
}


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
  return(c(Prob.hat, Prob.hat.Unse))		  							#Output: a vector of estimated relative abundance
}


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
  return(c(Prob.hat, Prob.hat.Unse))									#Output: a vector of estimated detection probability
}



# Estimation of asymptotic diversity
# 
# \code{asyTD} Estimation of species detection distribution to obtain bootstrap s.e.
asyTD = function(data, datatype, q, nboot, conf) {
  
  if(datatype=="abundance"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile(data[[i]],q)
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(data[[i]])
        Abun.Mat <- rmultinom(nboot, sum(data[[i]]), Prob.hat)
        
        mt = apply(Abun.Mat, 2, function(xb) Diversity_profile(xb, q))
        if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
        error <- qnorm(1-(1-conf)/2) * 
          apply(mt, 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq, "s.e." = error/qnorm(1-(1-conf)/2),"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i], "Method" = "Asymptotic")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence_freq"){
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
          
          mt = apply(Abun.Mat, 2, function(yb) Diversity_profile.inc(yb, q))
          if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
          error <- qnorm(1-(1-conf)/2) * 
            apply(mt, 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq, "s.e." = error/qnorm(1-(1-conf)/2),"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i],"Method" = "Asymptotic")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }
  
  return(out)
}



# Observed diversity
# 
# \code{TDinfo} Estimation of species detection distribution to obtain bootstrap s.e.
obsTD = function(data, datatype, q, nboot, conf) {
 
  if(datatype=="abundance"){
    out <- lapply(1:length(data),function(i){
      dq <- Diversity_profile_MLE(data[[i]],q)
      if(nboot > 1){
        Prob.hat <- EstiBootComm.Ind(data[[i]])
        Abun.Mat <- rmultinom(nboot, sum(data[[i]]), Prob.hat)
        
        mt = apply(Abun.Mat, 2, function(xb) Diversity_profile_MLE(xb, q))
        if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
        error <- qnorm(1-(1-conf)/2) * 
          apply(mt, 1, sd, na.rm=TRUE)
        
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq, "s.e." = error/qnorm(1-(1-conf)/2),"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i], "Method" = "Empirical")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }else if(datatype=="incidence_freq"){
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
          
          mt = apply(Abun.Mat, 2, function(yb) Diversity_profile_MLE.inc(yb, q))
          if (!is.matrix(mt)) mt = matrix(mt, nrow = 1)
          error <- qnorm(1-(1-conf)/2) * 
            apply(mt, 1, sd, na.rm=TRUE)
        }
      } else {error = NA}
      out <- data.frame("Order.q" = q, "qD" = dq, "s.e." = error/qnorm(1-(1-conf)/2),"qD.LCL" = dq - error, "qD.UCL" = dq + error,
                        "Assemblage" = names(data)[i],"Method" = "Empirical")
      out$qD.LCL[out$qD.LCL<0] <- 0
      out
    })
    out <- do.call(rbind,out)
  }
  
  return(out)
}



# Data information
# 
# \code{obsTD} Estimation of species detection distribution to obtain bootstrap s.e.
TDinfo = function(data, datatype) {
  Fun.abun <- function(x){
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x==k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x>0)
    f0.hat <- ifelse(f2==0, (n-1)/n*f1*(f1-1)/2, (n-1)/n*f1^2/2/f2)
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
    Q0.hat <- ifelse(Q2==0, (nT-1)/nT*Q1*(Q1-1)/2, (nT-1)/nT*Q1^2/2/Q2)
    A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
    Chat <- round(1 - Q1/U*A,4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  
  if(datatype == "abundance"){
    if(inherits(data, c("numeric", "integer"))){
      out <- matrix(Fun.abun(data), nrow=1)
    }else if(inherits(data, "list")){
      out <- do.call("rbind", lapply(data, Fun.abun))
    } else if(inherits(data, c("matrix", "data.frame"))){
      out <- t(apply(as.matrix(data), 2, Fun.abun))  
    }
    if(nrow(out) > 1 | inherits(data, 'list')){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage", "n", "S.obs", "SC", paste("f",1:10, sep=""))
    }
    out = as.data.frame(out)
    
  }else if(datatype == "incidence_freq"){
    if(inherits(data, c("numeric", "integer"))){
      out <- matrix(Fun.ince(data), nrow=1)
    }else if(inherits(data, "list")){
      out <- do.call("rbind", lapply(data, Fun.ince))
    } else if(inherits(data, c("matrix", "data.frame"))){
      out <- t(apply(as.matrix(data), 2, Fun.ince))  
    }
    if(nrow(out) > 1 | inherits(data, 'list')){
      out <- data.frame(site=rownames(out), out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
      rownames(out) <- NULL
    }else{
      out <- data.frame(site="site.1", out)
      colnames(out) <-  c("Assemblage","T", "U", "S.obs", "SC", paste("Q",1:10, sep=""))
    }
    out = as.data.frame(out)
  }
}
