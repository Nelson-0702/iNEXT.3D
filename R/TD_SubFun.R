# TD.m.est -------------------------------------------------------------------
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


# TD.m.est_inc -------------------------------------------------------------------
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


# iNEXT.Ind -------------------------------------------------------------------
# iNterpolation and EXTrapolation of abundance-based Hill number
# 
# \code{iNEXT.Ind} Estimation of interpolation and extrapolation of abundance-based Hill number with order q
# 
# @param Spec a vector of species abundances
# @param q a numerical vector of the order of Hill number
# @param m a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Sam}}
# @examples
# data(spider)
# # q = 0 with specific endpoint
# iNEXT.Ind(spider$Girdled, q=0, endpoint=500)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Ind(spider$Girdled, q=1, m=c(1, 10, 20, 50, 100, 200, 400, 600), se=FALSE)
iNEXT.Ind <- function(Spec, q=0, m=NULL, endpoint=2*sum(Spec), knots=40, se=TRUE, nboot=200, conf=0.95, unconditional_var = TRUE)
{
  qtile <- qnorm(1-(1-conf)/2)
  n <- sum(Spec)		  	#sample size
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
  C.hat <- Chat.Ind(Spec, m)
  #====unconditional====
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Ind(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[round(Dq.hat_unc$m) == n] = "Observed"
  }
  
  if(se==TRUE & nboot > 1 & length(Spec) > 1) {
    Prob.hat <- EstiBootComm.Ind(Spec)
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    ses_m <- apply(matrix(apply(Abun.Mat,2 ,function(x) TD.m.est(x, m, q)),
                          nrow = length(Dq.hat)),1,sd, na.rm=TRUE)
    
    ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)),nrow = length(m)),
                        1, sd, na.rm=TRUE)
    if(unconditional_var){
      ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(x) invChat.Ind(x, q,unique(Dq.hat_unc$goalSC))$qD),
                            nrow = nrow(Dq.hat_unc)),1,sd, na.rm=TRUE)
    }
  } else {
    ses_m <- rep(NA,length(Dq.hat))
    ses_C_on_m <- rep(NA,length(m))
    if(unconditional_var){
      ses_C <- rep(NA,nrow(Dq.hat_unc))
    }
  }
  out_m <- cbind("m"=rep(m,length(q)), "qD"=Dq.hat, 
                 "s.e."=ses_m, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), "SC.s.e." = ses_C_on_m,
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$m<n, "Rarefaction", ifelse(out_m$m==n, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(m))
  id_m <- match(c("m", "Method", "Order.q", "qD", "s.e.", "qD.LCL", "qD.UCL", "SC", "SC.s.e.", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc, 's.e.' = ses_C,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("goalSC","SC","m", "Method", "Order.q", "qD", "s.e.", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}


# iNEXT.Sam -------------------------------------------------------------------
# iNterpolation and EXTrapolation of incidence-based Hill number
# 
# \code{iNEXT.Sam} Estimation of interpolation and extrapolation of incidence-based Hill number with order q
# 
# @param Spec a vector of species incidence-based frequency, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param q a numerical vector of the order of Hill number
# @param t a integer vector of rarefaction/extrapolation sample size, default is NULL. If m is not be specified, then the program will compute sample units due to endpoint and knots.
# @param endpoint a integer of sample size that is the endpoint for rarefaction/extrapolation. Default is double the original sample size.
# @param knots a number of knots of computation, default is 40
# @param se calculate bootstrap standard error and 95% confidence interval; default is TRUE
# @param nboot the number of bootstrap resampling times, default is 200
# @return a list of interpolation and extrapolation Hill number with specific order q (qD) and sample coverage (SC)
# @seealso \code{\link{iNEXT.Ind}}
# @examples
# data(ant)
# # q = 0 with specific endpoint
# iNEXT.Sam(ant$h50m, q=0, endpoint=100)
# # q = 1 with specific sample size m and don't calculate standard error
# iNEXT.Sam(ant$h500m, q=1, t=round(seq(10, 500, length.out=20)), se=FALSE)
iNEXT.Sam <- function(Spec, t=NULL, q=0, endpoint=2*max(Spec), knots=40, se=TRUE, nboot=200, conf=0.95, unconditional_var = TRUE)
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
  C.hat <- Chat.Sam(Spec, t)
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[round(Dq.hat_unc$nt) == nT] = "Observed"
  }
  
  if(se==TRUE & nboot > 1 & length(Spec) > 2){
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
      
      ses_C_on_m <- apply(matrix(apply(Abun.Mat, 2, function(y) Chat.Sam(y, t)),nrow = length(t)),
                          1, sd, na.rm=TRUE)
      if(unconditional_var){
        ses_C <- apply(matrix(apply(Abun.Mat,2 ,function(y) invChat.Sam(y, q,unique(Dq.hat_unc$goalSC))$qD),
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
  
  out_m <- cbind("nt"=rep(t,length(q)), "qD"=Dq.hat, 
                 "s.e."=ses_m, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), "SC.s.e." = ses_C_on_m,
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$nt<nT, "Rarefaction", ifelse(out_m$nt==nT, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(t))
  id_m <- match(c("nt", "Method", "Order.q", "qD", "s.e.", "qD.LCL", "qD.UCL", "SC", "SC.s.e.", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc,'s.e.'=ses_C,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("goalSC","SC","nt", "Method", "Order.q", "qD", "s.e.", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
    out_C <- out_C[, id_C]
    out_C$qD.LCL[out_C$qD.LCL<0] <- 0
  }else{
    out_C <- NULL
  }
  return(list(size_based = out_m, coverage_based = out_C))
}


# Diversity_profile -------------------------------------------------------------------
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


# Diversity_profile.inc -------------------------------------------------------------------
Diversity_profile.inc <- function(data,q){
  nT = data[1]
  Yi = data[-1]
  Yi <- Yi[Yi!=0]
  U <- sum(Yi)
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  Sobs <- length(Yi)
  A <- AA.inc(data)
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


# Diversity_profile_MLE -------------------------------------------------------------------
Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}


# Diversity_profile_MLE.inc -------------------------------------------------------------------
Diversity_profile_MLE.inc <- function(data,q){
  Yi = data[-1]
  U = sum(Yi)
  Yi <- Yi[Yi!=0]
  ai <- Yi/U
  qD = qTD_MLE(q,ai)
  qD[which(q==1)] <- exp(-sum(ai*log(ai)))
  return(qD)
}


# AA.inc -------------------------------------------------------------------
AA.inc <- function(data){
  nT = data[1]
  U <- sum(data[-1])
  data = data[-1]
  Yi = data[data!=0]
  Q1 <- sum(Yi==1)
  Q2 <- sum(Yi==2)
  if(Q2>0 & Q1>0){
    A <- 2*Q2/((nT-1)*Q1+2*Q2)
  }
  else if(Q2==0 & Q1>1){
    A <- 2/((nT-1)*(Q1-1)+2)
  }
  else{
    A <- 1
  }
  return(A)
}



