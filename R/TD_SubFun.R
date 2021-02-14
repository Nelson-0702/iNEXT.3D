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
  # RFD_m2 <- RTD(ifi, n, n-2, qs)
  # whincr <- which(RFD_m != RFD_m2)
  # Dn1 <- obs; Dn1[whincr] <- obs + (obs - RFD_m)^2/(RFD_m - RFD_m2)
  #asymptotic value
  asy <- Diversity_profile(x,qs)
  #beta
  beta <- rep(0,length(qs))
  # beta0plus <- which(asy != obs)
  # beta[beta0plus] <- (Dn1[beta0plus]-obs[beta0plus])/(asy[beta0plus]-obs[beta0plus])
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
      if(m == round(m)) { mRTD[-1,mRTD[1,]==m] 
      } else { (ceiling(m)-m)*mRTD[-1,mRTD[1,]==floor(m)]+(m-floor(m))*mRTD[-1,mRTD[1,]==ceiling(m)] }
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
  out_m <- cbind("m"=rep(m,length(q)), "qD"=Dq.hat, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), 
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$m<n, "Rarefaction", ifelse(out_m$m==n, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(m))
  id_m <- match(c("m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("goalSC","SC","m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
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
  #====unconditional====
  # if(unconditional_var){
  #   goalSC <- unique(round(C.hat,4))
  #   goalSC[goalSC==1] <- 0.9999
  #   Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
  # }
  if(unconditional_var){
    goalSC <- unique(C.hat)
    Dq.hat_unc <- unique(invChat.Sam(x = Spec,q = q,C = goalSC))
    Dq.hat_unc$Method[Dq.hat_unc$t == nT] = "Observed"
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
  
  out_m <- cbind("t"=rep(t,length(q)), "qD"=Dq.hat, "qD.LCL"=Dq.hat-qtile*ses_m,
                 "qD.UCL"=Dq.hat+qtile*ses_m,"SC"=rep(C.hat,length(q)), 
                 "SC.LCL"=C.hat-qtile*ses_C_on_m, "SC.UCL"=C.hat+qtile*ses_C_on_m)
  out_m <- data.frame(out_m)
  out_m$Method <- ifelse(out_m$t<nT, "Rarefaction", ifelse(out_m$t==nT, "Observed", "Extrapolation"))
  out_m$Order.q <- rep(q,each = length(t))
  id_m <- match(c("t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL"), names(out_m), nomatch = 0)
  out_m <- out_m[, id_m]
  out_m$qD.LCL[out_m$qD.LCL<0] <- 0
  out_m$SC.LCL[out_m$SC.LCL<0] <- 0
  out_m$SC.UCL[out_m$SC.UCL>1] <- 1
  
  if(unconditional_var){
    out_C <- cbind(Dq.hat_unc,'qD.LCL' = Dq.hat_unc$qD-qtile*ses_C,
                   'qD.UCL' = Dq.hat_unc$qD+qtile*ses_C) 
    id_C <- match(c("goalSC","SC","t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(out_C), nomatch = 0)
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
  B <- sapply(q,function(q) ifelse(A==1,0,(Q1/nT)*(1-A)^(-nT+1)*(A^(q-1)-sum(sapply(c(0:(nT-1)),function(r) choose(q-1,r)*(A-1)^r)))))
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


# ggiNEXT.iNEXT -------------------------------------------------------------------

#' @export
#' @rdname ggiNEXT
ggiNEXT.iNEXT <- function(x, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7",  "#0072B2", "#D55E00"))
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("None", "Order.q", "Assemblage", "Both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")
  
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="Order.q") color.var <- "Assemblage"
  if(facet.var=="Assemblage") color.var <- "Order.q"
  
  options(warn = -1)
  z <- fortify(x, type=type)
  options(warn = 0)
  if(!('y.lwr' %in% names(z))) { se <- FALSE }
  datatype <- unique(z$datatype)
  if(color.var=="None"){
    if(levels(factor(z$Order.q))>1 & length(unique(z$Assemblage))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple assemblages and orders, change setting as Both")
      color.var <- "Both"
      z$col <- z$shape <- paste(z$Assemblage, z$Order.q, sep="-")
      
    }else if(length(unique(z$Assemblage))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple assemblages, change setting as Assemblage")
      color.var <- "Assemblage"
      z$col <- z$shape <- z$Assemblage
    }else if(levels(factor(z$Order.q))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as Order.q")
      color.var <- "Order.q"
      z$col <- z$shape <- factor(z$Order.q)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="Order.q"){     
    z$col <- z$shape <- factor(z$Order.q)
  }else if(color.var=="Assemblage"){
    if(length(unique(z$Assemblage))==1){
      warning("invalid color.var setting, the iNEXT object do not consist multiple assemblages, change setting as Order.q")
      z$col <- z$shape <- factor(z$Order.q)
    }
    z$col <- z$shape <- z$Assemblage
  }else if(color.var=="Both"){
    if(length(unique(z$Assemblage))==1){
      warning("invalid color.var setting, the iNEXT object do not consist multiple assemblages, change setting as Order.q")
      z$col <- z$shape <- factor(z$Order.q)
    }
    z$col <- z$shape <- paste(z$Assemblage, z$Order.q, sep="-")
  }
  zz=z
  z$Method[z$Method=="Observed"]="Rarefaction"
  z$lty <- factor(z$Method, levels = c("Rarefaction", "Extrapolation"))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$Method=="Observed"),]
  
  g <- ggplot(z, aes_string(x="x", y="y", colour="col")) + 
    geom_point(aes_string(shape="shape"), size=5, data=data.sub)+
    scale_colour_manual(values=cbPalette)
  
  
  g <- g + geom_line(aes_string(linetype="lty"), lwd=1.5) +
    guides(linetype=guide_legend(title="Method"),
           colour=guide_legend(title="Guides"), 
           fill=guide_legend(title="Guides"), 
           shape=guide_legend(title="Guides")) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          legend.key.width = unit(1.2,"cm")) 
  
  if(type==2L) {
    g <- g + labs(x="Number of sampling units", y="Sample coverage")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Sample coverage")
  }else if(type==3L|type==4L) {
    g <- g + labs(x="Sample coverage", y="Species diversity")
  }else {
    g <- g + labs(x="Number of sampling units", y="Species diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Species diversity")
  }
  
  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.2)+
    scale_fill_manual(values=cbPalette)
  
  
  if(facet.var=="Order.q"){
    if(length(levels(factor(z$Order.q))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")      
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(~Order.q, nrow=1, labeller = odr_grp)
      if(color.var=="Both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$Order.q))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
      if(type==2){
        g <- g + theme(strip.background = element_blank(),strip.text.x = element_blank())
        
      }
    }
  }
  
  if(facet.var=="Assemblage"){
    if(length(unique(z$Assemblage))==1) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple assemblages")
    }else{
      g <- g + facet_wrap(~Assemblage, nrow=1)
      if(color.var=="Both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$Order.q)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(facet.var=="Both"){
    if(length(levels(factor(z$Order.q))) == 1 | length(unique(z$Assemblage))==1){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple assemblages or orders.")
    }else{
      odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2")) 
      g <- g + facet_wrap(Assemblage~Order.q,labeller = labeller(Order.q = odr_grp)) 
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$Assemblage))), byrow=TRUE),
                         fill=guide_legend(title="Guides"))
      }
    }
  }
  
  if(grey){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      guides(linetype=guide_legend(title="Method"), 
             colour=guide_legend(title="Guides"), 
             fill=guide_legend(title="Guides"), 
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank())
  }
  g <- g + theme(legend.box = "vertical")
  return(g)
  
}


# ggiNEXT.default -------------------------------------------------------------------
#' @export
#' @rdname ggiNEXT
ggiNEXT.default <- function(x, ...){
  stop(
    "iNEXT doesn't know how to deal with data of class ",
    paste(class(x), collapse = "/"),
    call. = FALSE
  )
}


# fortify.iNEXT -------------------------------------------------------------------
#' Fortify method for classes from the iNEXT package.
#'
#' @name fortify.iNEXT
#' @param model \code{iNEXT} to convert into a dataframe.
#' @param data not used by this method
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' @param ... not used by this method
#' @import ggplot2
#' @export
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggplot2::fortify(out1, type=1)
fortify.iNEXT <- function(model, data = model$iNextEst, type = 1, ...) {
  datatype <- ifelse(names(model$DataInfo)[2]=="n","abundance","incidence")
  z <- data
  # if(class(z) == "list"){
  #   if(datatype=='abundance'){
  #     id_match <- match(c("Assemblage","m", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC"), names(z$coverage_based), nomatch = 0)  
  #   }else if (datatype=='incidence'){
  #     id_match <- match(c("Assemblage","t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC"), names(z$coverage_based), nomatch = 0)  
  #   }
  #   z$coverage_based <- cbind(z$coverage_based[,id_match],SC.LCL=NA,SC.UCL=NA)
  #   z <- data.frame(do.call("rbind", z), base=rep(names(z), each = sapply(z, nrow)))
  #   rownames(z) <- NULL
  # }else{
  #   z$site <- ""
  # }
  
  # if(ncol(z)==6) {
  #   warning("invalid se setting, the iNEXT object do not consist confidence interval")
  #   se <- FALSE
  # }else if(ncol(z)>6) {
  #   se <- TRUE
  # }
  
  if(is.na(z$size_based$qD.LCL[1])) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else{
    se <- TRUE
  }
  
  if(type==1L) {
    z <- z$size_based
    z$x <- z[,2]
    z$y <- z$qD
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }else if(type==2L){
    z <- z$size_based
    if(length(unique(z$Order.q))>1){
      z <- subset(z, Order.q==unique(z$Order.q)[1])
    }
    z$x <- z[,2]
    z$y <- z$SC
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$SC.LCL
      z$y.upr <- z$SC.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }else if(type==3L){
    z <- z$coverage_based
    z$x <- z$SC
    z$y <- z$qD
    z$datatype <- datatype
    z$plottype <- type
    if(se){
      z$y.lwr <- z$qD.LCL
      z$y.upr <- z$qD.UCL
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
    }else{
      data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
    }
  }
  # else if(type==4L){
  #   z <- z$size_based
  #   z$x <- z$SC
  #   z$y <- z$qD
  #   z$datatype <- datatype
  #   z$plottype <- type
  #   if(se){
  #     z$y.lwr <- z$qD.LCL
  #     z$y.upr <- z$qD.UCL
  #     data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y","y.lwr","y.upr")]
  #   }else{
  #     data <- z[,c("datatype","plottype","Assemblage","Method","Order.q","x","y")]
  #   }
  # }
  data
}


# plot.iNEXT -------------------------------------------------------------------
#' Plotting iNEXT object
#' 
#' \code{plot.iNEXT}: Plotting method for objects inheriting from class "iNEXT"
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).                 
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param show.legend a logical variable to display legend.
#' @param show.main a logical variable to display title.
#' @param col a vector for plotting color
#' @param ... arguments to be passed to methods, such as graphical parameters (\code{\link{par}}).
#' @examples
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' plot(x=out1, type=1)
#' plot(x=out1, type=2)
#' plot(x=out1, type=3)
#' 

#' @export
plot.iNEXT <- function(x, type=1, se=TRUE, show.legend=TRUE, show.main=TRUE, col=NULL,...){
  
  if(class(x) != "iNEXT")
    stop("invalid object class")
  TYPE <-  c(1, 2, 3)
  # SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  
  type <- pmatch(type, 1:3)
  
  y <- method <- site <- shape <- y.lwr <- y.upr <- NULL
  site <<- NULL
  
  # z <- x$iNextEst
  # if(class(z) == "list"){
  #   z <- data.frame(do.call("rbind", z), site=rep(names(z), sapply(z, nrow)))
  #   rownames(z) <- NULL
  # }else{
  #   z$site <- ""
  #   z$site <- factor(z$site)
  # }
  
  z <- fortify(x, type=type)
  
  
  if("y.lwr" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if("y.lwr" %in% names(z) & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  if(type==1L) {
    #z$x <- z[,1]
    #z$y <- z$qD
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Species diversity"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
  }else if(type==2L){
    if(length(unique(z$Order.q))>1){
      z <- subset(z, Order.q==unique(z$Order.q)[1])
    }
    # z$x <- z[,1]
    # z$y <- z$SC
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Sample coverage"
    # if(se){
    #   z$y.lwr <- z[,8]
    #   z$y.upr <- z[,9]
    # }
  }else if(type==3L){
    # z$x <- z$SC
    # z$y <- z$qD
    if(!is.null(xlab)) xlab <- "Sample coverage"
    if(!is.null(ylab)) ylab <- "Species diversity"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  conf.reg=function(x,LCL,UCL,...) {
    x.sort <- order(x)
    x <- x[x.sort]
    LCL <- LCL[x.sort]
    UCL <- UCL[x.sort]
    polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
  }
  
  SITE <- unique(z$Assemblage)
  ORDER <- unique(z$Order.q)
  
  if(is.null(col)){
    col <- gg_color_hue(length(SITE))
  }else{
    col <- rep(col,length(SITE))[1:length(SITE)]
  }
  pch <- (16+1:length(SITE))%%25
  
  for(j in 1:length(ORDER)){
    if(se==TRUE){
      tmp.sub <- subset(z, Order.q==ORDER[j])
      tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage, Order.q=tmp.sub$Order.q,
                          Method=tmp.sub$Method, 
                          x=tmp.sub$x, y=tmp.sub$y,
                          y.lwr=tmp.sub$y.lwr, y.upr=tmp.sub$y.upr)
      
      plot(y.upr~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }else{
      tmp.sub <- subset(z, Order.q==ORDER[j])
      
      tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage, Order.q=tmp.sub$Order.q,
                          Method=tmp.sub$Method, 
                          x=tmp.sub$x, y=tmp.sub$y)
      
      plot(y~x, data=tmp.j, type="n", xlab="", ylab="", ...)
    }
    
    for(i in 1:length(SITE)){
      tmp <- subset(tmp.j, Assemblage==SITE[i])
      if(se==TRUE){
        conf.reg(x=tmp$x, LCL=tmp$y.lwr, UCL=tmp$y.upr, border=NA, col=adjustcolor(col[i], 0.25))
      }
      lines(y~x, data=subset(tmp, Method=="Rarefaction"), lty=1, lwd=2, col=col[i])
      lines(y~x, data=subset(tmp, Method=="Extrapolation"), lty=2, lwd=2, col=col[i])
      points(y~x, data=subset(tmp, Method=="Observed"), pch=pch[i], cex=2, col=col[i])
      
    }
    if(show.legend==TRUE){
      if(type==3L){
        legend("topleft", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }else{
        legend("bottomright", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
      }
    }
    title(xlab=xlab, ylab=ylab)
    if(show.main==TRUE) title(main=paste("Order q =", ORDER[j]))
    par(ask=TRUE)
  }
  par(ask=FALSE)
}


# print.iNEXT -------------------------------------------------------------------
#' Printing iNEXT object
#' 
#' \code{print.iNEXT}: Print method for objects inheriting from class "iNEXT"
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param ... additional arguments.
#' @export
print.iNEXT <- function(x, ...){
  site.n <- nrow(x$DataInfo)
  order.n <- paste(unique(x$iNextEst$size_based$Order.q), collapse = ", ")
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep="")
  cat("$class: iNEXT\n\n")
  cat("$DataInfo: basic data information\n")
  print(x$DataInfo)
  cat("\n")
  cat("$iNextEst: diversity estimates with rarefied and extrapolated samples.\n")
  if(class(x$iNextEst)=="data.frame"){
    y <- x$iNextEst
    m <- quantile(y[,1], type = 1)
    res <- y[y[,1]%in%m,]
  }else{
    res <- lapply((x$iNextEst), function(y){
      Assemblages <- unique(x$iNextEst$size_based$Assemblage)
      tmp <- lapply(1:length(Assemblages),function(i){
        y_each <- subset(y, Assemblage==Assemblages[i])
        m <- quantile(y_each[,2], type = 1)
        y_each[y_each[,2]%in%m,]
      })
      do.call(rbind,tmp)
    })
  }
  print(res)
  cat("\n")
  cat("$AsyEst: asymptotic diversity estimates along with related statistics.\n")
  print(x$AsyEst)
  cat("\n")
  cat("NOTE: Only show five estimates, call iNEXT.object$iNextEst. to show complete output.\n")
  return(invisible())
}


