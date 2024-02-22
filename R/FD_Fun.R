data_transform <- function(data, dij, tau, datatype, integer = FALSE, truncate = TRUE, filt_zero = TRUE){
  if(datatype == 'abundance'){
    if (filt_zero) {
      dij <- dij[data>0,data>0]
      data <- data[data>0] 
    }
    out <- lapply(tau,function(tau_){
      dij_ <- dij
      if(tau_==0){
        dij_[dij_>0] <- 1
        a <- as.vector((1 - dij_/1) %*% data )
      }else{
        dij_[which(dij_>tau_,arr.ind = T)] <- tau_
        a <- as.vector((1 - dij_/tau_) %*% data )
      }
      if (filt_zero) {
        data <- data[a!=0]
        a <- a[a!=0]
      }
      a[a<1] <- 1
      if (integer) a <- round(a)
      v <- data/a
      cbind(a,v)
    })
  }else if (datatype == 'incidence_freq'){
    nT = data[1]
    data <- data[-1]
    if (filt_zero) {
      dij <- dij[data>0,data>0]
      data <- data[data>0]
    }
    out <- lapply(tau,function(tau_){
      dij_ <- dij
      if(tau_==0){
        dij_[dij_>0] <- 1
        a <- as.vector((1 - dij_/1) %*% data )
      }else{
        dij_[which(dij_>tau_,arr.ind = T)] <- tau_
        a <- as.vector((1 - dij_/tau_) %*% data )
      }
      if (filt_zero) {
        data <- data[a!=0]
        a <- a[a!=0]
      }
      a[a<1] <- 1
      if (truncate) a[a>nT] <- nT
      if (integer) a <- round(a)
      v <- data/a
      cbind(a,v)
    })
  }
  out_a <- matrix(sapply(out, function(x) x[,1]),ncol = length(tau))
  out_v <- matrix(sapply(out, function(x) x[,2]),ncol = length(tau))
  colnames(out_a) <- colnames(out_v) <- paste0('tau_',round(tau,3))
  
  output = list(ai = out_a, vi = out_v)
  output
}


FD.m.est = function(ai_vi, m, q, nT, ai_vi_MLE){
  EFD = function(m,qs,obs,asy,beta,av){
    m = m-nT
    out <- sapply(1:length(qs), function(i){
      if( qs[i] != 2 ) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( qs[i] == 2 & beta[i] != 0 ){
        V_bar^2/sum( (av[,2])*((1/(nT+m))*(av[,1]/nT)+((nT+m-1)/(nT+m))*(av[,1]*(av[,1]-1)/(nT*(nT-1)))) )
      }else if( qs[i] == 2 & beta[i] == 0 ){
        asy[i]
      } 
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])/nT
  asy <- FD_est(ai_vi,q,nT,ai_vi_MLE)
  obs <- FD_mle(ai_vi,q)
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    av = cbind(ai = ai_vi$ai[,i], vi = ai_vi$vi[,i])
    RFD_m = RFD(av, nT, nT-1, q, V_bar)
    beta <- rep(0,length(q))
    #asymptotic value; observed value
    asy_i <- asy[,i];obs_i <- obs[,i]
    asy_i <- sapply(1:length(q), function(j){
      max(asy_i[j],obs_i[j])
    })
    RFD_m[RFD_m > obs_i] = obs_i[RFD_m > obs_i]
    beta0plus <- which( asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus]-RFD_m[beta0plus])/(asy_i[beta0plus]-RFD_m[beta0plus])
    
    if (sum(m < nT) != 0) {
      int.m = sort(unique(c(floor(m[m<nT]), ceiling(m[m<nT]))))
      mRFD = rbind(int.m, sapply(int.m, function(k) RFD(av,nT,k,q,V_bar)))
      
      if (0 %in% int.m) mRFD[,mRFD[1,] == 0] = 0
    }
    
    sapply(m, function(mm){
      if(mm<nT){
        if(mm == round(mm)) { mRFD[-1,mRFD[1,] == mm] 
        } else { (ceiling(mm) - mm)*mRFD[-1, mRFD[1,] == floor(mm)] + (mm - floor(mm))*mRFD[-1, mRFD[1,] == ceiling(mm)] }
      }else if(mm==nT){
        obs_i
      }else if(mm==Inf){
        asy_i
      }else{
        EFD(m = mm,qs = q,obs = obs_i,asy = asy_i,beta = beta,av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out,ncol = ncol(ai_vi$ai))
}


iNextFD = function(datalist, dij, q = c(0,1,2), datatype, tau, nboot, conf = 0.95, m){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  length_tau <- length(tau)
  
  if(datatype=="abundance"){
    
    out <- lapply(1:length(datalist), function(i){
      x <- datalist[[i]]
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      qFDm <- FD.m.est(ai_vi = data_aivi,m = m[[i]],q = q,nT = n, ai_vi_MLE = data_aivi_MLE) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          qFDm_b <- FD.m.est(ai_vi = Boot_aivi,m = m[[i]],q = q,nT = n,ai_vi_MLE = Boot_aivi_MLE) %>%
            as.numeric()
          covm_b = Coverage(Boot.X[,B], datatype, m[[i]])
          return(c(qFDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(c(qFDm,covm)))
      }
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,length(q)*length(tau))
      orderq <- rep(rep(q,each = length(m[[i]])),length(tau))
      threshold <- rep(tau,each = length(q)*length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_cov <- rep(ses_cov,each = length(q)*length(tau))
      ses_fd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      covm <- rep(covm,length(q)*length(tau))
      data.frame(Assemblage = sites[i], Order.q=orderq, m=rep(m[[i]],length(q)*length(tau)), Method=method, 
                 qFD=qFDm, s.e.=ses_fd, qFD.LCL=qFDm-qtile*ses_fd, qFD.UCL=qFDm+qtile*ses_fd,
                 SC=covm, SC.s.e.=ses_cov, SC.LCL=covm-qtile*ses_cov, SC.UCL=covm+qtile*ses_cov,
                 Tau = threshold) %>% 
        arrange(Tau, Order.q)
    }) %>% do.call(rbind, .)
    
  }else if(datatype=="incidence_freq"){
    
    out <- lapply(1:length(datalist), function(i){
      x <- datalist[[i]]
      nT=x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      qFDm <- FD.m.est(ai_vi = data_aivi,m = m[[i]],q = q,nT = nT,ai_vi_MLE = data_aivi_MLE) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          qFDm_b <- FD.m.est(ai_vi = Boot_aivi,m = m[[i]],q = q,nT = nT,ai_vi_MLE = Boot_aivi_MLE) %>% as.numeric()
          covm_b = Coverage(Boot.X, datatype, m[[i]])
          return(c(qFDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(c(qFDm,covm)))
      }
      method <- ifelse(m[[i]]>nT,'Extrapolation',ifelse(m[[i]]<nT,'Rarefaction','Observed'))
      method <- rep(method,length(q)*length(tau))
      orderq <- rep(rep(q,each = length(m[[i]])),length(tau))
      threshold <- rep(tau,each = length(q)*length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_cov <- rep(ses_cov,each = length(q)*length(tau))
      ses_fd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      covm <- rep(covm,length(q)*length(tau))
      data.frame(Assemblage = sites[i], Order.q=orderq, m=rep(m[[i]],length(q)*length(tau)), Method=method,
                 qFD=qFDm, s.e.=ses_fd, qFD.LCL=qFDm-qtile*ses_fd, qFD.UCL=qFDm+qtile*ses_fd,
                 SC=covm, SC.s.e.=ses_cov, SC.LCL=covm-qtile*ses_cov, SC.UCL=covm+qtile*ses_cov,
                 Tau = threshold) %>% 
        arrange(Tau,Order.q)
    }) %>% do.call(rbind, .)
  }
  
  out$qFD.LCL[out$qFD.LCL<0] <- 0
  out$SC.LCL[out$SC.LCL<0] <- 0
  out$SC.UCL[out$SC.UCL>1] <- 1
  return(out)
}


invChatFD <- function(datalist, dij, q, datatype, level, nboot, conf = 0.95, tau){
  qtile <- qnorm(1-(1-conf)/2)
  
  if(datatype=='abundance'){
    out <- lapply(datalist,function(x_){
      data_aivi <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype)
      
      est <- invChatFD_abu(ai_vi = data_aivi,data_ = x_,q = q,Cs = level,tau = tau, ai_vi_MLE = data_aivi_MLE)
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x_,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        n=sum(x_)
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          
          invChatFD_abu(ai_vi = Boot_aivi,data_ = Boot.X[,B],q = q,Cs = level,tau = tau, ai_vi_MLE = Boot_aivi_MLE)$qFD
        })
        
        if (length(q) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(s.e.=ses,qFD.LCL=qFD-qtile*ses,qFD.UCL=qFD+qtile*ses) 
    }) %>% do.call(rbind,.)
  }else if(datatype=='incidence_freq'){
    out <- lapply(datalist,function(x_){
      nT=x_[1]
      data_aivi <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x_,dij = dij,tau = tau,datatype = datatype)
      
      est <- invChatFD_inc(ai_vi = data_aivi,data_ = x_,q = q,Cs = level,tau = tau, ai_vi_MLE = data_aivi_MLE)
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x_,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          
          invChatFD_inc(ai_vi = Boot_aivi,data_ = Boot.X,q = q,Cs = level,tau = tau, ai_vi_MLE = Boot_aivi_MLE)$qFD
        })
        
        if (length(q) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(s.e.=ses,qFD.LCL=qFD-qtile*ses,qFD.UCL=qFD+qtile*ses) 
    }) %>% do.call(rbind,.)
  }
  Assemblage = rep(names(datalist), each = length(q)*length(level)*length(tau))
  out <- out %>% mutate(Assemblage = Assemblage) %>% select(
    Assemblage, Order.q, SC, m, Method, qFD, s.e., qFD.LCL, qFD.UCL, Tau
  )
  rownames(out) <- NULL
  out
}


invChatFD_abu <- function(ai_vi, data_, q, Cs, tau, ai_vi_MLE){
  n <- sum(data_)
  refC = Coverage(data_, 'abundance', n)
  f <- function(m, cvrg) abs(Coverage(data_, 'abundance', m) - cvrg)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, cvrg = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      if (cvrg == 0) mm = 0
    }else if (refC <= cvrg) {
      f1 <- sum(data_ == 1)
      f2 <- sum(data_ == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 0 & f2 > 0) {
        A <- 0
      }
      if (f1 == 1 & f2 == 0) {
        A <- 0
      }
      if (f1 == 0 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A==0,0,(log(n/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
    }
    mm
  })
  # mm[mm < 1] <- 1
  SC <- Cs
  out <- FD.m.est(ai_vi = ai_vi,m = mm,q = q,nT = n,ai_vi_MLE = ai_vi_MLE)
  out <- as.vector(out)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q)*length(tau))
  m <- rep(mm,length(q)*length(tau))
  order <- rep(rep(q,each = length(mm)),length(tau))
  SC <- rep(SC,length(q)*length(tau))
  threshold <- rep(tau,each = length(q)*length(mm))
  method[SC == refC] = "Observed"
  data.frame(m = m,Method = method,Order.q = order,
             qFD = out,SC=SC, Tau = threshold)
}


invChatFD_inc <- function(ai_vi, data_, q, Cs, tau, ai_vi_MLE){
  n <- data_[1]
  refC = Coverage(data_, 'incidence_freq', n)
  f <- function(m, cvrg) abs(Coverage(data_, 'incidence_freq', m) - cvrg)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, cvrg = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
      if (cvrg == 0) mm = 0
    }else if (refC <= cvrg) {
      f1 <- sum(data_ == 1)
      f2 <- sum(data_ == 2)
      U <- sum(data_[-1])
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 0 & f2 > 0) {
        A <- 0
      }
      if (f1 == 1 & f2 == 0) {
        A <- 0
      }
      if (f1 == 0 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A==0,0,(log(U/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
    }
    mm
  })
  # mm[mm < 1] <- 1
  SC <- Cs
  out <- FD.m.est(ai_vi = ai_vi,m = mm,q = q,nT = n,ai_vi_MLE = ai_vi_MLE)
  out <- as.vector(out)
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,length(q)*length(tau))
  m <- rep(mm,length(q)*length(tau))
  order <- rep(rep(q,each = length(mm)),length(tau))
  SC <- rep(SC,length(q)*length(tau))
  threshold <- rep(tau,each = length(q)*length(mm))
  method[SC == refC] = "Observed"
  data.frame(m = m,Method = method,Order.q = order,
             qFD = out,SC=SC, Tau = threshold)
}


FD_mle <- function(ai_vi, q){
  v_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    a <- ai_vi$ai[,i]
    v = ai_vi$vi[,i]
    sapply(q,function(qq){
      if(qq==1){
        exp(sum(-v*a/v_bar*log(a/v_bar)))
      }else{
        (sum(v*(a/v_bar)^qq))^(1 / (1-qq))
      }
    })
  })
  
  matrix(out,nrow = length(q),ncol = ncol(ai_vi$ai))
}


FDtable_mle <- function(datalist, dij, tau, q, datatype, nboot = 30, conf = 0.95){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(datatype=='abundance'){
    out <- lapply(datalist, function(x){
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      
      # emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      dmin <- min(dij[x>0,x>0][lower.tri(dij[x>0,x>0])])
      if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
        TDq = Diversity_profile_MLE(x, q)
        data_aivi.trun = list('ai' = data_aivi$ai[,tau > dmin], 'vi' = data_aivi$vi[,tau > dmin])
        emp <- c(rep(TDq, sum(tau <= dmin)),
                 FD_mle(ai_vi = data_aivi.trun,q = q) %>% as.numeric())
      } else if (sum(tau > dmin)!= 0) {
        emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      } else {
        TDq = Diversity_profile_MLE(x, q)
        emp <- rep(TDq, length(tau))
      }
      
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        ses <- sapply(1:nboot, function(B){
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          
          if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
            TDqb = Diversity_profile_MLE(Boot.X[,B], q)
            data_aivi.trun.b = list('ai' = Boot_aivi$ai[,tau > dmin], 'vi' = Boot_aivi$vi[,tau > dmin])
            c(rep(TDqb, sum(tau <= dmin)), 
              FD_mle(ai_vi = data_aivi.trun.b,q = q) %>% as.numeric())
          } else if (sum(tau > dmin)!= 0) {
            FD_mle(ai_vi = Boot_aivi,q = q) %>% as.numeric()
          } else {
            TDqb = Diversity_profile_MLE(Boot.X[,B], q)
            rep(TDqb, length(tau))
          }
        })
        
        if (length(q) == 1 & length(tau) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,ses,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype == 'incidence_freq'){
    out <- lapply(datalist, function(x){
      nT = x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      
      # emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      dmin <- min(dij[x[-1]>0,x[-1]>0][lower.tri(dij[x[-1]>0,x[-1]>0])])
      if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
        TDq = Diversity_profile_MLE.inc(x, q) %>% as.numeric()
        data_aivi.trun = list('ai' = data_aivi$ai[,tau > dmin], 'vi' = data_aivi$vi[,tau > dmin])
        emp <- c(rep(TDq, sum(tau <= dmin)),
                 FD_mle(ai_vi = data_aivi.trun,q = q) %>% as.numeric())
      } else if (sum(tau > dmin)!= 0) {
        emp <- FD_mle(ai_vi = data_aivi,q = q) %>% as.numeric()
      } else {
        TDq = Diversity_profile_MLE.inc(x, q) %>% as.numeric()
        emp <- rep(TDq, length(tau))
      }
      
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X,dij = dij_boot,tau = tau,datatype = datatype)
          
          # FD_mle(ai_vi = Boot_aivi,q = q) %>% as.numeric()
          if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
            TDqb = Diversity_profile_MLE.inc(Boot.X, q) %>% as.numeric()
            data_aivi.trun.b = list('ai' = Boot_aivi$ai[,tau > dmin], 'vi' = Boot_aivi$vi[,tau > dmin])
            c(rep(TDqb, sum(tau <= dmin)),
              FD_mle(ai_vi = data_aivi.trun.b,q = q) %>% as.numeric())
          } else if (sum(tau > dmin)!= 0) {
            FD_mle(ai_vi = Boot_aivi,q = q) %>% as.numeric()
          } else {
            TDqb = Diversity_profile_MLE.inc(Boot.X, q) %>% as.numeric()
            rep(TDqb, length(tau))
          }
          
        })
        
        if (length(q) == 1 & length(tau) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,ses,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
    ### to be added
  }
  sites_tmp <- rep(sites,each = length(q)*length(tau))
  tau_tmp <- rep(rep(tau,each = length(q)),length(sites))
  Output <- data.frame(Assemblage = sites_tmp, Order.q = rep(q,length(tau)*length(sites)), qFD = out[,1],
                       s.e. = out[,2], qFD.LCL = out[,3], qFD.UCL = out[,4],
                       Method='Observed', Tau = tau_tmp)
  Output
}


FD_est = function(ai_vi, q, nT, ai_vi_MLE){ # ai_vi is array containing two elements: ai and vi
  V_bar <- sum(ai_vi$ai[,1]*ai_vi$vi[,1])/nT
  
  Sub <- function(q,FD_obs,nT,f1,f2,h1,h2,A,av,avtab,deltas){
    if(q==0){
      ans <- FD_obs+FDq0(nT,f1,f2,h1,h2,A)
    }else if(q==1){
      h_est_2 <- FDq1_1(nT,h1,A)
      h_est_1 <- av %>% filter(ai<=(nT-1)) %>% mutate(diga = digamma(nT)-digamma(ai)) %>%
        apply(., 1, prod) %>% sum(.)/nT
      ans <- V_bar*exp((h_est_1+h_est_2)/V_bar)
    }else if(q==2){
      ans <- FDq2(as.matrix(avtab),nT)*V_bar^2
    }else{
      k <- 0:(nT-1)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(h1==0|A==1,0,(h1*((1-A)^(1-nT))/nT)*(round(A^(q-1)-sum(choose(q-1,k)*(A-1)^k),12)))
      ans <- ((a+b)/(V_bar^q))^(1/(1-q))
    }
    return(ans)
  }
  
  out <- sapply(1:ncol(ai_vi$ai), function(i){
    
    if ( !((sum(ai_vi$ai[,i] == 1) == 0) | (sum(ai_vi$ai[,i] == 2) == 0 & sum(ai_vi$ai[,i] == 1) == 1)) ) {
      
      av = tibble(ai = ai_vi$ai[,i], vi = ai_vi$vi[,i])
      FD_obs <- sum(av[,2])
      f1 <- sum(av[,1]==1); h1 <- ifelse(f1>0,sum(av[av[,1]==1,2]),0)
      f2 <- sum(av[,1]==2); h2 <- ifelse(f2>0,sum(av[av[,1]==2,2]),0)
      if(f2 > 0){
        A = 2*f2/((nT-1)*f1+2*f2)
      }else if(f2 == 0 & f1 > 0){
        A = 2/((nT-1)*(f1-1)+2)
      }else{
        A = 1
      }
      if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
        avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
        deltas <- sapply(0:(nT-1), function(k){
          del_tmp <- avtab[avtab[,1]<=(nT-k),,drop=FALSE]
          delta(del_tmp,k,nT)
        })
      }else{
        deltas <- 0
        avtab <- av %>% group_by(ai, vi) %>% summarise(n_group = n()) %>% as.matrix()
      }
      c(sapply(q, function(qq) Sub(qq,FD_obs,nT,f1,f2,h1,h2,A,av,avtab,deltas)))
      
    } else FD_mle(list('ai' = ai_vi_MLE$ai[, i, drop = FALSE], 'vi' = ai_vi_MLE$vi[, i, drop = FALSE]), q) %>% c
  }) 
  
  matrix(out,ncol = ncol(ai_vi$ai))
}


FDtable_est <- function(datalist, dij, tau, q, datatype, nboot = 30, conf = 0.95){#change final list name
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(datatype=="abundance"){
    out <- lapply(datalist, function(x){
      n=sum(x)
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      
      # est <- FD_est(ai_vi = data_aivi,q = q,nT = n)
      dmin <- min(dij[x>0,x>0][lower.tri(dij[x>0,x>0])])
      
      if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
        
        TDq = Diversity_profile(x, q)
        data_aivi.trun = list('ai' = data_aivi$ai[,tau > dmin], 'vi' = data_aivi$vi[,tau > dmin])
        data_aivi_MLE.trun = list('ai' = data_aivi_MLE$ai[,tau > dmin], 'vi' = data_aivi_MLE$vi[,tau > dmin])
        
        est <- c(rep(TDq, sum(tau <= dmin)),
                 FD_est(ai_vi = data_aivi.trun, q = q, nT = n, ai_vi_MLE = data_aivi_MLE.trun) %>% as.numeric())
        
      } else if (sum(tau > dmin)!= 0) {
        est <- FD_est(ai_vi = data_aivi, q = q, nT = n, ai_vi_MLE = data_aivi_MLE) %>% as.numeric()
        
      } else {
        TDq = Diversity_profile(x, q)
        est <- rep(TDq, length(tau))
      }
      
      
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, n, p_hat)
        
        ses <- sapply(1:nboot, function(B){
          dmin <- min(dij_boot[Boot.X[,B]>0,Boot.X[,B]>0][lower.tri(dij_boot[Boot.X[,B]>0,Boot.X[,B]>0])])
          
          Boot_aivi <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X[,B],dij = dij_boot,tau = tau,datatype = datatype)
          
          # FD_est(ai_vi = Boot_aivi,q = q,nT = n) %>% as.numeric()
          
          if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
            
            TDqb = Diversity_profile(Boot.X[,B], q)
            data_aivi.trun.b = list('ai' = Boot_aivi$ai[,tau > dmin], 'vi' = Boot_aivi$vi[,tau > dmin])
            data_aivi_MLE.trun.b = list('ai' = Boot_aivi_MLE$ai[,tau > dmin], 'vi' = Boot_aivi_MLE$vi[,tau > dmin])
            
            c(rep(TDqb, sum(tau <= dmin)),
              FD_est(ai_vi = data_aivi.trun.b, q = q, nT = n, ai_vi_MLE = data_aivi_MLE.trun.b) %>% as.numeric())
            
          } else if (sum(tau > dmin)!= 0) {
            FD_est(ai_vi = Boot_aivi,q = q,nT = n, ai_vi_MLE = Boot_aivi_MLE) %>% as.numeric()
            
          } else {
            TDqb = Diversity_profile(Boot.X[,B], q)
            rep(TDqb, length(tau))
          }
        })
        
        if (length(q) == 1 & length(tau) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      output <- cbind(est,ses,est-qtile*ses,est+qtile*ses)
      output[output[,2]<0,2] <- 0
      # list(estimates = output,info = est_info$info)
      output
    }) 
    
  }else if(datatype=="incidence_freq"){
    
    out <- lapply(datalist, function(x){
      nT=x[1]
      data_aivi <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype, integer = TRUE)
      data_aivi_MLE <- data_transform(data = x,dij = dij,tau = tau,datatype = datatype)
      
      # est_info <- FD_est(ai_vi = data_aivi,q = q,nT = nT)
      dmin <- min(dij[x[-1]>0,x[-1]>0][lower.tri(dij[x[-1]>0,x[-1]>0])])
      
      if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
        
        TDq = Diversity_profile.inc(x, q) %>% as.numeric()
        data_aivi.trun = list('ai' = data_aivi$ai[,tau > dmin], 'vi' = data_aivi$vi[,tau > dmin])
        data_aivi_MLE.trun = list('ai' = data_aivi_MLE$ai[,tau > dmin], 'vi' = data_aivi_MLE$vi[,tau > dmin])
        
        est <- c(rep(TDq, sum(tau <= dmin)),
                 FD_est(ai_vi = data_aivi.trun, q = q, nT = nT, ai_vi_MLE = data_aivi_MLE.trun) %>% as.numeric())
        
      } else if (sum(tau > dmin)!= 0) {
        est <- FD_est(ai_vi = data_aivi, q = q, nT = nT, ai_vi_MLE = data_aivi_MLE) %>% as.numeric()
        
      } else {
        TDq = Diversity_profile.inc(x, q) %>% as.numeric()
        est <- rep(TDq, length(tau))
      }
      
      if(nboot>1){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        ses <- sapply(1:nboot, function(B){
          Boot.X <- c(nT,rbinom(n = p_hat,size = nT,prob = p_hat))
          Boot_aivi <- data_transform(data = Boot.X, dij = dij_boot,tau = tau,datatype = datatype, integer = TRUE)
          Boot_aivi_MLE <- data_transform(data = Boot.X, dij = dij_boot,tau = tau,datatype = datatype)
          
          dmin <- min(dij_boot[Boot.X[-1]>0,Boot.X[-1]>0][lower.tri(dij_boot[Boot.X[-1]>0,Boot.X[-1]>0])])
          
          # FD_est(ai_vi = Boot_aivi,q = q,nT = nT) %>% as.numeric()
          if (sum(tau <= dmin) != 0 & sum(tau > dmin) != 0) {
            
            TDqb = Diversity_profile.inc(Boot.X, q) %>% as.numeric()
            data_aivi.trun.b = list('ai' = Boot_aivi$ai[,tau > dmin], 'vi' = Boot_aivi$vi[,tau > dmin])
            data_aivi_MLE.trun.b = list('ai' = Boot_aivi_MLE$ai[,tau > dmin], 'vi' = Boot_aivi_MLE$vi[,tau > dmin])
            
            c(rep(TDqb, sum(tau <= dmin)),
              FD_est(ai_vi = data_aivi.trun.b, q = q, nT = nT, ai_vi_MLE = data_aivi_MLE.trun.b) %>% as.numeric())
            
          } else if (sum(tau > dmin)!= 0) {
            FD_est(ai_vi = Boot_aivi, q = q, nT = nT, ai_vi_MLE = Boot_aivi_MLE) %>% as.numeric()
            
          } else {
            TDqb = Diversity_profile.inc(Boot.X, q) %>% as.numeric()
            rep(TDqb, length(tau))
          }
        })
        
        if (length(q) == 1 & length(tau) == 1) ses = matrix(ses, nrow = 1)
        ses = ses %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      output <- cbind(est,ses,est-qtile*ses,est+qtile*ses)
      output[output[,3]<0,3] <- 0
      # list(estimates = output,info = est_info$info)
      output
    })
  }
  # info <- lapply(out, function(x) x[[2]]) %>% 
  #   do.call(rbind,.) %>% as_tibble %>% 
  #   mutate(Assemblage = rep(names(datalist),each = length(tau)),
  #          tau = rep(tau,length(datalist))) %>% 
  #   select(Assemblage,tau,nT,S.obs,f1,f2,h1,h2)
  
  
  Estoutput <- out %>% do.call(rbind,.) #%>% mutate(Assemblage = rep(names(datalist),each = length(q)))
  sites_tmp <- rep(sites,each = length(q)*length(tau))
  tau_tmp <- rep(rep(tau,each = length(q)),length(sites))
  Estoutput <- data.frame(Assemblage = sites_tmp, Order.q = rep(q,length(tau)*length(sites)), 
                          qFD = Estoutput[,1],
                          s.e. = Estoutput[,2], qFD.LCL = Estoutput[,3], qFD.UCL = Estoutput[,4],
                          Method = "Asymptotic", Tau = tau_tmp)
  Estoutput$qFD.LCL[Estoutput$qFD.LCL<0] = 0
  # return(list(Estoutput = Estoutput, info = info))
  return(Estoutput)
}


AUCtable_iNextFD <- function(datalist, dij, q = c(0,1,2), datatype, tau=NULL,
                         nboot=0, conf=0.95, m, FDcut_number = NULL) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    # tau <- seq(0,1,length.out = FDcut_number)
    tau <- seq(1e-8,1,length.out = FDcut_number)
  }
  AUC <- iNextFD(datalist,dij,q,datatype,tau,nboot = 0,m = m) %>%
    group_by(Assemblage,Order.q,m) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC),Method = unique(Method)) %>% ungroup %>% 
    mutate(qFD = (AUC_L+AUC_R)/2) %>% select(Assemblage,Order.q,m,qFD,SC,Method)
  if(datatype == 'abundance'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Assemblage_ <- rep(sites[[i]],length(q)*length(m[[i]]))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        m_boot <- lapply(1:nboot, function(b) m[[i]])
        ses <- iNextFD(Boot.X,dij_boot,q,datatype,tau,nboot = 0,m = m_boot) %>% 
          group_by(Assemblage,Order.q,m) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q,m) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(Assemblage = Assemblage_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(Assemblage,Order.q,m) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }else if (datatype == 'incidence_freq'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Assemblage_ <- rep(sites[[i]],length(q)*length(m[[i]]))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>% 
          split(., col(.))
        m_boot <- lapply(1:nboot, function(b) m[[i]])
        ses <- iNextFD(Boot.X,dij_boot,q,datatype,tau,nboot = 0,m = m_boot) %>% 
          group_by(Assemblage,Order.q,m) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q,m) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(Assemblage = Assemblage_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(Assemblage,Order.q,m) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }
  
  AUC <- left_join(x = AUC, y = ses, by = c('Assemblage','Order.q','m')) %>% mutate(
    s.e. = AUC_se, qFD.LCL = qFD - AUC_se * qtile, qFD.UCL = qFD + AUC_se * qtile,
    SC.s.e. = SC_se, SC.LCL = SC - SC_se * qtile, SC.UCL = SC + SC_se * qtile) %>% 
    select(Assemblage,Order.q,m,Method,qFD,s.e.,qFD.LCL,qFD.UCL,SC,SC.s.e.,SC.LCL,SC.UCL)
  AUC$qFD.LCL[AUC$qFD.LCL<0] <- 0
  AUC$SC.LCL[AUC$SC.LCL<0] <- 0
  AUC$SC.UCL[AUC$SC.UCL>1] <- 1
  AUC
}


AUCtable_invFD <- function(datalist, dij, q = c(0,1,2), datatype, level, nboot = 0, conf = 0.95, tau = NULL, FDcut_number = NULL){
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  if(is.null(tau)){
    # tau <- seq(0,1,length.out = FDcut_number)
    tau <- seq(1e-8,1,length.out = FDcut_number)
  }
  AUC <- invChatFD(datalist,dij,q,datatype,level,nboot = 0,tau = tau) %>%
    group_by(Assemblage,Order.q,SC) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC),m = mean(m),Method = unique(Method)) %>% 
    ungroup %>% mutate(qFD = (AUC_L+AUC_R)/2) %>% select(Assemblage,Order.q,m,Method,qFD,SC)
  if(datatype == 'abundance'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Assemblage_ <- rep(sites[[i]],length(q)*length(level))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- invChatFD(Boot.X,dij_boot,q,datatype,level,nboot = 0,tau = tau) %>% 
          group_by(Assemblage,Order.q,SC) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q,SC) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(Assemblage = Assemblage_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(Assemblage,Order.q,SC) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }else if(datatype == 'incidence_freq'){
    if(nboot>1){
      ses <- lapply(1:length(datalist),function(i){
        Assemblage_ <- rep(sites[[i]],length(q)*length(level))
        x <- datalist[[i]]
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>% 
          split(., col(.))
        ses <- invChatFD(Boot.X,dij_boot,q,datatype,level,nboot = 0,tau = tau) %>% 
          group_by(Assemblage,Order.q,SC) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau)),SC = mean(SC)) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q,SC) %>% 
          summarise(AUC_se = sd(AUC),SC_se = sd(SC)) %>% 
          ungroup %>% mutate(Assemblage = Assemblage_)
      }) %>% do.call(rbind,.) 
    }else{
      ses <- AUC %>% select(Assemblage,Order.q,SC) %>% mutate(AUC_se = NA, SC_se = NA)
    }
  }
  
  AUC <- left_join(x = AUC, y = ses, by = c('Assemblage','Order.q','SC')) %>% mutate(
    s.e. = AUC_se, qFD.LCL = qFD - AUC_se * qtile, qFD.UCL = qFD + AUC_se * qtile,
    SC.s.e. = SC_se, SC.LCL = SC - SC_se * qtile, SC.UCL = SC + SC_se * qtile) %>% 
    select(Assemblage, Order.q, SC, m, Method, qFD, s.e., qFD.LCL, qFD.UCL)
  AUC$qFD.LCL[AUC$qFD.LCL<0] <- 0
  AUC
}


AUCtable_mle <- function(datalist, dij, q = c(0,1,2), tau=NULL, datatype,
                         nboot=0, conf=0.95, FDcut_number = NULL) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    # tau <- seq(0,1,length.out = FDcut_number)
    tau <- seq(1e-8,1,length.out = FDcut_number)
  }
  #q_int <- c(0, 1, 2)
  
  AUC <- FDtable_mle(datalist,dij,tau,q,datatype,nboot = 0) %>%
    group_by(Assemblage,Order.q) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
    mutate(qFD = (AUC_L+AUC_R)/2) %>% select(Assemblage,Order.q,qFD)
  if(datatype=='abundance'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- FDtable_mle(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Assemblage,Order.q) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(sites,each = length(q)))
    }else{
      ses <- data.frame(Order.q = rep(q,length(datalist)), se = NA, Assemblage = rep(sites,each = length(q)))
    }
  }else if(datatype=='incidence_freq'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>%
          split(., col(.))
        ses <- FDtable_mle(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Assemblage,Order.q) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(sites,each = length(q)))
    }else{
      ses <- data.frame(Order.q = rep(q,length(datalist)), se = NA, Assemblage = rep(sites,each = length(q)))
    }
  }
  
  AUC <- left_join(x = AUC, y = ses, by = c('Assemblage','Order.q')) %>% mutate(s.e. = se, 
                                                                                qFD.LCL = qFD - se * qtile,
                                                                                qFD.UCL = qFD + se * qtile,
                                                                                Method = "Observed") %>% 
    select(-se)
  AUC$qFD.LCL[AUC$qFD.LCL<0] <- 0
  AUC = AUC %>% select(c('Assemblage', 'Order.q', 'qFD', 's.e.', 'qFD.LCL', 'qFD.UCL', 'Method'))
  AUC = data.frame(AUC)
  
  AUC
}


AUCtable_est <- function(datalist, dij, q = c(0,1,2), tau=NULL, datatype,
                         nboot=0, conf=0.95, FDcut_number = NULL) {
  qtile <- qnorm(1-(1-conf)/2)
  sites <- names(datalist)
  # dmin <- min(dij[dij>0])
  # dmax <- max(dij)
  # if(is.null(tau)){
  #   tau <- seq(dmin,dmax,length.out = knots)
  # }
  if(is.null(tau)){
    # tau <- seq(0,1,length.out = FDcut_number)
    tau <- seq(1e-8,1,length.out = FDcut_number)
  }
  #q_int <- c(0, 1, 2)
  AUC <- FDtable_est(datalist,dij,tau,q,datatype,nboot = 0) %>%
    group_by(Assemblage,Order.q) %>% 
    summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
              AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
    mutate(qFD = (AUC_L+AUC_R)/2) %>% select(Assemblage,Order.q,qFD)
  if(datatype=='abundance'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X = rmultinom(nboot, sum(x), p_hat) %>% split(., col(.))
        ses <- FDtable_est(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Order.q,Assemblage) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(sites,each = length(q)))
    }else{
      ses <- data.frame(Order.q = rep(q,length(datalist)), se = NA, Assemblage = rep(sites,each = length(q)))
    }
  }else if(datatype=='incidence_freq'){
    if(nboot>1){
      ses <- lapply(datalist,function(x){
        BT <- EstiBootComm.Func(data = x,distance = dij,datatype = datatype)
        p_hat = BT[[1]]
        dij_boot = BT[[2]]
        Boot.X <- sapply(1:nboot,function(b) c(x[1],rbinom(n = p_hat,size = x[1],prob = p_hat))) %>%
          split(., col(.))
        ses <- FDtable_est(Boot.X,dij_boot,tau,q,datatype,nboot = 0) %>% 
          group_by(Order.q,Assemblage) %>% 
          summarise(AUC_L = sum(qFD[seq_along(qFD[-1])]*diff(tau)),
                    AUC_R = sum(qFD[-1]*diff(tau))) %>% ungroup %>% 
          mutate(AUC = (AUC_L+AUC_R)/2) %>% group_by(Order.q) %>% summarise(se = sd(AUC)) %>% 
          ungroup
      }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(sites,each = length(q)))
    }else{
      ses <- data.frame(Order.q = rep(q,length(datalist)), se = NA, Assemblage = rep(sites,each = length(q)))
    }
  }
  AUC <- left_join(x = AUC, y = ses, by = c('Assemblage','Order.q')) %>% mutate(s.e. = se, 
                                                                                qFD.LCL = qFD - se * qtile,
                                                                                qFD.UCL = qFD + se * qtile,
                                                                                Method = "Asymptotic") %>% 
    select(-se)
  AUC$qFD.LCL[AUC$qFD.LCL<0] <- 0
  AUC = AUC %>% select(c('Assemblage', 'Order.q', 'qFD', 's.e.', 'qFD.LCL', 'qFD.UCL', 'Method'))
  AUC = data.frame(AUC)
  
  return(AUC)
}


EstiBootComm.Func = function(data, distance, datatype){
  if (datatype=="incidence_freq") {
    n <- data[1]
    data <- data[-1]
    u=sum(data)
  } else if (datatype=="abundance") {
    n = sum(data)
    data=data
  }
  distance = as.matrix(distance)
  dij = distance[data!=0, data!=0]
  
  X = data[data>0]
  f1 <- sum(X == 1) ; f2 <- sum(X == 2)
  f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
  if (datatype=="abundance") {
    C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
    W <- (1 - C1)/sum(X/n*(1-X/n)^n)
    Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)
  } else if (datatype=="incidence_freq") {
    C1 = ifelse(f2>0, 1-f1/u*(n-1)*f1/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/u/((n-1)*(f1-1)+2))
    W <- (1 - C1)/sum(X/u*(1-X/n)^n)
    Prob.hat.Unse <- rep(u/n*(1-C1)/f0.hat, f0.hat)
  }
  
  Prob.hat <- X/n*(1-W*(1-X/n)^n)
  Prob <- c(Prob.hat, Prob.hat.Unse)
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  if (datatype=="abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype=="incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat==0) {
    d=dij
  } else if (f0.hat==1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(list("pi" = Prob,"dij" = d))
}


