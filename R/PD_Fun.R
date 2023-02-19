PhD.m.est = function(ai,Lis, m, q, nt, reft, cal){
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  if(sum(m>nt)>0){
    #Extrapolation
    RPD_m <- RPD(ai,Lis,nt,nt-1,q)
    obs <- RPD(ai, Lis, nt,nt, q)
    EPD = function(m,obs,asy){
      m = m-nt
      out <- sapply(1:ncol(Lis), function(i){
        asy_i <- asy[,i];obs_i <- obs[,i];RPD_m_i <- RPD_m[,i]
        Li <- Lis[,i];t_bar <- t_bars[i]
        asy_i <- sapply(1:length(q), function(j){
          max(asy_i[j],obs_i[j])
        })
        beta <- rep(0,length(q))
        beta0plus <- which(asy_i != obs_i)
        beta[beta0plus] <-(obs_i[beta0plus]-RPD_m_i[beta0plus])/(asy_i[beta0plus]-RPD_m_i[beta0plus])
        outq <- sapply(1:length(q), function(i){
          if( q[i]!=2 ) {
            obs_i[i]+(asy_i[i]-obs_i[i])*(1-(1-beta[i])^m)
          }else if( q[i] == 2 ){
            1/sum( (Li/(t_bar)^2)*((1/(nt+m))*(ai/nt)+((nt+m-1)/(nt+m))*(ai*(ai-1)/(nt*(nt-1)))) )
          }
        })
        outq
      })
      return(out)
    }
    # obs <- PD.qprofile(aL = aL, q = Q, cal="PD" ,datatype = datatype , nforboot = nforboot, splunits = splunits)
    #asymptotic value
    asy <- matrix(PhD.q.est(ai = ai,Lis = Lis,q = q, nt = nt, reft = reft, cal = 'PD'),nrow = length(q),ncol = length(t_bars))
  }else if (sum(m==nt)>0){
    obs <- RPD(ai, Lis, nt,nt, q)
  }
  
  if (sum(m < nt) != 0) {
    int.m = sort(unique(c(floor(m[m<nt]), ceiling(m[m<nt]))))
    mRPD = lapply(int.m, function(k) RPD(ai = ai,Lis = Lis,n = nt,m = k,q = q))
    names(mRPD) = int.m
  }
  
  if (cal == 'PD'){
    out <- sapply(m, function(mm){
      if(mm<nt){
        if(mm == round(mm)) { ans <- mRPD[names(mRPD) == mm][[1]]
        } else { ans <- (ceiling(mm) - mm)*mRPD[names(mRPD) == floor(mm)][[1]] + (mm - floor(mm))*mRPD[names(mRPD) == ceiling(mm)][[1]] }
      }else if(mm==nt){
        ans <- obs
      }else if(mm==Inf){
        ans <- asy
      }else{
        ans <- EPD(m = mm,obs = obs,asy = asy)
      }
      return(as.numeric(ans))
    })
  } else if (cal == 'meanPD') {
    out <- sapply(m, function(mm){
      if(mm<nt){
        if(mm == round(mm)) { ans <- mRPD[names(mRPD) == mm][[1]]
        } else { ans <- (ceiling(mm) - mm)*mRPD[names(mRPD) == floor(mm)][[1]] + (mm - floor(mm))*mRPD[names(mRPD) == ceiling(mm)][[1]] }
      }else if(mm==nt){
        ans <- obs
      }else if(mm==Inf){
        ans <- asy
      }else{
        ans <- EPD(m = mm,obs = obs,asy = asy)
      }
      if (inherits(ans, c('numeric', 'integer'))) ans = matrix(ans, ncol = length(reft))
      ans <- sapply(1:length(reft), function(i){
        ans[,i]/reft[i]
      })
      as.numeric(ans)
    })
  }
  out <- matrix(out,ncol = length(m))
  return(out)
}

inextPD = function(datalist, datatype, phylotr, q,reft, m, cal, nboot, conf=0.95, unconditional_var=TRUE){
  # m is a list
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      #====conditional on m====
      qPDm <- PhD.m.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,m = m[[i]],
                        q = q,nt = n,reft = reft,cal = cal) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- invChatPD_abu(x = x,ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                 q = q,Cs = goalSC,n = n,reft = reft,cal = cal)
        qPD_unc$Method[round(qPD_unc$m) == n] = "Observed"
      }
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT)
        Li_b <- Boots$Li
        refinfo <- colnames(Li_b)
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        colnames(Li_b) <- refinfo
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        if(unconditional_var){
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,reft = reft,cal = cal) %>% as.numeric()
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]])
            qPD_unc_b <- invChatPD_abu(x = ai_B[isn0&tgroup_B=="Tip"],
                                       ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                       q = q,Cs = goalSC,n = n,reft = reft,cal = cal)$qPD
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b,qPD_unc_b))
          }) %>% apply(., 1, sd)
        }else{
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,reft = reft,cal = cal) %>% as.numeric()
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]])
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b))
          }) %>% apply(., 1, sd)
        }
      }else{
        if(unconditional_var){
          ses <- rep(NA,length(c(qPDm,covm,qPD_unc$qPD)))
        }else{
          ses <- rep(NA,length(c(qPDm,covm)))
        }
      }
      
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm)+1):(length(qPDm)+length(covm))]
      m_ <- rep(m[[i]],each = length(q)*length(reft))
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,each = length(q)*length(reft))
      orderq <- rep(q,length(reft)*length(m[[i]]))
      SC_ <- rep(covm,each = length(q)*length(reft))
      SC.se <- rep(ses_cov,each = length(q)*length(reft))
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ <- rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- data.frame(Assemblage = nms[i], m=m_,Method=method,Order.q=orderq,
                          qPD=qPDm,s.e. = ses_pd,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                          SC=SC_,SC.s.e.=SC.se,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                          Reftime = reft_,
                          Type=cal) %>%
        arrange(Reftime,Order.q,m)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    s.e.=ses_pd_unc, Type=cal,
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','SC','m', 'Method', 'Order.q', 'qPD', 's.e.', 'qPD.LCL','qPD.UCL','Reftime',
                        'Type'), names(out_C), nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reftime,Order.q,m)
        out_C$qPD.LCL[out_C$qPD.LCL<0] <- 0
      }else{
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,colSums(.)>0]
      n <- ncol(x)
      #====conditional on m====
      qPDm <- PhD.m.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,m = m[[i]],
                        q = q,nt = n,reft = reft,cal = cal) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- invChatPD_inc(x = rowSums(x),ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                 q = q,Cs = goalSC,n = n,reft = reft,cal = cal)
        qPD_unc$Method[round(qPD_unc$nt) == n] = "Observed"
      }
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        refinfo <- colnames(Li_b)
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        colnames(Li_b) <- refinfo
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        if(unconditional_var){
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,reft = reft,cal = cal) %>% as.numeric()
            covm_b = Coverage(c(n, ai_B[isn0&tgroup_B=="Tip"]), "incidence_freq", m[[i]])
            qPD_unc_b <- invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],
                                       ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                       q = q,Cs = goalSC,n = n,reft = reft,cal = cal)$qPD
            # btime <- Sys.time()
            # print(paste0("Est boot sample",B,": ",btime-atime))
            return(c(qPDm_b,covm_b,qPD_unc_b))
          }) %>% apply(., 1, sd)
        }else{
          ses <- sapply(1:nboot, function(B){
            # atime <- Sys.time()
            ai_B <- Boots$boot_data[,B]
            isn0 <- ai_B>0
            qPDm_b <-  PhD.m.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                 m=m[[i]],q=q,nt = n,reft = reft,cal = cal) %>% as.numeric()
            covm_b = Coverage(c(n, ai_B[isn0&tgroup_B=="Tip"]), "incidence_freq", m[[i]])
            return(c(qPDm_b,covm_b))
          }) %>% apply(., 1, sd)
        }
      }else{
        if(unconditional_var){
          ses <- rep(NA,length(c(qPDm,covm,qPD_unc$qPD)))
        }else{
          ses <- rep(NA,length(c(qPDm,covm)))
        }
      }
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm)+1):(length(qPDm)+length(covm))]
      m_ <- rep(m[[i]],each = length(q)*length(reft))
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      method <- rep(method,each = length(q)*length(reft))
      orderq <- rep(q,length(reft)*length(m[[i]]))
      SC_ <- rep(covm,each = length(q)*length(reft))
      SC.se <- rep(ses_cov,each = length(q)*length(reft))
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ = rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- data.frame(Assemblage = nms[i], nt=m_,Method=method,Order.q=orderq,
                          qPD=qPDm,s.e.=ses_pd,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                          SC=SC_,SC.s.e.=SC.se,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                          Reftime = reft_,
                          Type=cal) %>%
        arrange(Reftime,Order.q,nt)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    s.e. = ses_pd_unc, Type=cal,
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','SC','nt', 'Method', 'Order.q', 'qPD', 's.e.', 'qPD.LCL','qPD.UCL','Reftime',
                        'Type'), names(out_C), nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reftime,Order.q,nt)
        out_C$qPD.LCL[out_C$qPD.LCL<0] <- 0
      }else{
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }
  if(unconditional_var){
    ans <- list(size_based = lapply(Estoutput, function(x) x$size_based) %>% do.call(rbind,.) %>% as_tibble,
                coverage_based = lapply(Estoutput, function(x) x$coverage_based) %>% do.call(rbind,.) %>% as_tibble)
  }else{
    ans <- list(size_based = do.call(rbind,lapply(Estoutput, function(x) x$size_based)))
  }
  return(ans)
}


invChatPD <- function(datalist, datatype,phylotr, q, reft, cal,level, nboot, conf){
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=='abundance'){
    out <- lapply(datalist,function(x_){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x_,'abundance',refT = reft)
      x_ <- x_[x_>0]
      n <- sum(x_)
      est <- invChatPD_abu(x = x_,ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                           q = q,Cs = level, n = n, reft = reft, cal = cal)
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        refinfo <- colnames(Li_b)
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        colnames(Li_b) <- refinfo
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x_)+f0),rep("Inode",nrow(Li_b)-length(x_)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          # atime <- Sys.time()
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          # isn0 <- as.vector(aL_table_b[,1]>0)
          # Li_b_tmp <- Li_b[isn0,]
          # aL_table_b <- aL_table_b[isn0,]
          est_b <- invChatPD_abu(x = ai_B[isn0&tgroup_B=="Tip"],ai = ai_B[isn0],
                                 Lis = Li_b[isn0,,drop=F],q = q,Cs = level,
                                 n = n, reft = reft, cal = cal)$qPD
          
          return(est_b)
        }) %>% matrix(.,nrow = length(q)*length(reft)*length(level)) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(s.e.=ses,qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
    }) %>% do.call(rbind,.)
  }else if(datatype=='incidence_raw'){
    out <- lapply(datalist,function(x_){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x_,'incidence_raw',refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      x_ <- x_[rowSums(x_)>0,colSums(x_)>0]
      n <- ncol(x_)
      est <- invChatPD_inc(x = rowSums(x_),ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                           q = q,Cs = level, n = n, reft = reft, cal = cal)
      if(nboot>1){
        Boots <- Boots.one(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        refinfo <- colnames(Li_b)
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        colnames(Li_b) <- refinfo
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x_)+f0),rep("Inode",nrow(Li_b)-nrow(x_)-f0))
        #aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          est_b <- invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],ai = ai_B[isn0],
                                 Lis = Li_b[isn0,,drop=F],q = q,Cs = level,
                                 n = n, reft = reft, cal = cal)$qPD
          return(est_b)
        }) %>% matrix(.,nrow = length(q)*length(reft)*length(level)) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(s.e.=ses,qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
    }) %>% do.call(rbind,.)
  }
  Assemblage = rep(names(datalist), each = length(q)*length(reft)*length(level))
  out <- out %>% mutate(Assemblage = Assemblage, Type=cal)
  if(datatype=='abundance'){
    out <- out %>% select(Assemblage,SC,m,Method,Order.q,qPD,s.e.,qPD.LCL,qPD.UCL,
                          Reftime,Type)
  }else if(datatype=='incidence_raw'){
    out <- out %>% select(Assemblage,SC,nt,Method,Order.q,qPD,s.e.,qPD.LCL,qPD.UCL,
                          Reftime,Type)
  }
  out$qPD.LCL[out$qPD.LCL<0] <- 0
  rownames(out) <- NULL
  out
}


invChatPD_abu <- function(x,ai,Lis, q, Cs, n, reft, cal){
  #x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
  refC = Coverage(x, 'abundance', n)
  f <- function(m, C) abs(Coverage(x, 'abundance', m) - C)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
    }else if (refC <= cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if(f1 > 1 & f2 == 0) {
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
  mm[which(round(mm) - n <= 1)] = round(mm[which(round(mm) - n <= 1)]) 
  SC <- Cs
  out <- as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,reft=reft,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reftime = rep(rep(reft,each = length(q)),length(Cs))
  data.frame(m = m,Method = method,Order.q = order,
             qPD = out,SC=SC,Reftime = Reftime)
}


invChatPD_inc <- function(x,ai,Lis, q, Cs, n, reft, cal){
  refC = Coverage(c(n, x), 'incidence_freq', n)
  f <- function(m, C) abs(Coverage(c(n, x), 'incidence_freq', m) - C)
  mm <- sapply(Cs, function(cvrg){
    if (refC > cvrg) {
      opt <- optimize(f, C = cvrg, lower = 0, upper = n)
      mm <- opt$minimum
    }else if (refC <= cvrg) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }else if(f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }else if (f1 == 0 & f2 > 0) {
        A <- 0
      }else if(f1 == 1 & f2 == 0) {
        A <- 0
      }else if(f1 == 0 & f2 == 0) {
        A <- 0
      }
      mm <- ifelse(A==0,0,(log(U/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
    }
    mm
  })
  mm[mm < 1] <- 1
  mm[which(round(mm) - n <= 1)] = round(mm[which(round(mm) - n <= 1)]) 
  SC <- Cs
  out <-  as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,reft = reft,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reftime = rep(rep(reft,each = length(q)),length(Cs))
  data.frame(nt = m,Method = method,Order.q = order,
             qPD = out,SC=SC, Reftime = Reftime)
}



EmpPD <- function(datalist, datatype, phylotr, q, reft, cal, nboot, conf){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      # emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,reft = reft,cal = cal,nt = n) %>% c()
      first.Inode = subset(aL$treeNabu, tgroup == 'Tip')$branch.length %>% min
      if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
        TDq = Diversity_profile_MLE(x, q)
        emp <- cbind(sapply(reft[first.Inode >= reft], function(i) TDq*ifelse(cal == 'meanPD', 1, i)),
                     PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT[,first.Inode < reft,drop = F], q=q,reft = reft[first.Inode < reft],cal = cal,nt = n) %>% t) %>% t %>% c
      } else if (sum(first.Inode < reft) != 0) {
        emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,reft = reft,cal = cal,nt = n) %>% c()
      } else {
        TDq = Diversity_profile(x, q)
        emp <- sapply(reft, function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% t %>% c}
      
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          
          # out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal, nt = n) %>% c()
          if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
            TDqb = Diversity_profile_MLE(Boots$boot_data[1:(length(x)+f0),B], q)
            out_b <- cbind(sapply(reft[first.Inode >= reft], function(i) TDqb*ifelse(cal == 'meanPD', 1, i)),
                          PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,first.Inode < reft,drop=F],q = q,reft = reft[first.Inode < reft],cal = cal,nt = n) %>% t) %>% t %>% c
          } else if (sum(first.Inode < reft) != 0) {
            out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal, nt = n) %>% c()
          } else {
            TDqb = Diversity_profile_MLE(Boots$boot_data[1:(length(x)+f0),B], q)
            out_b <- sapply(reft, function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% t %>% c}
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,ses,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      # emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q,reft = reft,cal = cal,nt = n) %>% c()
      first.Inode = subset(aL$treeNabu, tgroup == 'Tip')$branch.length %>% min
      if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
        TDq = Diversity_profile_MLE.inc(c(ncol(x), rowSums(x)), q)
        emp <- cbind(sapply(reft[first.Inode >= reft], function(i) TDq*ifelse(cal == 'meanPD', 1, i)),
                     PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT[,first.Inode < reft,drop = F], q=q,reft = reft[first.Inode < reft],cal = cal,nt = n) %>% t) %>% t %>% c
      } else if (sum(first.Inode < reft) != 0) {
        emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,reft = reft,cal = cal,nt = n) %>% c()
      } else {
        TDq = Diversity_profile_MLE.inc(c(ncol(x), rowSums(x)), q)
        emp <- sapply(reft, function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% t %>% c}
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft,
                           BLs = aL$BLbyT,splunits = n)
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          
          # out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal,nt = n) %>% c()
          if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
            TDqb = Diversity_profile_MLE.inc(c(ncol(x), Boots$boot_data[1:(nrow(x)+f0),B]), q)
            out_b <- cbind(sapply(reft[first.Inode >= reft], function(i) TDqb*ifelse(cal == 'meanPD', 1, i)),
                          PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,first.Inode < reft,drop=F],q = q,reft = reft[first.Inode < reft],cal = cal,nt = n) %>% t) %>% t %>% c
          } else if (sum(first.Inode < reft) != 0) {
            out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal, nt = n) %>% c()
          } else {
            TDqb = Diversity_profile_MLE.inc(c(ncol(x), Boots$boot_data[1:(nrow(x)+f0),B]), q)
            out_b <- sapply(reft, function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% t %>% c}
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,ses,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  Output <- tibble(Order.q = rep(rep(q, each=length(reft)),length(datalist)),
                   qPD = out[,1], s.e. = out[,2],qPD.LCL = out[,3], qPD.UCL = out[,4],
                   Assemblage = rep(nms, each=length(reft)*length(q)),
                   Method='Empirical',
                   Reftime = rep(reft,length(q)*length(datalist)),
                   Type=cal) %>%
    arrange(Reftime)
  return(Output)
}


PD.Tprofile=function(ai,Lis, q, reft, cal, nt) {
  isn0 <- ai>0
  ai <- ai[isn0]
  Lis <- Lis[isn0,,drop=F]
  #q can be a vector
  t_bars <- as.numeric(ai %*% Lis/nt)
  
  pAbun <- ai/nt
  
  out <- sapply(1:length(t_bars), function(i) {
    sapply(q, function(j){
      if(j==1) as.numeric(-(pAbun/t_bars[i]*log(pAbun/t_bars[i])) %*% Lis[,i]) %>% exp()
      else as.numeric((pAbun/t_bars[i])^j %*% Lis[,i]) %>% .^(1/(1-j))
    }) %>% matrix(.,ncol = length(q))
  }) %>% t()
  if(cal=="meanPD"){
    out <- sapply(1:length(reft), function(i){
      out[i,]/reft[i]
    }) %>% matrix(.,ncol = length(reft)) %>% t()
  }
  out
}


asymPD <- function(datalist, datatype, phylotr, q,reft, cal,nboot, conf){#change final list name
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  tau_l <- length(reft)
  if(datatype=="abundance"){
    Estoutput <- lapply(datalist,function(x){
      x <- x[x>0]
      n <- sum(x)
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      
      # est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
      first.Inode = subset(aL$treeNabu, tgroup == 'Tip')$branch.length %>% min
      if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
        TDq = Diversity_profile(x, q)
        est <- c(sapply(reft[first.Inode >= reft], function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% as.vector(),
                 PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT[,first.Inode < reft,drop = F],q = q,nt = n,reft = reft[first.Inode < reft],cal = cal))
      } else if (sum(first.Inode < reft) != 0) {
        est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
      } else {
        TDq = Diversity_profile(x, q)
        est <- sapply(reft, function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% as.vector()}
      
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        # tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        # aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          
          # outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
            TDqb = Diversity_profile(Boots$boot_data[1:(length(x)+f0),B], q)
            outb <- c(sapply(reft[first.Inode >= reft], function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% as.vector(),
                      PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,first.Inode < reft,drop=F],q = q,nt = n,reft = reft[first.Inode < reft],cal = cal))
          } else if (sum(first.Inode < reft) != 0) {
            outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          } else {
            TDqb = Diversity_profile(Boots$boot_data[1:(length(x)+f0),B], q)
            outb <- sapply(reft, function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% as.vector()}
          
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est, s.e. = ses,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reftime = rep(reft,each = length(q)))
      est
    })
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(datalist,function(x){
      x <- x[rowSums(x)>0,colSums(x)>0]
      n <- ncol(x)
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x,datatype,refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      # est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
      first.Inode = subset(aL$treeNabu, tgroup == 'Tip')$branch.length %>% min
      if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
        TDq = Diversity_profile.inc(c(ncol(x), rowSums(x)), q)
        est <- c(sapply(reft[first.Inode >= reft], function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% as.vector(),
                 PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT[,first.Inode < reft,drop = F],q = q,nt = n,reft = reft[first.Inode < reft],cal = cal))
      } else if (sum(first.Inode < reft) != 0) {
        est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
      } else {
        TDq = Diversity_profile.inc(c(ncol(x), rowSums(x)), q)
        est <- sapply(reft, function(i) TDq*ifelse(cal == 'meanPD', 1, i)) %>% as.vector()}
      
      
      if(nboot!=0){
        Boots <- Boots.one(phylo = phylotr,aL = aL$treeNabu,datatype = datatype,nboot = nboot,
                           splunits = n,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        Li_b <- sapply(1:length(reft),function(l){
          tmp <- Li_b[,l]
          tmp[tmp>reft[l]] <- reft[l]
          tmp
        })
        f0 <- Boots$f0
        # tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        # aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          
          # outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          if (sum(first.Inode >= reft) != 0 & sum(first.Inode < reft) != 0) {
            TDqb = Diversity_profile.inc(c(ncol(x), Boots$boot_data[1:(nrow(x)+f0),B]), q)
            outb <- c(sapply(reft[first.Inode >= reft], function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% as.vector(),
                      PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,first.Inode < reft,drop=F],q = q,nt = n,reft = reft[first.Inode < reft],cal = cal))
          } else if (sum(first.Inode < reft) != 0) {
            outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          } else {
            TDqb = Diversity_profile.inc(c(ncol(x), Boots$boot_data[1:(nrow(x)+f0),B]), q)
            outb <- sapply(reft, function(i) TDqb*ifelse(cal == 'meanPD', 1, i)) %>% as.vector()}
          
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est, s.e. = ses,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reftime = rep(reft,each = length(q)))
      est
    })
  }
  Estoutput <- do.call(rbind,Estoutput) %>%
    mutate(Assemblage = rep(names(datalist),each = length(q)*tau_l),Method = 'Asymptotic',
           Type=cal) %>%
    select(Order.q,qPD,s.e.,qPD.LCL,qPD.UCL,Assemblage, 
           Method,Reftime,Type) %>%
    arrange(Reftime)
  Estoutput$qPD.LCL[Estoutput$qPD.LCL<0] = 0
  return(Estoutput)
}


PhD.q.est = function(ai,Lis, q, nt, reft, cal){
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  S <- length(ai)
  if(1 %in% q){
    ai_h1_I <- ai<=(nt-1)
    h1_pt2 <- rep(0,S)
    ai_h1 <- ai[ai_h1_I]
    h1_pt2[ai_h1_I] <- tibble(ai = ai) %>% .[ai_h1_I,] %>% mutate(diga = digamma(nt)-digamma(ai)) %>%
      apply(., 1, prod)/nt
  }
  if(2 %in% q){
    q2_pt2 <- unlist(ai*(ai-1)/nt/(nt-1))
  }
  if(sum(abs(q-round(q))!=0)>0 | max(q)>2) {
    deltas_pt2 <- sapply(0:(nt-1), function(k){
      ai_delt_I <- ai<=(nt-k)
      deltas_pt2 <- rep(0,S)
      deltas_pt2[ai_delt_I] <- delta_part2(ai = ai[ai_delt_I],k = k,n = nt)
      deltas_pt2
    }) %>% t() # n x S matrix of delta (2nd part)
  }
  Sub <- function(q,f1,f2,A,g1,g2,PD_obs,t_bar,Li){
    if(q==0){
      ans <- PD_obs+PDq0(nt,f1,f2,g1,g2)
    }else if(q==1){
      h2 <- PDq1_2(nt,g1,A)
      h1 <- sum(Li*h1_pt2)
      h <- h1+h2
      ans <- t_bar*exp(h/t_bar)
    }else if(q==2){
      #ans <- PDq2(as.matrix(tmpaL),nt,t_bar)
      ans <- t_bar^2/sum(Li*q2_pt2)
    }else{
      # timea <- Sys.time()
      k <- 0:(nt-1)
      deltas <- as.numeric(deltas_pt2 %*% Li)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(g1==0|A==1,0,(g1*((1-A)^(1-nt))/nt)*(A^(q-1)-round(sum(choose(q-1,k)*(A-1)^k), 12)))
      ans <- ((a+b)/(t_bar^q))^(1/(1-q))
      # timeb <- Sys.time()
      # print(timeb-timea)
    }
    return(ans)
  }
  est <- sapply(1:ncol(Lis),function(i){
    Li = Lis[,i]
    I1 <- which(ai==1&Li>0);I2 <- which(ai==2&Li>0)
    f1 <- length(I1);f2 <- length(I2)
    A <- ifelse(f2 > 0, 2*f2/((nt-1)*f1+2*f2), ifelse(f1 > 0, 2/((nt-1)*(f1-1)+2), 1))
    t_bar <- t_bars[i]
    PD_obs <- sum(Li)
    g1 <- sum(Li[I1])
    g2 <- sum(Li[I2])
    est <- sapply(q, function(q_) Sub(q = q_,f1 = f1, f2 = f2, A = A,g1 = g1,g2 = g2,PD_obs = PD_obs,t_bar = t_bar,Li = Li))
  })
  if(cal=='PD'){
    est <- as.numeric(est)
  }else if (cal=='meanPD'){
    est <- as.numeric(sapply(1:length(reft), function(i){
      est[,i]/reft[i]
    }))
  }
  return(est)
}



Boots.one = function(phylo, aL, datatype, nboot,reft, BLs, splunits = NULL){
  if(datatype=='abundance'){
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    names(data) <- rownames(BLs)[1:length(data)]
    n <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/n)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- (1-c) / sum((data/n)*(1- (data/n) )^n)
    
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (1-c) / f0 , f0 )
    if(length(p_hat0)>0) names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g0_hat <- sapply(1:length(reft), function(i){
      Li = BLs[,i]
      I1 <- which(aL$branch.abun==1&Li>0);I2 <- which(aL$branch.abun==2&Li>0)
      f1 <- length(I1);f2 <- length(I2)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      g0_hat <- ifelse(f1==0,0, 
                       ifelse(g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) ))
      if (is.na(g0_hat)) g0_hat = 0
      g0_hat
    })
    # g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    # g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    # g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    ###Notice that the species order of pL_b doesn't change even that of data changes. (property of phyBranchAL_Abu)
    pL_b <- phyBranchAL_Abu(phylo, p_hat, datatype,refT = reft[1])
    pL_b$treeNabu$branch.length <- pL_b$BLbyT[,1]
    pL_b <- pL_b$treeNabu
    pL_b[length(p_hat)+1,"branch.abun"] <- 1
    Li <- BLs
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li),byrow = T), Li[-(1:length(data)),,drop=F])
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.abun"]))
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }else if(datatype=='incidence_raw'){
    n <- splunits
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    names(data) <- rownames(BLs)[1:length(data)]
    u <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/u)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/u)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- u/n*(1-c) / sum((data/n)*(1- (data/n) )^n)
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (u/n) * (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    if(length(p_hat0)>0) names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g0_hat <- sapply(1:length(reft), function(i){
      Li = BLs[,i]
      I1 <- which(aL$branch.abun==1&Li>0);I2 <- which(aL$branch.abun==2&Li>0)
      f1 <- length(I1);f2 <- length(I2)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      g0_hat <- ifelse(f1==0,0, 
                       ifelse(g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) ))
      if (is.na(g0_hat)) g0_hat = 0
      g0_hat
    })
    # g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    # g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    # g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    pL_b <- phy_BranchAL_IncBootP(phylo = phylo, pdata = p_hat, refT = reft[1])
    pL_b <- pL_b$treeNincBP
    pL_b[length(p_hat)+1,"branch.incBP"] <- 1
    #pL_b$treeNincBP$branch.length <- pL_b$BLbyT[,1] # delete for now since length is a list instead of a matrix.
    #data_iB <- unlist(aL$branch.abun)
    #pL_b <- (data_iB/n) * (1-lambda*(1- (data_iB/n) )^n)
    Li <- BLs
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li),byrow = T),
                Li[-(1:length(data)),,drop=F])
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.incBP"]))
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }
  return(list(boot_data=boot_data,Li = Li, f0 = f0))
}




datainf <- function(data, datatype, phylotr, reft){
  if(datatype == "abundance"){
    new <- phyBranchAL_Abu(phylotr,data,datatype,reft)
    #new$treeNabu$branch.length <- new$BLbyT[,1]
    data <- data[data>0]
    ai <- new$treeNabu$branch.abun
    Lis <- new$BLbyT
    a1 <- sapply(1:ncol(Lis),function(i){
      Li = Lis[,i]
      I1 <- which(ai==1&Li>0);I2 <- which(ai==2&Li>0)
      f1 <- length(I1);f2 <- length(I2)
      PD_obs <- sum(Li)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      c(f1,f2,PD_obs,g1,g2)
    }) %>% matrix(nrow = 5) %>% t()
    a1 <- tibble('n' = sum(data),'S.obs' = length(data),'PD.obs' = a1[,3],
                 'f1*' = a1[,1],'f2*' = a1[,2], 'g1' = a1[,4],'g2' = a1[,5])
  }else if(datatype == 'incidence_raw'){
    new <- phyBranchAL_Inc(phylotr,data,datatype,reft)
    #new$treeNabu$branch.length <- new$BLbyT[,1]
    data <- data[rowSums(data)>0,colSums(data)>0,drop=F]
    ai <- new$treeNabu$branch.abun
    Lis <- new$BLbyT
    a1 <- sapply(1:ncol(Lis),function(i){
      Li = Lis[,i]
      I1 <- which(ai==1);I2 <- which(ai==2)
      f1 <- length(I1);f2 <- length(I2)
      PD_obs <- sum(Li)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      c(f1,f2,PD_obs, g1, g2)
    }) %>% matrix(nrow = 5) %>% t()
    a1 <- tibble('T' = ncol(data),'S.obs' = nrow(data),'PD.obs' = a1[,3],
                 'Q1*' = a1[,1],'Q2*' = a1[,2], 'R1' = a1[,4],'R2' = a1[,5])
  }
  return(a1)
  
}


