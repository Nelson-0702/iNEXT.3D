# PhDObs ------------------------------------------------------------------
datainf <- function(data, datatype, phylotr,reft){
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
    a1 <- tibble('nT' = ncol(data),'S.obs' = nrow(data),'PD.obs' = a1[,3],
                 'Q1*' = a1[,1],'Q2*' = a1[,2], 'R1' = a1[,4],'R2' = a1[,5])
  }
  return(a1)
  
}
color_nogreen <- function(n) {
  all <- c("red", "blue", "orange", "purple", "pink", "cyan", "brown", "yellow")
  all[1:n]
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
EmpPD <- function(datalist,datatype, phylotr, q, reft, cal, nboot, conf){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,reft = reft,cal = cal,nt = n) %>%
        c()
      
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
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal, nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q,reft = reft,cal = cal,nt = n) %>%
        c()
      
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
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,reft = reft,cal = cal,nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  Output <- tibble(Order.q = rep(rep(q, each=length(reft)),length(datalist)),
                   qPD = out[,1],qPD.LCL = out[,2], qPD.UCL = out[,3],
                   Assemblage = rep(nms, each=length(reft)*length(q)),
                   Method='Empirical',
                   Reftime = rep(reft,length(q)*length(datalist)),
                   Type=cal) %>%
    arrange(Reftime)
  return(Output)
}
EmpPD2 <- function(datalist,datatype, phylotr, q, reft, cal, nboot, conf){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  reft_complete <- reft
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[.>0]
      #divide reference time into two parts.
      first_parent <- aL$treeNabu %>% filter(tgroup=='Tip') %>% select(branch.length) %>% min
      reft_taxo <- reft[reft <= first_parent]
      reft <- reft[reft > first_parent]
      if(length(reft_taxo)>0) {
        aL$BLbyT <- aL$BLbyT[,-(1:length(reft_taxo))]
        ans_taxo <- TD.Tprofile(x=x,q=q,datatype=datatype,nboot=nboot,conf=conf,cal=cal,reft_taxo=reft_taxo)
        ans_taxo <- tibble(qPD = ans_taxo[,1],qPD.LCL = ans_taxo[,2], qPD.UCL = ans_taxo[,3],
                           Order.q = rep(q, each=length(reft_taxo)),Reftime = rep(reft_taxo,length(q)))
      }
      
      n <- sum(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, q=q,cal = cal,nt = n) %>%
        c()
      
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
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,cal = cal, nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output <- tibble(qPD = output[,1],qPD.LCL = output[,2], qPD.UCL = output[,3],
                       Order.q = rep(q, each=length(reft)),Reftime = rep(reft,length(q)))
      if(length(reft_taxo)>0) {
        output <- rbind(ans_taxo,output) %>% 
          arrange(Order.q,Reftime) %>% select(qPD,qPD.LCL,qPD.UCL) %>% as.matrix()
      }
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      #divide reference time into two parts.
      first_parent <- aL$treeNabu %>% filter(tgroup=='Tip') %>% select(branch.length) %>% min
      reft_taxo <- reft[reft <= first_parent]
      reft <- reft[reft > first_parent]
      if(length(reft_taxo)>0) {
        aL$BLbyT <- aL$BLbyT[,-(1:length(reft_taxo))]
        ans_taxo <- TD.Tprofile(x=x,q=q,datatype=datatype,nboot=nboot,conf=conf,cal=cal,reft_taxo=reft_taxo)
        ans_taxo <- tibble(qPD = ans_taxo[,1],qPD.LCL = ans_taxo[,2], qPD.UCL = ans_taxo[,3],
                           Order.q = rep(q, each=length(reft_taxo)),Reftime = rep(reft_taxo,length(q)))
      }
      
      n <- ncol(x)
      emp <- PD.Tprofile(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q,cal = cal,nt = n) %>%
        c()
      
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
          out_b <- PD.Tprofile(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q=q,cal = cal,nt = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output <- tibble(qPD = output[,1],qPD.LCL = output[,2], qPD.UCL = output[,3],
                       Order.q = rep(q, each=length(reft)),Reftime = rep(reft,length(q)))
      if(length(reft_taxo)>0) {
        output <- rbind(ans_taxo,output) %>% 
          arrange(Order.q,Reftime) %>% select(qPD,qPD.LCL,qPD.UCL) %>% as.matrix()
      }
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }
  Output <- tibble(Assemblage = rep(nms, each=length(reft_complete)*length(q)),
                   Order.q = rep(rep(q, each=length(reft_complete)),length(datalist)),
                   qPD = out[,1],qPD.LCL = out[,2], qPD.UCL = out[,3],
                   Reftime = rep(reft_complete,length(q)*length(datalist)),
                   Method='Empirical',
                   Type=cal) %>%
    arrange(Reftime)
  return(Output)
}


# Plott ------------------------------------------------------------------
Plott <- function(out){
  fort <- out
  fort$Order.q <- factor(paste0("q = ", fort$Order.q),levels = paste0("q = ", unique(fort$Order.q)))
  Assemblage <- unique(fort$Assemblage)
  ylab_ <- paste0(unique(fort$Method)," ",unique(fort$Type))
  if(length(Assemblage)==1){
    p2 <- ggplot(fort, aes(x=Reftime, y=qPD)) + theme_bw() +geom_line(size=1.5,aes(color=Order.q))+
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Order.q),linetype = 0,alpha=0.2)+
      xlab("Reference time")+ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
  }else{
    p2 <- ggplot(fort, aes(x=Reftime, y=qPD, color=Assemblage, linetype=Assemblage)) + theme_bw() + geom_line(size=1.5)  +
      geom_ribbon(aes(ymin=qPD.LCL,ymax=qPD.UCL,fill=Assemblage,alpha=0.2),linetype = 0,alpha=0.2)+
      scale_color_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
      scale_fill_manual(values = color_nogreen(length(unique(fort$Assemblage))))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      facet_wrap(~Order.q, scales = "free")+xlab("Reference time")+ylab(ylab_)
  }
  return(p2)
}


# asymPD ------------------------------------------------------------------
asymPD <- function(datalist, datatype, phylotr, q,reft, cal,nboot, conf){#change final list name
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  tau_l <- length(reft)
  if(datatype=="abundance"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[x>0]
      n <- sum(x)
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      #aL$treeNabu$branch.length <- aL$BLbyT[,1]
      #aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
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
          outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reftime = rep(reft,each = length(q)))
      est
    })
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[rowSums(x)>0,colSums(x)>0]
      n <- ncol(x)
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x,datatype,refT = reft)
      # aL$treeNabu$branch.length <- aL$BLbyT[,1]
      # aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      est <- PhD.q.est(ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,q = q,nt = n,reft = reft,cal = cal)
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
          outb <- PhD.q.est(ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],q = q,nt = n,reft = reft,cal = cal)
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,length(est))
      }
      est <- tibble(Order.q = rep(q,tau_l), qPD = est,
                    qPD.LCL = est - qtile*ses, qPD.UCL = est + qtile*ses,
                    Reftime = rep(reft,each = length(q)))
      est
    })
  }
  Estoutput <- do.call(rbind,Estoutput) %>%
    mutate(Assemblage = rep(names(datalist),each = length(q)*tau_l),Method = 'Asymptotic',
           Type=cal) %>%
    select(Order.q,qPD,qPD.LCL,qPD.UCL,Assemblage, 
           Method,Reftime,Type) %>%
    arrange(Reftime)
  Estoutput$qPD.LCL[Estoutput$qPD.LCL<0] = 0
  return(Estoutput)
}


# PhD.q.est ------------------------------------------------------------------
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


# Boots.one ------------------------------------------------------------------
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
    g0_hat <- sapply(1:length(reft), function(i){
      Li = BLs[,i]
      I1 <- which(aL$branch.abun==1&Li>0);I2 <- which(aL$branch.abun==2&Li>0)
      f1 <- length(I1);f2 <- length(I2)
      g1 <- sum(Li[I1])
      g2 <- sum(Li[I2])
      g0_hat <- ifelse(f1==0,0, 
                       ifelse(g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) ))
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


# inextPD ------------------------------------------------------------------
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
      covm = Coverage(x, datatype, m[[i]],n)
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- unique(invChatPD_abu(x = x,ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                        q = q,Cs = goalSC,n = n,reft = reft,cal = cal))
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
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            qPD_unc_b <- unique(invChatPD_abu(x = ai_B[isn0&tgroup_B=="Tip"],
                                              ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                              q = q,Cs = goalSC,n = n,reft = reft,cal = cal))$qPD
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
            covm_b <- Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
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
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ <- rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], m=m_,Method=method,Order.q=orderq,
                      qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                      SC=SC_,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                      Reftime = reft_,
                      Type=cal) %>%
        arrange(Reftime,Order.q,m)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    Type=cal,
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','goalSC','SC','m', 'Method', 'Order.q', 'qPD', 'qPD.LCL','qPD.UCL','Reftime',
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
      covm = Coverage(x, datatype, m[[i]], n)
      #====unconditional====
      if(unconditional_var){
        goalSC <- unique(covm)
        qPD_unc <- unique(invChatPD_inc(x = rowSums(x),ai = aL$treeNabu$branch.abun,Lis = aL$BLbyT,
                                        q = q,Cs = goalSC,n = n,reft = reft,cal = cal))
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
            covm_b = Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
            qPD_unc_b <- unique(invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],
                                              ai = ai_B[isn0],Lis = Li_b[isn0,,drop=F],
                                              q = q,Cs = goalSC,n = n,reft = reft,cal = cal))$qPD
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
            covm_b = Coverage(ai_B[isn0&tgroup_B=="Tip"], datatype, m[[i]],n)
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
      SC.LCL_ <- rep(covm-qtile*ses_cov,each = length(q)*length(reft))
      SC.UCL_ <- rep(covm+qtile*ses_cov,each = length(q)*length(reft))
      reft_ = rep(rep(reft,each = length(q)),length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], nt=m_,Method=method,Order.q=orderq,
                      qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
                      SC=SC_,SC.LCL=SC.LCL_,SC.UCL=SC.UCL_,
                      Reftime = reft_,
                      Type=cal) %>%
        arrange(Reftime,Order.q,nt)
      out_m$qPD.LCL[out_m$qPD.LCL<0] <- 0;out_m$SC.LCL[out_m$SC.LCL<0] <- 0
      out_m$SC.UCL[out_m$SC.UCL>1] <- 1
      if(unconditional_var){
        ses_pd_unc <- ses[-(1:(length(qPDm)+length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD-qtile*ses_pd_unc,qPD.UCL = qPD+qtile*ses_pd_unc,
                                    Type=cal,
                                    Assemblage = nms[i])
        id_C <- match(c('Assemblage','goalSC','SC','nt', 'Method', 'Order.q', 'qPD', 'qPD.LCL','qPD.UCL','Reftime',
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
    ans <- list(size_based = do.call(rbind,lapply(Estoutput, function(x) x$size_based)),
                coverage_based = do.call(rbind,lapply(Estoutput, function(x) x$coverage_based)))
  }else{
    ans <- list(size_based = do.call(rbind,lapply(Estoutput, function(x) x$size_based)))
  }
  return(ans)
}


# PhD.m.est ------------------------------------------------------------------
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
    mRPD = rbind(int.m, sapply(int.m, function(k) RPD(ai = ai,Lis = Lis,n = nt,m = k,q = q)))
  }
  
  if (cal == 'PD'){
    out <- sapply(m, function(mm){
      if(mm<nt){
        if(mm == round(mm)) { ans <- mRPD[-1,mRPD[1,] == mm] 
        } else { ans <- (ceiling(mm) - mm)*mRPD[-1, mRPD[1,] == floor(mm)] + (mm - floor(mm))*mRPD[-1, mRPD[1,] == ceiling(mm)] }
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
        if(mm == round(mm)) { ans <- mRPD[-1,mRPD[1,] == mm] 
        } else { ans <- (ceiling(mm) - mm)*mRPD[-1, mRPD[1,] == floor(mm)] + (mm - floor(mm))*mRPD[-1, mRPD[1,] == ceiling(mm)] }
      }else if(mm==nt){
        ans <- obs
      }else if(mm==Inf){
        ans <- asy
      }else{
        ans <- EPD(m = mm,obs = obs,asy = asy)
      }
      ans <- sapply(1:length(reft), function(i){
        ans[,i]/reft[i]
      })
      as.numeric(ans)
    })
  }
  out <- matrix(out,ncol = length(m))
  return(out)
}


# Coverage ------------------------------------------------------------------
Coverage = function(data, datatype, m, nt){
  if(datatype == "incidence_raw") datatype = "incidence"
  ifelse("matrix" %in% class(data) || "data.frame" %in% class(data), type <- "raw", type <- "numeric")
  ifelse(type == "raw", x <- rowSums(data), x <- data )
  if(type=="raw" || datatype=='incidence') u<-sum(x)
  x <- x[x>0]
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (nt - 1) / nt * f1 * (f1 - 1) / 2, (nt - 1) / nt * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, nt*f0.hat/(nt*f0.hat+f1), 1)
  Sub <- function(m){
    if(m < nt) {
      xx <- x[(nt-x)>=m]
      out <- 1-sum(xx / nt * exp(lgamma(nt-xx+1)-lgamma(nt-xx-m+1)-lgamma(nt)+lgamma(nt-m)))
    }
    if(m == nt) out <- 1-f1/nt*A
    if(m > nt) out <- 1-f1/nt*A^(m-nt+1)
    out
  }
  Sub2 <- function(m){
    if(m < nt) {
      xx <- x[(nt-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(nt-xx+1)-lgamma(nt-xx-m+1)-lgamma(nt)+lgamma(nt-m)))
    }
    if(m == nt) out <- 1-f1/u*A
    if(m > nt) out <- 1-f1/u*A^(m-nt+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype!='abundance', Sub2(i), Sub(i) )
  })
}


# RE_plot ------------------------------------------------------------------
RE_plot = function(data, type){
  #data <- as.data.frame(data)
  datatype <- ifelse(colnames(data$size_based[,2])=='m','abundance','incidence_raw')
  x <- ifelse(datatype=='incidence_raw', 'sampling units', "individuals")
  x_name <- colnames(data$size_based[,2])
  if(type == 1){
    data <- data$size_based
    id <- match(c(x_name,'Method','qPD','qPD.LCL','qPD.UCL','Assemblage','Order.q',
                  'Reftime'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- paste0("Number of ", x)
  }else if(type == 2){
    data <- data$size_based %>% filter(Order.q==1)
    id <- match(c(x_name,'Method','SC','SC.LCL','SC.UCL','Assemblage','Order.q',
                  'Reftime'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- paste0("Number of ", x)
    ylab_ <- "Sample Coverage"
  }else if(type == 3){
    data <- data$coverage_based
    id <- match(c('SC','Method','qPD','qPD.LCL','qPD.UCL','Assemblage','Order.q',
                  'Reftime'), names(data), nomatch = 0)
    output <- data[, id]
    xlab_ <- "Sample Coverage"
  }
  ylab_ <- unique(data$Type)
  title <- c("Sample-size-based sampling curve", "Sample completeness curve","Coverage-based sampling curve")
  title <- title[type]
  
  
  Assemblage <- unique(data$Assemblage)
  colnames(output) <- c("x", "Method", "y", "LCL", "UCL", "Assemblage","Order.q",'Reftime')
  output$Method <- as.character(output$Method)
  output$Assemblage <- as.character(output$Assemblage)
  output$Reftime <- round(output$Reftime,3)
  output$Reftime <- factor(paste0('Ref.time = ',output$Reftime),
                           levels = paste0('Ref.time = ',unique(output$Reftime)))
  output_obser <- output %>% filter(Method=="Observed")
  output$Method[output$Method=="Observed"] <- "Rarefaction"
  
  # odr_grp <- as_labeller(c(`0` = "q = 0", `1` = "q = 1",`2` = "q = 2"))
  # odr_grp <- list(
  #   '0'="q = 0",
  #   '1'="q = 1",
  #   '2'="q = 2"
  # )
  # refts <- unique(output$Reftime)
  # reft_grp <- sapply(refts, function(i){
  #   paste0('Ref.time = ',refts[i])
  # })
  # names(reft_grp) <- as.character(refts)
  # plot_labeller <- function(variable,value){
  #   if (variable=='Order.q') {
  #     return(odr_grp[value])
  #   } else {
  #     return(reft_grp[value])
  #   }
  # }
  mylab <- labeller(
    Order.q = c(`0` = "q = 0", `1` = "q = 1",`2` = 'q = 2')
  )
  if(length(Assemblage) == 1){
    outp <- ggplot(output, aes(x = x, y = y))+ theme_bw() +
      geom_ribbon(aes(ymin = LCL, ymax = UCL),fill="#F8766D",alpha=0.2)+geom_line(size=1.5, aes(x = x, y = y, linetype=Method),color="#F8766D")+
      geom_point(size=5, data=output_obser,color="#F8766D")+xlab(xlab_)+ylab(ylab_)+
      scale_linetype_manual(values = c("solid", "22"), name="Method",breaks=c("Rarefaction", "Extrapolation"), labels=c("Rarefaction", "Extrapolation"))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
      ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
    if(type!=3) outp <- outp + facet_wrap(Reftime~Order.q,scales = "free_y",labeller=mylab)
    else if (type==3) outp <- outp + facet_wrap(~Reftime,scales = "free_y")
    #if(length(unique(output$Reftime))==1) outp <- outp + theme(strip.background = element_blank(), strip.text.x = element_blank())
  }else{
    outp <- ggplot(output, aes(x = x, y = y, color=Assemblage)) + theme_bw() +
      geom_line(size=1.5, aes(x = x, y = y, color=Assemblage, linetype=Method))+
      scale_color_manual(values = color_nogreen(length(unique(output$Assemblage))))+
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Assemblage), alpha=0.2, colour=NA)+
      scale_fill_manual(values = color_nogreen(length(unique(output$Assemblage))))+
      geom_point(size=5, data=output_obser, aes(shape=Assemblage))+xlab(xlab_)+ylab(ylab_)+
      scale_linetype_manual(values = c("solid", "22"), name="Method",breaks=c("Rarefaction", "Extrapolation"), labels=c("Rarefaction", "Extrapolation"))+
      theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"),legend.box = "vertical")+
      ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
    if(type!=2)  outp <- outp + facet_wrap(Reftime~Order.q,scales = "free_y",labeller=mylab)
    else if (type == 2) outp <- outp + facet_wrap(~Reftime,scales = "free_y")
    #if(length(unique(output$Reftime))==1) outp <- outp + theme(strip.background = element_blank(), strip.text.x = element_blank())
  }
  return(outp)
}


# invChatPD ------------------------------------------------------------------
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
      est <- est %>% mutate(qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
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
          # atime <- Sys.time()
          ai_B <- Boots$boot_data[,B]
          isn0 <- ai_B>0
          # isn0 <- as.vector(aL_table_b[,1]>0)
          # Li_b_tmp <- Li_b[isn0,]
          # aL_table_b <- aL_table_b[isn0,]
          est_b <- invChatPD_inc(x = ai_B[isn0&tgroup_B=="Tip"],ai = ai_B[isn0],
                                 Lis = Li_b[isn0,,drop=F],q = q,Cs = level,
                                 n = n, reft = reft, cal = cal)$qPD
          return(est_b)
        }) %>% matrix(.,nrow = length(q)*length(reft)*length(level)) %>% apply(., 1, sd)
      }else{
        ses <- rep(NA,nrow(est))
      }
      est <- est %>% mutate(qPD.LCL=qPD-qtile*ses,qPD.UCL=qPD+qtile*ses)
    }) %>% do.call(rbind,.)
  }
  Assemblage = rep(names(datalist), each = length(q)*length(reft)*length(level))
  out <- out %>% mutate(Assemblage = Assemblage, Type=cal)
  if(datatype=='abundance'){
    out <- out %>% select(Assemblage,goalSC,SC,m,Method,Order.q,qPD,qPD.LCL,qPD.UCL,
                          Reftime,Type)
  }else if(datatype=='incidence_raw'){
    out <- out %>% select(Assemblage,goalSC,SC,nt,Method,Order.q,qPD,qPD.LCL,qPD.UCL,
                          Reftime,Type)
  }
  out$qPD.LCL[out$qPD.LCL<0] <- 0
  rownames(out) <- NULL
  out
}


# invChatPD_abu ------------------------------------------------------------------
invChatPD_abu <- function(x,ai,Lis, q, Cs, n, reft, cal){
  #x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
  refC = Coverage(x, 'abundance', n, n)
  f <- function(m, C) abs(Coverage(x, 'abundance', m, n) - C)
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
      }else if(f1 == 1 & f2 == 0) {
        A <- 1
      }else if(f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- ifelse(A==1,0,(log(n/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
    }
    mm
  })
  mm[mm < 1] <- 1
  mm[which(round(mm) - n <= 1)] = round(mm[which(round(mm) - n <= 1)]) 
  SC <- Coverage(x, 'abundance', mm, n)
  out <- as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,reft=reft,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  goalSC <- rep(Cs,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reftime = rep(rep(reft,each = length(q)),length(Cs))
  tibble(m = m,Method = method,Order.q = order,
         qPD = out,SC=SC,goalSC = goalSC,Reftime = Reftime)
}


# invChatPD_inc ------------------------------------------------------------------
invChatPD_inc <- function(x,ai,Lis, q, Cs, n, reft, cal){ # x is a matrix
  #x <- unlist(aL_table$branch.abun[aL_table$tgroup=="Tip"])
  refC = Coverage(x, 'incidence', n, n)
  f <- function(m, C) abs(Coverage(x, 'incidence', m, n) - C)
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
      }else if(f1 == 1 & f2 == 0) {
        A <- 1
      }else if(f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- ifelse(A==1,0,(log(U/f1) + log(1 - cvrg))/log(A) - 1)
      mm <- n + mm
    }
    mm
  })
  mm[mm < 1] <- 1
  mm[which(round(mm) - n <= 1)] = round(mm[which(round(mm) - n <= 1)]) 
  SC <- Coverage(x, 'incidence', mm, n)
  out <-  as.numeric(PhD.m.est(ai = ai,Lis = Lis,m = mm,q = q,nt = n,reft = reft,cal = cal))
  method <- ifelse(mm>n,'Extrapolation',ifelse(mm<n,'Rarefaction','Observed'))
  method <- rep(method,each = length(q)*ncol(Lis))
  m <- rep(mm,each = length(q)*ncol(Lis))
  order <- rep(q,ncol(Lis)*length(mm))
  SC <- rep(SC,each = length(q)*ncol(Lis))
  goalSC <- rep(Cs,each = length(q)*ncol(Lis))
  reft <- as.numeric(substr(colnames(Lis),start = 2,stop = nchar(colnames(Lis))))
  Reftime = rep(rep(reft,each = length(q)),length(Cs))
  tibble(nt = m,Method = method,Order.q = order,
         qPD = out,SC=SC,goalSC = goalSC, Reftime = Reftime)
}


# TD.Tprofile ------------------------------------------------------------------
#' @importFrom stats rmultinom
TD.Tprofile <- function(x,q,datatype="abundance",nboot=50,conf=0.95,cal,reft_taxo){
  qtile <- qnorm(1-(1-conf)/2)
  if(cal=='meanPD') reft_taxo <- rep(1,length(reft_taxo))
  reft_taxo_dummy <- rep(reft_taxo,length(q))
  if(datatype=="abundance"){
    dq <- rep(Diversity_profile_MLE(x,q),each = length(reft_taxo))*reft_taxo_dummy
    if(nboot>1){
      Prob.hat <- EstiBootComm.Ind(x)
      ses <- sapply(1:length(reft_taxo), function(l){
        Abun.Mat <- rmultinom(nboot, sum(x), Prob.hat)
        se <- apply(matrix(apply(Abun.Mat, 2, function(xb) Diversity_profile_MLE(xb,q))*reft_taxo[l],
                           nrow = length(q)), 1, sd, na.rm=TRUE)
        c(se)
      }) %>% matrix(.,nrow = length(q),ncol = length(reft_taxo)) %>% t %>% c
    }else{ses = rep(NA,length(reft_taxo)*length(q))}
  }else if(datatype=='incidence_raw'){
    nT <- ncol(x)
    x <- c(nT,rowSums(x))
    dq <- rep(Diversity_profile_MLE.inc(x,q),each = length(reft_taxo))*reft_taxo_dummy
    if(nboot>1){
      Prob.hat <- EstiBootComm.Sam(x)
      ses <- sapply(1:length(reft_taxo), function(l){
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          se = NA
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          se <- apply(matrix(apply(Abun.Mat, 2, function(yb) Diversity_profile_MLE.inc(yb,q))*reft_taxo[l],
                         nrow = length(q)), 1, sd, na.rm=TRUE)
        }
        se
      }) %>% matrix(.,nrow = length(q),ncol = length(reft_taxo)) %>% t %>% c
    
    }else{ses = rep(NA,length(reft_taxo)*length(q))}
  }
  output <- cbind(dq,dq-qtile*ses,dq+qtile*ses)
  output
}


# TD.q.est ------------------------------------------------------------------
TD.q.est <- function(x,q,datatype="abundance",nboot=50,conf=0.95,cal,reft_taxo){
  qtile <- qnorm(1-(1-conf)/2)
  if(cal=='meanPD') reft_taxo <- rep(1,length(reft_taxo))
  reft_taxo_dummy <- rep(reft_taxo,length(q))
  if(datatype=="abundance"){
    dq <- rep(Diversity_profile(x,q),each = length(reft_taxo))*reft_taxo_dummy
    if(nboot>1){
      Prob.hat <- EstiBootComm.Ind(x)
      ses <- sapply(1:length(reft_taxo), function(l){
        Abun.Mat <- rmultinom(nboot, sum(x), Prob.hat)
        se <- apply(matrix(apply(Abun.Mat, 2, function(xb) Diversity_profile(xb,q))*reft_taxo[l],
                           nrow = length(q)), 1, sd, na.rm=TRUE)
        c(se)
      }) %>% matrix(.,nrow = length(q),ncol = length(reft_taxo)) %>% t %>% c
    }else{ses = rep(NA,length(reft_taxo)*length(q))}
  }else if(datatype=='incidence_raw'){
    nT <- ncol(x)
    x <- c(nT,rowSums(x))
    dq <- rep(Diversity_profile.inc(x,q),each = length(reft_taxo))*reft_taxo_dummy
    if(nboot>1){
      Prob.hat <- EstiBootComm.Sam(x)
      ses <- sapply(1:length(reft_taxo), function(l){
        Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, nT, p)))
        Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)),ncol=nboot)
        tmp <- which(colSums(Abun.Mat)==nT)
        if(length(tmp)>0) Abun.Mat <- Abun.Mat[,-tmp]
        if(ncol(Abun.Mat)==0){
          se = NA
          warning("Insufficient data to compute bootstrap s.e.")
        }else{		
          se <- apply(matrix(apply(Abun.Mat, 2, function(yb) Diversity_profile.inc(yb,q))*reft_taxo[l],
                             nrow = length(q)), 1, sd, na.rm=TRUE)
        }
        se
      }) %>% matrix(.,nrow = length(q),ncol = length(reft_taxo)) %>% t %>% c
      
    }else{ses = rep(NA,length(reft_taxo)*length(q))}
  }
  output <- cbind(dq,dq-qtile*ses,dq+qtile*ses)
  output
}

