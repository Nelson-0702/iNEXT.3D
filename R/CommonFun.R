# Transform incidence raw data to incidence frequencies (iNEXT input format) 
# 
# \code{as.incfreq}: transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units).
# @param x a \code{data.frame} or \code{matirx} of species by sites presence-absence matrix.
# @return a \code{vector} of species incidence frequencies, the first entry of the input data must be total number of sampling units.
# @examples
# data(ciliates)
# as.incfreq(ciliates)
as.incfreq <- function(data, nT = NULL) {
  
  if (inherits(data, c("data.frame", "matrix"))) {
    if (sum(data > 1) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
    
    if(is.null(nT)) nT = ncol(data)
    if(inherits(nT, 'data.frame')) nT = unlist(nT)
    mydata = list()
    if(ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(nT)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp+nT[i]
    }
    if(is.null(names(nT))) {
      names(mydata) <- paste0("Assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
    data = lapply(mydata, function(i){
      out = c('nT' = ncol(i), rowSums(i))
      return(out)
    })
  } else if (inherits(data, "list")) 
    data <- lapply(data, function(i) {
      if (sum(i > 1) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
      c('nT' = ncol(i), rowSums(i))}
      )
  if (length(data) == 1) data = data[[1]]
  return(data)
}


# Abundance-based or Incidence-based sample coverage
# 
# \code{Coverage} Estimation of abundance-based or incidence-based sample coverage function
# 
# @param x a vector of species abundances, a vector of species incidence-based frequency, or a matrix/data.frame of species incidence-raw data, the first entry is the total number of sampling units, followed by the speceis incidences abundances.
# @param datatype a selection for 'abundance', 'incidence_freq' and 'incidence_raw'.
# @param m a integer vector of rarefaction/extrapolation sample size or sample units
# @return a vector of estimated sample coverage function
# @export
Coverage = function(data, datatype, m){
  if (!(datatype %in% c('abundance', 'incidence_freq', 'incidence_raw')))
    stop("Invalid Coverage datatype", call. = FALSE)
  
  if (datatype == 'incidence_raw') {data = as.incfreq(data); datatype = 'incidence_freq'}
  n <- ifelse(datatype=='incidence_freq', data[1], sum(data) )
  if(datatype == "incidence_freq"){
    x <- data[-1]
    u<-sum(x)
  }else if(datatype == "abundance"){
    x <- data
  }
  x <- x[x>0]
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    if(m < n) {
      if (m == round(m)) {
        xx <- x[(n-x)>=m]
        out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
      } else {
        cbym = rbind(c(floor(m), ceiling(m)), 
                     sapply(c(floor(m), ceiling(m)), function(k) {
                       xx <- x[(n-x)>=k]
                       if (k == n) 1-f1/n*A else 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-k+1)-lgamma(n)+lgamma(n-k)))
                     }))
        out <- (ceiling(m) - m)*cbym[-1, cbym[1,] == floor(m)] + (m - floor(m))*cbym[-1, cbym[1,] == ceiling(m)] 
      }
      }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  Sub2 <- function(m){
    if(m < n) {
      if (m == round(m)) {
        xx <- x[(n-x)>=m]
        out <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
      } else {
        cbym = rbind(c(floor(m), ceiling(m)), 
                     sapply(c(floor(m), ceiling(m)), function(k) {
                       xx <- x[(n-x)>=k]
                       if (k == n) 1-f1/u*A else 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-k+1)-lgamma(n)+lgamma(n-k)))
                     }))
        out <- (ceiling(m) - m)*cbym[-1, cbym[1,] == floor(m)] + (m - floor(m))*cbym[-1, cbym[1,] == ceiling(m)] 
      }
    }
    if(m == n) out <- 1-f1/u*A
    if(m > n) out <- 1-f1/u*A^(m-n+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype!='abundance', Sub2(i), Sub(i) )
  })
}



#' Printing iNEXT3D object
#' 
#' \code{print.iNEXT3D}: Print method for objects inheriting from class "iNEXT3D"
#' @param x an \code{iNEXT3D} object computed by \code{iNEXT3D}.
#' @param ... additional arguments.
#' @return a list of three objects (see \code{iNEXT3D} for more details) with simplified outputs and notes.
#' @export
print.iNEXT3D <- function(x, ...){
  site.n <- nrow(x[[1]])
  order.n <- paste(unique(x[[2]]$size_based$Order.q), collapse = ", ")
  cat("Compare ", site.n, " assemblages with Hill number order q = ", order.n,".\n", sep = "")
  cat("$class: iNEXT3D\n\n")
  cat(names(x)[1], ": basic data information\n", sep = "")
  print(x[[1]])
  cat("\n")
  cat(names(x)[2],": diversity estimates with rarefied and extrapolated samples.\n", sep = "")
  cat("$size_based (LCL and UCL are obtained for fixed size.)\n")
  cat("\n")
  
  res <- lapply(x[[2]], function(y){
    Assemblages <- unique(x[[2]]$size_based$Assemblage)
    tmp <- lapply(1:length(Assemblages),function(i){
      # y_each <- subset(y, Assemblage==Assemblages[i])
      y_each <- y[y$Assemblage==Assemblages[i],]
      m <- quantile(unlist(y_each[,3]), type = 1)
      y_each[unlist(y_each[,3]) %in% m,]
    })
    do.call(rbind,tmp)
  })
  
  print(data.frame(res[[1]]))
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage in each order q; call iNEXT3D.object$", names(x)[2],
      "$size_based to view complete output.\n", sep = "")
  cat("\n")
  cat("$coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)\n")
  cat("\n")
  print(data.frame(res[[2]]))
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage in each order q; call iNEXT3D.object$", names(x[2]), 
      "$coverage_based to view complete output.\n", sep = "")
  cat("\n")
  cat(names(x)[3], ": asymptotic diversity estimates along with related statistics.\n", sep = "")
  print(x[[3]])
  return(invisible())
}


# check datatype and transform incidence_raw to incidence_freq
# 
# \code{check.datatype}
# 
# @param data input data
# @param datatype data type
# @param nT the vector of sampling units for each assemblage
# @param to.datalist a binary choice whether transform data to datalist
# @param raw.to.inci a binary choice whether transform incidence raw data to incidence frequency data
# @return a list of datatype, matrix data, and nT
# @export

check.datatype <- function(data, datatype, nT = nT, to.datalist = FALSE, raw.to.inci = TRUE) {
  
  if(datatype == "incidence") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.") #datatype = "incidence_freq" or   
  # if(datatype == "incidence_freq") stop("iNEXT.3D can only accept 'datatype = incidence_raw'.") #datatype = "incidence_freq" or   
  
  DATATYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, DATATYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, DATATYPE)
  
  if(sum(nT <= 0) != 0) stop("Number of sampling units should be a positive value.", call. = FALSE)
  
  if (datatype == "incidence_raw" & raw.to.inci == TRUE) {
    
    if (!inherits(data, "list")) 
      data = as.incfreq(data, nT = nT) else if (length(data) != 1)
        data = as.incfreq(data, nT = nT) else {
        tmp = names(data)
        data = list(as.incfreq(data, nT = nT))
        names(data) = tmp
      }
    datatype = "incidence_freq"
      }
  
  if (datatype == "incidence_raw") {
    
    if (inherits(data, c("data.frame", "matrix"))) {
      data = as.matrix(data)
      if (is.null(nT)) nT = c('Assemblage_1' = ncol(data))
      if (is.null(names(nT))) names(nT) = paste0("Assemblage_", 1:length(nT))
    }
    
    if (inherits(data, c("numeric", "integer", "double"))) {
      data = as.matrix(data)
      nT = c('Assemblage_1' = 1)
    }
    
    if (inherits(data, "list")) {
      data = lapply(data, function(i) data.frame(i))
      data2 = lapply(data, function(i) {
        i$species = rownames(i)
        return(i) 
      })
      nT = as.vector(sapply(data, ncol))
      names(nT) = if (is.null(data)) paste0("Assemblage_", 1:length(data)) else names(data)
      
      data = data2[[1]]
      if (length(data2) > 1) {
        for(i in 2:length(data2)){
          data = full_join(data, data2[[i]], by = "species")
        }
      }
      data[is.na(data)] = 0
      rownames(data) = data$species
      data = data %>% select(-species)
      
    }
    
    if (ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
    if (sum(data > 1) != 0) stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
  }
  
  if (datatype != "incidence_raw") {
    
    if(inherits(data, "list")){
      
      if(length(data) == 1){
        
        dat = as.matrix(data[[1]])
        if (is.null(names(data))) colnames(dat) = "Assemblage_1" else colnames(dat) = names(data)
        data = dat
        
      } else {
        region_names = if (is.null(names(data))) paste0("Assemblage_", 1:length(data)) else names(data)
        
        data2 = lapply(data, function(x) {
          if (is.null(names(x)) & datatype == 'abundance') names(x) = paste('Species', 1:length(x), sep = '')
          if (is.null(names(x)) & datatype == 'incidence_freq') names(x) = c('nT', paste('Species', 1:(length(x)-1), sep = ''))
          
          x = as.matrix(x)
          x = data.frame('species' = rownames(x), x)
          
          return(x)
        })
        datam = data2[[1]]
        for(i in 2:length(data2)){
          datam = data.frame(full_join(datam, data2[[i]], by = "species"))
        }
        datam[is.na(datam)] = 0
        datam = column_to_rownames(datam, var = "species")
        names(datam) = region_names
        
        if (is.null(names(data[[1]]))) rownames(datam) = NULL
        data = datam
      }
      
    } else if (inherits(data, c("numeric", "integer", "double"))) {
      data = as.matrix(data)
      colnames(data) = 'Assemblage_1'
    }
    
    data = as.matrix(data)
    
    if (datatype == "incidence_freq") nT = data[1,]
    
    if (to.datalist == TRUE) {
      datalist <- lapply(1:ncol(data), function(i)  x <- data[,i])
      names(datalist) = colnames(data)
      data = datalist
    }
    
  }
  
  if(inherits(nT, 'data.frame')) nT = unlist(nT)
  if(datatype != "abundance" & sum(nT <= 3) != 0) stop("Number of sampling units of some assemblages is too less. Please add more sampling units data.", call. = FALSE)
  
  return(list(datatype, data, nT))
}


# check dist and transform data matrix to data list
# 
# \code{check.dist}
# 
# @param data input a data matrix
# @param datatype data type
# @param distM a symmetric distance matrix
# @param threshold a value between zero and one
# @return a list of threshold, distance matrix, and datalist
# @export

check.dist <- function(data, datatype, distM, threshold) {
  distM = as.matrix(distM)
  
  if (datatype == 'incidence_freq'){
    nT <- data[1,]
    data <- data[-1,,drop =FALSE]
  }
  
  if (is.null(rownames(data)) | is.null(rownames(distM))){
    warning('The species names are not provided in data or distance matrix.', call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <- paste0('Species', 1:nrow(data))
  } else {
    if (sum(rownames(data) %in% rownames(distM)) != nrow(data))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  
  distM = distM[rownames(distM) %in% rownames(data), colnames(distM) %in% rownames(data)]
  
  if (nrow(data) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix", call. = FALSE)
  
  order_sp <- match(rownames(data),rownames(distM))
  
  distM <- distM[order_sp, order_sp]
  distM <- distM[rowSums(data)>0, rowSums(data)>0]
  data <- data[rowSums(data)>0, , drop=FALSE]
  
  if(datatype == 'incidence_freq'){
    data <- rbind(nT, data)
  }
  
  if (is.null(threshold)) {
    
    if (datatype == 'abundance') {
      
      tmp = rowSums(data)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
      
    } else if (datatype == 'incidence_freq') {
      
      tmp = rowSums(data)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
      
    }
    threshold <- sum ( (tmp %*% t(tmp) ) * distM)
    
  } else if(sum(threshold<0) > 0 | sum(threshold>1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean.", call. = FALSE)
  }
  
  datalist <- lapply(1:ncol(data), function(i)  x <- data[,i])
  names(datalist) = colnames(data)
  
  
  return(list(threshold, distM, datalist))
}


# check tree and transform data matrix to data list
# 
# \code{check.tree}
# 
# @param data input a data matrix
# @param datatype data type
# @param tree a phylo tree
# @param reftime a positive value
# @return a list of reftime, tree, and datalist
# @export

check.tree <- function(data, datatype, tree, reftime, nT) {
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  
  # if(sum(duplicated(tree$tip.label))>0)
  #   stop("The phylo tree should not contains duplicated tip labels, please remove them.", call. = FALSE)
  
  if( is.null(rownames(data)) )
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  pool.name <- rownames(data)
  if (sum(pool.name %in% tree$tip.label) != length(pool.name))
    stop("Data and tree tip label contain unmatched species", call. = FALSE)
  
  tip <- tree$tip.label[-match(pool.name, tree$tip.label)]
  mytree <- drop.tip(tree, tip)
  
  # if (is.null(reftime)) reftime <- ifelse(is.ultrametric(mytree), get.rooted.tree.height(mytree), max(ape::node.depth.edgelength(mytree))) else 
  #   reftime <- reftime
  if (is.null(reftime)) reftime <- max(ape::node.depth.edgelength(mytree)) else reftime <- reftime
  
  reftime <- sort(unique(reftime))
  
  if (sum(reftime<=0) > 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE) }
  
  if (datatype == "incidence_raw") {
    
    datalist = list()
    
    ntmp <- 0
    for(i in 1:length(nT)){
      datalist[[i]] <- data[,(ntmp+1):(ntmp+nT[i])]
      ntmp <- ntmp + nT[i]
    }
    
    names(datalist) = names(nT)
    
  } else if (datatype == 'abundance') {
    
    datalist <- lapply(1:ncol(data), function(i)  x <- data[,i])
    names(datalist) = colnames(data)
    
  }
  
  return(list(reftime, mytree, datalist))
}

check.q <- function(q) {
  
  if(!inherits(q, "numeric"))
    stop("invalid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q", call. = FALSE)
    q <- q[q >= 0]
  }
  
  return(q)
}

check.conf <- function(conf) {
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  return(conf)
}

check.nboot <- function(nboot) {
  
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  return(nboot)
}

check.base <- function(base) {
  
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) stop("invalid datatype")
  if (pmatch(base, BASE) == -1) stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  return(base)
}

check.PDtype <- function(type) {
  
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of phylogenetic diversity type, please use either PD or meanPD.", call. = FALSE)
  
  return(type)
}

check.size <- function(data, datatype, size, endpoint, knots) {
  
  if (length(knots) != length(data)) knots <- rep(knots,length(data))
  
  if (is.null(size)) {
    
    if (is.null(endpoint)) {
      
      if (datatype == "abundance") {
        endpoint <- sapply(data, function(x) 2*sum(x))
      } else if (datatype == "incidence_freq"){
        endpoint <- sapply(data, function(x) 2*x[1])
      } else if (datatype == "incidence_raw"){
        endpoint <- sapply(data, function(x) 2*ncol(x))
      }
      
    } else {
      
      if(length(endpoint) != length(data)){
        endpoint <- rep(endpoint, length(data))
      }
      
    }
    
    size <- lapply(1:length(data), function(i){
      
      if (datatype == "abundance") {
        ni <- sum(data[[i]])
      } else if (datatype == "incidence_freq") {
        ni <- data[[i]][1]
      } else if (datatype == "incidence_raw") {
        ni <- ncol(data[[i]])
      }
      
      if(endpoint[i] <= ni){
        mi <- floor(seq(1,endpoint[i], length.out = knots[i]))
      }else{
        mi <- floor(c(seq(1, ni, length.out = floor(knots[i]/2)), seq(ni+1, endpoint[i], length.out = knots[i]-floor(knots[i]/2))))
      }
      
      if(sum(mi < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
      unique(mi)
    })
    
  } else {
    
    if (inherits(size, c("numeric", "integer", "double"))) {
      size <- list(size = size)
    }
    
    if (length(size) != length(data)) size <- lapply(1:length(data), function(x) size[[1]])
    size <- lapply(1:length(data),function(i){
      if(datatype == "abundance") {
        ni <- sum(data[[i]])
      } else if (datatype == "incidence_freq") {
        ni <- data[[i]][1]
      } else if(datatype == "incidence_raw"){
        ni <- ncol(data[[i]])
      }
      
      if ( (sum(size[[i]] == ni) == 0) & (sum(size[[i]] > ni) != 0) & (sum(size[[i]] < ni) != 0) ) 
        mi <- sort(c(ni,size[[i]])) else mi <- sort(size[[i]])
      
      if(sum(mi < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
      unique(mi)
    })
  }
  
  
  return(size)
}

check.level <- function(data, datatype, base, level) {
  
  if (is.null(level) & base == 'size') {
    
    if (datatype == "abundance") {
      level <- sapply(data, function(x) 2*sum(x))
    } else if(datatype == "incidence_freq") {
      level <- sapply(data, function(x) 2*x[1])
    } else if(datatype == "incidence_raw") {
      level <- sapply(data, function(x) 2*ncol(x))
    }
    
    level <- min(level)
    
  } else if (is.null(level) & base == 'coverage') {
    
    if (datatype=='abundance') {
      
      level <- sapply(data,function(x) {
        ni <- sum(x)
        Coverage(data = x, datatype = datatype, m = 2*ni)
      })
      
    } else if (datatype=='incidence_freq') {
      
      level <- sapply(data,function(x){
        ni <- x[1]
        Coverage(data = x,datatype = datatype,m = 2*ni)
      })
      
    } else if (datatype=='incidence_raw') {
      
      level <- sapply(data,function(x) {
        ni <- ncol(x)
        Coverage(data = x, datatype = datatype, m = 2*ni)
      })
      
    }
    
    level <- min(level)
  }
  
  if(base == "size" & sum(level < 0) > 0) stop("Sample size (or number of sampling units) cannot be a negative value.", call. = FALSE)
  if(base == "coverage" & sum(level < 0 | level > 1) > 0) stop("The sample coverage values should be between zero and one.", call. = FALSE)  
  
  return(level)
}

check.FDcut_number <- function(FDcut_number) {
  
  if (!inherits(FDcut_number, "numeric"))
    stop("invalid class of FD cut number, FDcut_number should be a postive value.", call. = FALSE)
  if (FDcut_number < 2)
    stop("invalid FDcut_number, FDcut_number should be a postive value larger than one.", call. = FALSE)
  if (length(FDcut_number) > 1)
    stop("FDcut_number only accept a value instead of a vector.", call. = FALSE)
  
  return(FDcut_number)
}

