# iNEXT.3D -------------------------------------------------------------------
#' Printing iNEXT.3D object
#' 
#' \code{print.iNEXT3D}: Print method for objects inheriting from class "iNEXT3D"
#' @param x an \code{iNEXT3D} object computed by \code{\link{iNEXT3D}}.
#' @param ... additional arguments.
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
      m <- quantile(unlist(y_each[,2]), type = 1)
      y_each[unlist(y_each[,2]) %in% m,]
    })
    do.call(rbind,tmp)
  })
  
  print(res[[1]])
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$", names(x)[2],
      "$size_based to view complete output.\n", sep = "")
  cat("\n")
  cat("$coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)\n")
  cat("\n")
  print(res[[2]])
  cat("\n")
  cat("NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$", names(x[2]), 
      "$coverage_based to view complete output.\n", sep = "")
  cat("\n")
  cat(names(x)[3], ": asymptotic diversity estimates along with related statistics.\n", sep = "")
  print(x[[3]])
  return(invisible())
}


# check.datatype -------------------------------------------------------------------
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
  if(datatype == "incidence") stop('Please try datatype = "incidence_freq" or datatype = "incidence_raw".')  
  DATATYPE <- c("abundance", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, DATATYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, DATATYPE)
  
  if (datatype == "incidence_raw" & raw.to.inci == TRUE) {data = as.incfreq(data, nT = nT); datatype = "incidence_freq"}
  
  if (datatype == "incidence_raw") {
    
    if (class(data)[1] == "numeric"|class(data)[1] == "integer"|class(data)[1] == "double") {
      data = as.matrix(data)
      nT = c('Assemblage_1' = 1)
    }
    
    if (class(data)[1] == "list") {
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
      data = column_to_rownames(data, var = "species")
      
    }
    
    if (ncol(data) != sum(nT)) stop("Number of columns does not euqal to the sum of nT (number of sampling units for each assemblage).", call. = FALSE)
  }
  
  if (datatype != "incidence_raw") {
    
    if(class(data)[1] == "list"){
      
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
      
    } else if (class(data)[1] == "numeric"|class(data)[1] == "integer"|class(data)[1] == "double") {
      data = as.matrix(data)
      colnames(data) = 'Assemblage_1'
    }
    
    data = as.matrix(data)
    
    if (to.datalist == TRUE) {
      datalist <- lapply(1:ncol(data), function(i)  x <- data[,i])
      names(datalist) = colnames(data)
      data = datalist
    }
    
  }
  
  return(list(datatype, data, nT))
}


# check.dist -------------------------------------------------------------------
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


# check.tree -------------------------------------------------------------------
# check tree and transform data matrix to data list
# 
# \code{check.dist}
# 
# @param data input a data matrix
# @param datatype data type
# @param distM a symmetric distance matrix
# @param threshold a value between zero and one
# @return a list of reftime, tree, and datalist
# @export

check.tree <- function(data, datatype, tree, reftime, nT) {
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  
  if( is.null(rownames(data)) )
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
  
  data <- data[rowSums(data)>0,,drop=FALSE]
  
  pool.name <- rownames(data)
  tip <- tree$tip.label[-match(pool.name, tree$tip.label)]
  mytree <- drop.tip(tree, tip)
  
  H_max = ifelse(is.ultrametric(mytree), get.rooted.tree.height(mytree), max(ape::node.depth.edgelength(mytree)))
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
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

# check.q -------------------------------------------------------------------
check.q <- function(q) {
  
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  if(min(q) < 0){
    warning("ambigous of order q, we only compute postive q", call. = FALSE)
    q <- q[q >= 0]
  }
  
  return(q)
}

# check.conf -------------------------------------------------------------------
check.conf <- function(conf) {
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  return(conf)
}

# check.nboot -------------------------------------------------------------------
check.nboot <- function(nboot) {
  
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  return(nboot)
}

# check.base -------------------------------------------------------------------
check.base <- function(base) {
  
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE))) stop("invalid datatype")
  if (pmatch(base, BASE) == -1) stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  
  return(base)
}

# check.PDtype -------------------------------------------------------------------
check.PDtype <- function(type) {
  
  divtype <- c("PD", "meanPD")
  if(is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.", call. = FALSE)
  
  return(type)
}

# check.size -------------------------------------------------------------------
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
      unique(mi)
    })
    
  } else {
    
    if (class(size) == "numeric" | class(size) == "integer" | class(size) == "double") {
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
        
      unique(mi)
    })
  }
  
  return(size)
}

# check.level -------------------------------------------------------------------
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
  
  return(level)
}


