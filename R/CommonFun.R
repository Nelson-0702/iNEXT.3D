# as.incfreq -------------------------------------------------------------------
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
      names(mydata) <- paste0("assemblage",1:length(nT))
    }else{
      names(mydata) = names(nT)
    }
    data = lapply(mydata, function(i){
      out = c('nT' = ncol(i), rowSums(i))
      return(out)
    })
  } else if (inherits(data, "list")) 
    data <- lapply(data, function(i) c('nT' = ncol(i), rowSums(i)))
  if (length(data) == 1) data = data[[1]]
  return(data)
}


# Coverage -------------------------------------------------------------------
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
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  Sub2 <- function(m){
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/u*A
    if(m > n) out <- 1-f1/u*A^(m-n+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype!='abundance', Sub2(i), Sub(i) )
  })
}


