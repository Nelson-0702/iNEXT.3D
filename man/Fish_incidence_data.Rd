\name{Fish_incidence_data}
\alias{Fish_incidence_data}
\docType{data}
\title{
Incidence data (datatype = "incidence_raw")
}
\description{
  This dataset includes incidence-based occurrence records of estuarine fishes collected at Bridgwater Bay in UK’s Bristol Channel (Chao et al. 2021). There are two time periods (“2013-2015” and “2016-2018”). Each time period is designated as an assemblage. The purpose was to compare 3D diversity of the two time periods. In each time period, species incidence was recorded in 36 sampling units in each assemblage; each sampling unit represents a sampling date. This dataset is a list of two data frames with each representing species-by-sampling-unit raw incidence (detection/non-detection) matrix. \cr
}
\usage{data(Fish_incidence_data)}
\format{
  \code{Fish_incidence_data} is a list of two species-by-sampling-unit data frames for 88 species and 36 sampling units. Each element in the data frame is 1 for a detection, and 0 for a non-detection. \cr
  A list of 2: \cr
  $ 2013-2015: int [1:88, 1:36] 0 0 0 0 0 0 0 0 0 0 ... \cr
  $ 2016-2018: int [1:88, 1:36] 1 0 0 0 0 0 0 0 0 0 ... \cr
 }
\references{
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M. and Magurran, A.
E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.
}
 
\keyword{datasets}
