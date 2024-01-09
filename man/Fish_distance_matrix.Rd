\name{Fish_distance_matrix}
\alias{Fish_distance_matrix}
\docType{data}
\title{
Species pairwise distance matrix for species in the dataset \code{Fish_incidence_data}
}
\description{
  This dataset includes a species pairwise symmetric distance matrix for 88 species listed in the \code{Fish_incidence_data} dataset. Each element in the matrix represents a Gower distance between two species computed from species traits. \cr
}
\usage{data(Fish_distance_matrix)}
\format{
  \code{Fish_distance_matrix} is a 88 x 88 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one. \cr
  $ Agonus_cataphractus: num 0 0.326 0.280 0.055 0.093 ... \cr
  $ Alosa_fallax: num 0.326 0 0.380 0.336 0.311 ... \cr
  $ Ammodytes_marinus: num 0.280 0.380 0 0.243 0.364 ... \cr
  ....... \cr
 }
\keyword{datasets}
