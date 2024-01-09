\name{Brazil_rainforest_distance_matrix}
\alias{Brazil_rainforest_distance_matrix}
\docType{data}
\title{
Species pairwise distance matrix for species in the dataset \code{Brazil_rainforest_abun_data}
}
\description{
  This dataset includes a tree species pairwise symmetric distance matrix for 425 species listed in the \code{Brazil_rainforest_abun_data} dataset. Each element in the matrix represents a Gower distance between two species computed from species traits. \cr
}
\usage{data(Brazil_rainforest_distance_matrix)}
\format{
  \code{Brazil_rainforest_distance_matrix} is a 425 x 425 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one. \cr
  $ Carpotroche_brasiliensis: num 0 0.522 0.522 0.253 0.263 ... \cr
  $ Astronium_concinnum: num 0.522 0 0 0.525 0.509 ... \cr
  $ Astronium_graveolens: num 0.522 0 0 0.525 0.509 ... \cr
  ....... \cr
 }

\keyword{datasets}
