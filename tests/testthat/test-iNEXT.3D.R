context("iNEXT3D")
test_that("iNEXT3D for abundance-based data", {
  # Test input by a demo data
  data("Brazil_rainforest_data")
  out <- iNEXT3D(Brazil_rainforest_data$data, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "n")
  expect_equal(nrow(out$TDInfo), length(Brazil_rainforest_data$data))
  
  # Test input by a vector
  x <- Brazil_rainforest_data$data$Edge
  out <- iNEXT3D(x, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "n")
  expect_equal(nrow(out$TDInfo), 1)
  
  # Test input by a data.frame
  x <- data.frame(a = c(10,20,30,40,50,0,0), b = c(11,22,0,0,33,44,55))
  out <- iNEXT3D(x, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "n")
  expect_equal(nrow(out$TDInfo), ncol(x))
  
})

test_that("iNEXT3D for sampling-unit-based incidence frequencies data", {
  # Test input by a demo data
  data("fish_incidence_data")
  out <- iNEXT3D(fish_incidence_data$data, q = 0, datatype = "incidence_raw")
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "T")
  expect_equal(nrow(out$TDInfo), length(fish_incidence_data$data))
  
  # # Test input by a vector
  # out <- iNEXT3D(c(ncol(fish_incidence_data$data$`2013-2015`), rowSums(fish_incidence_data$data$`2013-2015`)), q = 0, datatype = "incidence_freq")
  # expect_is(out, "iNEXT3D")
  # expect_output(str(out), "List of 3")
  # expect_equal(names(out$TDInfo)[2], "T")
  # expect_equal(nrow(out$TDInfo), 1)
})


# test_that("iNEXT3D for species by sampling-units incidence matrix", {
#   # Test input by a demo data
#   data(data.inc)
#   options(warn=-1)
#   out <- iNEXT3D(data.inc$data, q = 0, datatype = "incidence_raw", nT = data.inc$nT)
#   expect_is(out, "iNEXT3D")
#   expect_output(str(out), "List of 3")
#   expect_equal(names(out$TDInfo)[2], "T")
#   expect_equal(nrow(out$TDInfo), length(data.inc$nT))
#   
#   # Test input by a data.frame
#   # x <- ciliates$EtoshaPan
#   # # expect_equal(class(x), "matrix")
#   # out <- iNEXT3D(x, q=0, datatype="incidence_raw")
#   # expect_is(out, "iNEXT3D")
#   # expect_output(str(out), "List of 3")
#   # expect_equal(names(out$TDInfo)[2], "T")
#   # expect_equal(nrow(out$TDInfo), 1)
# })


