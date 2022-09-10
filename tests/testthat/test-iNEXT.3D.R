context("iNEXT.3D")
test_that("iNEXT.3D for abundance-based data", {
  # Test input by a demo data
  data(dunes)
  out <- iNEXT3D(dunes$data, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), length(dunes))
  
  # Test input by a vector
  x <- dunes$data$EM
  out <- iNEXT3D(x, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), 1)
  
  # Test input by a data.frame
  x <- data.frame(a = c(10,20,30,40,50,0,0), b = c(11,22,0,0,33,44,55))
  out <- iNEXT3D(x, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), ncol(x))
  
})

test_that("iNEXT.3D for sampling-unit-based incidence frequencies data", {
  # Test input by a demo data
  data(ant)
  out <- iNEXT3D(ant, q = 0, datatype = "incidence_freq")
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), length(ant))
  
  # Test input by a vector
  out <- iNEXT3D(ant$h50m, q = 0, datatype = "incidence_freq")
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), 1)
})


test_that("iNEXT.3D for species by sampling-units incidence matrix", {
  # Test input by a demo data
  data(data.inc)
  options(warn=-1)
  out <- iNEXT3D(data.inc$data, q = 0, datatype = "incidence_raw", nT = data.inc$nT)
  expect_is(out, "iNEXT.3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), length(data.inc$nT))
  
  # Test input by a data.frame
  # x <- ciliates$EtoshaPan
  # # expect_equal(class(x), "matrix")
  # out <- iNEXT.3D(x, q=0, datatype="incidence_raw")
  # expect_is(out, "iNEXT.3D")
  # expect_output(str(out), "List of 3")
  # expect_equal(names(out$DataInfo)[2], "T")
  # expect_equal(nrow(out$DataInfo), 1)
})


