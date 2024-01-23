context("iNEXT3D")
test_that("iNEXT3D for abundance-based data", {
  
  ## TD
  # Test input by a demo data
  data("Brazil_rainforest_abun_data")
  out <- iNEXT3D(Brazil_rainforest_abun_data$Interior, datatype = "abundance", q = 0, size = 1)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(nrow(out$TDInfo), 1)
  
  # Test input by a vector
  x <- Brazil_rainforest_abun_data$Edge
  out <- iNEXT3D(x, q = 0, datatype = "abundance")
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "n")
  expect_equal(nrow(out$TDInfo), 1)
  
  
  ## PD
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_phylo_tree")
  out <- iNEXT3D(Brazil_rainforest_abun_data, diversity = 'PD', datatype = "abundance", q = 1, size = 1, 
                 PDtree = Brazil_rainforest_phylo_tree)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(nrow(out$PDInfo), length(Brazil_rainforest_abun_data))
  expect_equal(nrow(out$PDAsyEst), 3*length(Brazil_rainforest_abun_data))
  
  
  ## FD (single tau)
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_distance_matrix")
  out <- iNEXT3D(Brazil_rainforest_abun_data, diversity = 'FD', datatype = "abundance", q = 2, size = 1, 
                 FDdistM = Brazil_rainforest_distance_matrix, nboot = 0, FDtype = "tau_values", FDtau = 0.1)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(nrow(out$FDInfo), length(Brazil_rainforest_abun_data))
  
  
  ## FD (AUC)
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_distance_matrix")
  out <- iNEXT3D(Brazil_rainforest_abun_data, diversity = 'FD', datatype = "abundance", q = 2, size = 1, 
                 FDdistM = Brazil_rainforest_distance_matrix, nboot = 0)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(nrow(out$FDInfo), length(Brazil_rainforest_abun_data))
})

test_that("iNEXT3D for sampling-unit-based incidence raw data", {
  
  ## TD
  # Test input by a demo data
  data("Fish_incidence_data")
  out <- iNEXT3D(Fish_incidence_data$`2013-2015`, q = 0, datatype = "incidence_raw", size = 2)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$TDInfo)[2], "T")
  expect_equal(names(out$TDiNextEst$size_based)[3], "mT")
  expect_equal(nrow(out$TDInfo), 1)
  
  
  ## PD
  data("Fish_incidence_data")
  data("Fish_phylo_tree")
  out <- iNEXT3D(Fish_incidence_data, diversity = 'PD', q = 1, datatype = "incidence_raw", 
                 size = 2, PDtree = Fish_phylo_tree, nboot = 0)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$PDInfo)[2], "T")
  expect_equal(names(out$PDiNextEst$size_based)[3], "mT")
  expect_equal(nrow(out$PDInfo), length(Fish_incidence_data))
  
  
  ## FD (single tau)
  data("Fish_incidence_data")
  data("Fish_distance_matrix")
  out <- iNEXT3D(Fish_incidence_data, diversity = 'FD', q = 1, datatype = "incidence_raw", size = 2,
                 FDdistM = Fish_distance_matrix, FDtype = "tau_values", FDtau = 0.1)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$FDInfo)[2], "T")
  expect_equal(names(out$FDiNextEst$size_based)[3], "mT")
  expect_equal(nrow(out$FDInfo), length(Fish_incidence_data))
  
  
  ## FD (AUC)
  data("Fish_incidence_data")
  data("Fish_distance_matrix")
  out <- iNEXT3D(Fish_incidence_data, diversity = 'FD', q = 2, datatype = "incidence_raw", size = 2,
                 FDdistM = Fish_distance_matrix, nboot = 0)
  expect_is(out, "iNEXT3D")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$FDInfo)[2], "T")
  expect_equal(names(out$FDiNextEst$size_based)[3], "mT")
  expect_equal(nrow(out$FDInfo), length(Fish_incidence_data))
})

test_that("ObsAsy3D for abundance-based data", {
  
  ## TD
  # Test input by a demo data
  data("Brazil_rainforest_abun_data")
  out <- ObsAsy3D(Brazil_rainforest_abun_data$Edge, diversity = 'TD', 
                  datatype = "abundance", q = 1)
  expect_equal(nrow(out), 2)
  
  ## PD
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_phylo_tree")
  out <- ObsAsy3D(Brazil_rainforest_abun_data, diversity = 'PD', q = 2, PDreftime = 2, 
                  datatype = "abundance", nboot = 10, PDtree = Brazil_rainforest_phylo_tree)
  expect_equal(nrow(out), 4)
  
  
  ## FD (single tau)
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_distance_matrix")
  out <- ObsAsy3D(Brazil_rainforest_abun_data, diversity = 'FD', q = 1, 
                  datatype = "abundance", nboot = 10, FDdistM = Brazil_rainforest_distance_matrix, 
                  FDtype = 'tau_values', FDtau = 0.1)
  expect_equal(nrow(out), 4)
  
  
  ## FD (AUC)
  data("Brazil_rainforest_abun_data")
  data("Brazil_rainforest_distance_matrix")
  out <- ObsAsy3D(Brazil_rainforest_abun_data, diversity = 'FD', q = 2, 
                  datatype = "abundance", nboot = 10, FDdistM = Brazil_rainforest_distance_matrix, 
                  FDtype = 'AUC', FDcut_number = 30, method = "Observed")
  expect_equal(nrow(out), 2)
})

test_that("ObsAsy3D for sampling-unit-based incidence raw data", {
  
  ## TD
  # Test input by a demo data
  data("Fish_incidence_data")
  out <- ObsAsy3D(Fish_incidence_data$`2016-2018`, diversity = 'TD', q = 0,
                  datatype = "incidence_raw")
  expect_equal(nrow(out), 2)
  
  
  ## PD
  data("Fish_incidence_data")
  data("Fish_phylo_tree")
  out <- ObsAsy3D(Fish_incidence_data, diversity = 'PD', q = 1, datatype = "incidence_raw", 
                  PDtree = Fish_phylo_tree, nboot = 0)
  expect_equal(nrow(out), 4)
  
  
  ## FD (single tau)
  data("Fish_incidence_data")
  data("Fish_distance_matrix")
  out <- ObsAsy3D(Fish_incidence_data, diversity = 'FD', q = 1, datatype = "incidence_raw", 
                  FDdistM = Fish_distance_matrix, FDtype = "tau_values", FDtau = 0.1)
  expect_equal(nrow(out), 4)
  
  
  ## FD (AUC)
  data("Fish_incidence_data")
  data("Fish_distance_matrix")
  out <- ObsAsy3D(Fish_incidence_data, diversity = 'FD', q = 2, datatype = "incidence_raw", 
                  FDdistM = Fish_distance_matrix, FDtype = "AUC", nboot = 20, method = "Observed")
  expect_equal(nrow(out), 2)
})

