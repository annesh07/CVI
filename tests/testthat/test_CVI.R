test_that("Number of clusters from latent allocation:",{
  expect_equal(CVI(N = 50, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01, L20 = 0.01,
   X = matrix(1, nrow = 50, ncol = 2), W1 = 0.01, W2 = 0.01,
   L1 = matrix(0.01, nrow = 10, ncol = 2), L2 = matrix(0.01, nrow = 10, ncol = 1
   ), Plog = matrix(-3, nrow = 50, ncol = 10), maxit = 1000)$Clusters,1)})
