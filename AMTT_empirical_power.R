#input:      p@Dimension of OTU
#       n1, n2@Sample size
#         nsim@Recommended number of replications
#            B@Recommended permutation times
#            t@Recommended truncated parameter
#            r@signal strength parameter(𝑟 = 0, represents Type I errors, while any 𝑟 ≠ 0 represents statistical powers)
#         beta@Number of difference dimensions between mu 
#       sigmax@covariance structures
#
#
#
# output: alpha1@Type I error rates or power under significance level of α = 0.05
#         alpha2@Type I error rates or power under significance level of α = 0.01




p <- 100
n1 <- 100
n2 <- 100
n <- n1 + n2
nsim <- 10
B = 500
t = c(1 / p, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
palpha <- rep(NA, nsim)
  for (j in 1:nsim) {
    r <- 0
    beta <- 0.1
    sigmax <- Psi_covariance(p)
    mu <- Mu(p, beta, r, n1, n2)
    mu1 <- mu$mu1
    mu2 <- mu$mu2
    Z <- normalZ(n1, n2, mu1, mu2, sigmax)
    Z1 <- Z$Z1
    Z2 <- Z$Z2
    W1 = exp(t(Z1))
    W2 = exp(t(Z2))
    X1 = W1 / matrix(rowSums(W1), nrow = n1, ncol = p, byrow = FALSE)
    X2 = W2 / matrix(rowSums(W2), nrow = n2, ncol = p, byrow = FALSE)
    Ip <- diag(p)
    ones_column_vector <- matrix(1, nrow = p, ncol = 1)
    X <- as.matrix(rbind(X1, X2))
    Y <- t((Ip - (ones_column_vector %*% t(ones_column_vector)) / p) %*% t(log(X)))
    Y1_log <- Y[1:n1, ]
    Y2_log <- Y[(n1 + 1):n, ]
    AMTT <- amidat_parameterization(Y1_log, Y2_log,B, t )
    palpha[j] <-AMTT$P_final
    print(j)
  }

alpha1<- length(which(palpha< 0.05)) / nsim
alpha2<- length(which(palpha< 0.01)) / nsim
alpha1
alpha2


