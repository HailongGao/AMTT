####################################generate sample############################################# 
######banded covariance

banded_covariance<-function(p){
  A <- matrix(0, p, p)
  for (i in 1:p) {
    A[i, i] <- 1
    if (i > 1) A[i, i - 1] <- 0.3
    if (i < p) A[i, i + 1] <- 0.3
  }
  D<-diag(runif(p,1,3))
  D_sqrt <- eigen(D)$vectors %*% diag(sqrt(eigen(D)$values)) %*% t(eigen(D)$vectors)
  sigMat <- D_sqrt %*% A %*% D_sqrt
  return(sigMat)
}

#######autoregressive covariance
auto_covariance<-function(p){
  R <- diag(p)  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      R[i, j] <- 0.5^abs(i - j)  
      R[j, i] <- R[i, j]  
    }
  }
  
  D<-diag(p)
  D_sqrt <- eigen(D)$vectors %*% diag(sqrt(eigen(D)$values)) %*% t(eigen(D)$vectors)
  sigMat <- D_sqrt %*% R %*% D_sqrt
  return(sigMat)
}
#####sparse covariance
Psi_covariance<-function(p){
  q <- round(3 * p^(1/2))
  generate_diagonal_matrix <- function(B) {
    epsilon <- max(-eigen(B, only.values = TRUE)$values, 0) + 0.05
    A1 <- B + diag(epsilon, nrow = q)
    A2 <- diag(p - q)
    return(list(A1 = A1, A2 = A2))
  }
  
  B <- matrix(0, nrow = q, ncol = q)
  for (i in 1:nrow(B)) {
    for (j in 1:ncol(B)) {
      if (j <= i) {
        rand <- runif(1)
        if (rand <= 0.5) {
          B[i, j] <- 0
        } else {
          B[i, j] <- sample(c(runif(1, -1, -0.5), runif(1, 0.5, 1)), 1)
        }
      }
    }
  }
  B <- B + t(B) - diag(diag(B))
  A_matrices <- generate_diagonal_matrix(B)
  A1 <- A_matrices$A1
  A2 <- A_matrices$A2
  Psi <- matrix(0, nrow = p, ncol = p)
  Psi[1:nrow(A1), 1:ncol(A1)] <- A1
  Psi[(nrow(A1)+1):(nrow(A1)+nrow(A2)), (ncol(A1)+1):(ncol(A1)+ncol(A2))] <- A2
  return(Psi=Psi)
}

#######
Mu<-function(p, beta, r,n1,n2) {
  # input: p@Dimension of mu
  #        beta@Number of difference dimensions between mu (p^beta)
  # output: mu1, mu2    
  
  b <- round(p^beta)
  signal <-sample(c(1, -1), size = b, replace = TRUE, prob = c(0.5, 0.5))*rep((2*r*(1/n1+1/n2)*log(p))^(1/2),b)
  ind<-sample(1:p, length(signal), replace = FALSE)
  mu1 = rep(0,p)
  
  if (r == 0) {
    mu2 <- mu1
  } else {
    mu2=mu1
    mu2[ind] = mu1[ind] + signal
  }
  return(list(mu1 = mu1, mu2 = mu2))
}

######Multtivariate normal distribution  
normalZ <- function(n1, n2,mu1,mu2,sigmat) {
  # input: p@Dimension of OTU
  #        s@Number of difference dimensions between OTU (1,2,3)->(0.05p,0.1p,0.5p)
  #        a@Control hypothesis
  #        n1, n2@Sample size   
  #  
  # output: Z1, Z2     
  Z1 <-t(MASS::mvrnorm(n1, mu = mu1, Sigma = sigmat))
  Z2 <-t(MASS::mvrnorm(n2, mu = mu2, Sigma = sigmat))
  return(list(Z1 = Z1, Z2 = Z2))
}