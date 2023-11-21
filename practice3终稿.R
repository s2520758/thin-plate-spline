#Yingjie Yang 
#s2520758
#Obtain consistent random results so that each plot is the same
set.seed(10)
ff <- function(x) {exp(-(x[,1]-.3)^2/.2^2-(x[,2] - .3)^2/.3^2)*.5 +
    exp(-(x[,1]-.7)^2/.25^2 - (x[,2] - .8 )^2/.3^2)}
n <- 500
x <- matrix(runif(n*2),n,2)
y <- ff(x) + rnorm(n)*.1 ## generate example data to fit
getTPS <- function(x, y, k = 100) {
  n <- nrow(x)
  
  # Selecting x* points
  if (k < n) {
    xk <- x[sample(1:n, size = k), , drop = FALSE]#if k<n,pick random k from x
  } else {
    
    xk<- x
  }
  
  # Define the function eta
  eta <- function(r) {
    ifelse(r > 0, r^2 * log(r), 0)
  }
  
  # Generate the E_star and E matrices
  
  # Generate the E_star matrix (which will be k x k)
  E_star <- matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      E_star[i, j] <- eta(norm(as.matrix(xk[i,] - xk[j,]), type = "F"))
    }
  }
  
  # Generate E matrix (which will be n x k)
  E <- matrix(0, n, k)
  for (i in 1:n) {
    for (j in 1:k) {
      E[i, j] <- eta(norm(as.matrix(x[i,] - xk[j,]), type = "F"))
    }
  }
  
  # Generate the T matrix (which will be n x 3)
  T <- cbind(rep(1, n), x)
  # Generate the T_star matrix (which will be k x 3)
  T_star<- cbind(rep(1, k), xk)
  
  # Extracting Z 
  Z <- qr.Q(qr(T_star),complete=TRUE)[,-(1:3)]
  
  # Creating the X matrix
  X <- cbind(E %*% Z, T)
  
  # Creating S matrix(which will be k by k) 
  S <- matrix(0, nrow=k, ncol=k)
  S[1:(k-3),1:(k-3)] <-(t(Z) %*% E_star %*% Z)
  
  #return the list
  return(list(xk = xk, X = X, S = S,Z = Z,E=E,T=T))
}


fitTPS <- function(x, y, k=100, lsp=c(-5,5)) {
  # Using getTPS to get X,S,T
  tps_result <- getTPS(x, k) 
  X <- tps_result$X
  S <- tps_result$S
  T <- tps_result$T
  
  # Extract Q matrix
  Q <- qr.Q(qr(X))
  # Extract R matrix
  R <- qr.R(qr(X))
  # Creating lambda values 
  lambda_values <- exp(seq(lsp[1], lsp[2], length.out=100))
  #  make t(solve(R)) %*% S %*% solve(R) be symmetric
  A<- t(solve(R)) %*% S %*% solve(R)
  A <- (t(A)+A)*.5
  #Using eigen decomposition from A to get U and Lambda
  eigen_result <- eigen( A)
  U <- eigen_result$vectors
  Lambda <- diag(eigen_result$values)
  #Creating identity matrix
  I <- diag(length(diag(Lambda)))
  
  # Initialize gcv_scores and EDFs
  gcv_scores <- numeric(100)
  EDFs <- numeric(100)
  #Using loop to find the index of optimal lambda
  for (i in 1:100) {
    lambda <- lambda_values[i]
    
    # Calculate  EDF     
    EDF <- sum(diag(I + lambda * Lambda)^-1)
    # Calculate  beta_hat 
    beta_hat= solve(R) %*% U %*% solve(I + lambda * Lambda) %*% t(U)%*% t(Q) %*% y
    # Calculate  mu
    mu <- X %*% beta_hat
    # Calculate gcv_scores for each   lambda
    gcv_scores[i] <- sum((y - mu)^2) / (n - EDF)^2
    EDFs[i] <- EDF
  }
  
  # Find the index of optimal lambda 
  best_lambda_index <- which.min(gcv_scores)
  #Find the optimal lambda 
  optimal_lambda<-lambda_values[best_lambda_index]
  #Find beta at optimal lambda
  beta <- solve(R) %*% U %*% solve(I + optimal_lambda * Lambda) %*% t(U)%*% t(Q) %*% y
  #Find mu at optimal lambda
  mu <- X %*% beta
  #return an object of class tps
  result<- list(beta = beta, mu = mu,medf = EDFs[best_lambda_index],lambda = lambda_values,
                gcv = gcv_scores,edfs = EDFs)
  class(result) <- "tps"
  return(result)
}


plot.tps <- function(result, x, k=100) {
  result<-fitTPS (x, y, k, lsp=c(-5,5))
  #Geting   beta from fitTPS
  beta <- result$beta
  #Geting alpha from beta 
  alpha <- beta[(nrow(beta) - 2):nrow(beta), , drop = FALSE]
  #Geting delta_z from beta
  delta_z <- beta[1:(nrow(beta) - 3), , drop = FALSE]
  #Geting   Z,xk  from gettTPS
  Z <- getTPS(x, k)$Z
  xk <- getTPS(x, k)$xk
  #Geting delta according to formula
  delta <- Z %*% delta_z
  #Creating sample size
  m=50
  #Extract m numbers from each column of x to form x1, x2 respectively 
  x1 <- seq(min(x[, 1]), max(x[, 1]),length=m)
  x2 <- seq(min(x[,2]), max(x[,2]),length=m)
  xp<-cbind(x1,x2)
  #add"1" to make sure later a0 could multiply 1
  Tp<-cbind(rep(1,m),xp)
  #Redefine eta
  eta <- function(r) {
    ifelse(r > 0, r^2 * log(r), 0)
  }
  #Calculate Ep
  Ep <- matrix(0, nrow(xp), k)
  for (i in 1:nrow(xp)) {
    for (j in 1:k) {
      Ep[i, j] <- eta(norm(as.matrix(xp[i,] - xk[j,]), type = "F"))
    }
  }
  #according to formula calculate f(x)
  fx <- matrix(0, m, m)
  for (i in 1: m) {
    for (j in 1: m) {
      fx[i, j] <- Tp[i, ] %*% alpha + Ep[i, ] %*% delta + Ep[j, ] %*% delta
    }
  }
  
  persp(x1, x2, matrix(fx,m,m), theta = 30, phi = 30)
}
#call plot.tps
plot.tps (result, x, k=100) 

