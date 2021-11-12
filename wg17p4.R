bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
    # Purpose:
    # Input:  theta - vector of intial values for optimization parameters
    #         f     - objective function to minimise
    #                 first argument is vector of optimization parameters
    #                 second argument is logical, should the gradients of the 
    #                   objective function w.r.t. parameters be computed
    #         tol   - convergence tolerance
    #         maxit - max number of BFGS iterations to try
    #
    #
    # Output: f     - scalar value of objective func at the minimum
    #         theta - vector of values of parameters at minimum
    #         iter  - number of iterations it took to reach maximum
    #         g     - gradient vector at the minimum
    #         H     - approx Hessian matrix at the minimum

  # Set input theta as our initial theta
  theta0 <- theta
  print('Initial theta')
  print(theta0)


  # Finite differencing to find gradiant vector of f
  f_dtheta0 <- th0 <- theta0
  eps <- 1e-7 ## step length for finite differencing
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
    f_dtheta0[i] <- (f(th1) - f(th0))/eps ## approximate -dl/dth[i]
  }
  print('grad vector for f at intial theta')
  print(f_dtheta0)


  # Initial B matrix
  B0 <- diag(length(theta))
  
  # Quasi Newton step from theta0 to theta1 
  delta <- drop(-B0 %*% f_dtheta0)
  print('HEREEE')
  print(theta0)
  print(delta)


  # New theta value 
  theta1 <- theta0 + delta
  print('new theta val')
  print(theta1)


  # Calculating gradient vector for f at new value of theta
  f_dtheta1 <- th0 <- theta1
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
    f_dtheta1[i] <- (f(th1) - f(th0))/eps ## approximate -dl/dth[i]
  }
  print('grad vec of f at new theta')
  print(f_dtheta1)



  # Initial s and y vectors
  s0 <- theta1 - theta0
  y0 <- f_dtheta1 - f_dtheta0
  print('initial s')
  print(s0)
  print('initial y')
  print(y0)

  
  # Initial rho vector
  rho_inv_0 <- drop(t(s0) %*% y0)
  rho_0 <- 1/ rho_inv_0
  print('intial rho')
  print(rho_0)


  # To calc next B matrix
  expression1 <- (diag(length(theta)) - rho_0*(s0 %*% t(y0)))
  expression2 <- (diag(length(theta)) - rho_0*(y0 %*% t(s0)))
  expression3 <- rho_0 * s0 %*% t(s0)
  B1 <- expression1 %*% expression2 %*% expression3
  print('new B')
  print(B1)



  #outputs <- list(f, theta, iter, g, H)
  #return(outputs)
  
}

g <- function(theta){
  return(theta[1]^2 + 2*theta[2])
}

g(c(1,2))
initial <- c(1,2)
bfgs(initial, g)