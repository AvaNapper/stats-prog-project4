require(debug)

finite_diff <- function(theta, f) {
  # Purpose: Find gradiant vector of f at point theta using finite differentiation.
  # Input:  theta - Vector, the location where we want to differentiate f
  #         f     - Function to differentiate
  
  # Output: The value of f' at point theta.

  #TODO: Should our initial value for the delta be this or 0?
  f_dtheta <- theta
  
  # https://stackoverflow.com/questions/2619543/how-do-i-obtain-the-machine-epsilon-in-r
  # For nice functions, the square root of machine precision yields the optimum value for epsilon
  # Too big and our gradient is incorrect, too small and the computer cannot store the value.
  eps <- sqrt(.Machine$double.eps)
  
  for (i in 1:length(theta)) { ## loop over parameters/dimensions
    th1 <- theta
    th1[i] <- th1[i] + eps ## move eps in dimension i
    f_dtheta[i] <- (f(th1, ...) - f(theta, ...))/eps ## approximate -df/dth[i]
  }
  
  f_dtheta
}

bfgs <- function(theta, f, ..., tol=1e-5, fscale=1, maxit=100) {
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

  # Initial B matrix
  B0 <- diag(length(theta0))
  B <- B0

  for (j in 1:maxit){

    f_dtheta0 <- finite_diff(theta0, f)

    # Quasi Newton step from theta0 to theta1 
    delta <- drop(-B %*% f_dtheta0)


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
    sk <- theta1 - theta0
    yk <- f_dtheta1 - f_dtheta0

  
    # Initial rho vector
    rho_inv_k <- drop(t(sk) %*% yk)
    rho_k <- 1/ rho_inv_k
    print('intial rho')
    print(rho_k)


    # To calc next B matrix
    expression1 <- (diag(length(theta)) - rho_k*(sk %*% t(yk)))
    expression2 <- (diag(length(theta)) - rho_k*(yk %*% t(sk)))
    expression3 <- rho_k * sk %*% t(sk)
    B <- expression1 %*% B %*% expression2 + expression3
    print('new B')
    print(B)

    theta0 <- theta1
  }



  #outputs <- list(f(theta1), theta1, iter, g, H)
  print('Final theta')
  print(theta1)
  #return(outputs)
  
}

# Test function Simon provided (Used to check the Quasi method)
rb <- function(theta,getg=FALSE,k=100) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ##


bfgs(initial, rb, maxit = 40)