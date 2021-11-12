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



  # Finite differencing to find gradiant vector of f
  fd <- th0 <- theta
  eps <- 1e-7
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
    fd[i] <- (f(th1) - f(th0))/eps ## approximate -dl/dth[i]
  }
  
  return(f, theta, iter, g, H)
  
}