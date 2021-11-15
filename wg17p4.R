require(debug)

finite_diff <- function(theta, f) {
  # Purpose: Find gradient vector of f at point theta using finite differentiation.
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
    f_dtheta[i] <- (f(th1) - f(theta))/eps ## approximate -df/dth[i]
  }
  
  f_dtheta
}

get_gradient <- function(f_theta) {
  # Purpose: Get the gradient of function f at point theta.
  
  # Do we need to calculate the gradient?
  f_grad <- attr(f_theta, "gradient")
  if(is.null(f_grad)) {
    f_dtheta <- finite_diff(theta, f)
  } else {
    f_dtheta <- f_grad
  }
  
  f_dtheta
}

has_converged <- function(f_val, f_grad, fscale, tol) {
  val <- max(abs(f_grad)) < (abs(f_val)+fscale)*tol
  if (is.na(val)) {
    FALSE
  } else {
    val
  }
}

second_wolfe <- function(delta, f, theta) {
  theta_delta <- theta + delta
  theta_delta_grad <- finite_diff(theta_delta, f)
  left_side <- t(theta_delta_grad) %*% delta
  right_side <- 0.9* t(finite_diff(theta, f)) %*% delta
  
  left_side >= right_side
}

reduces_obj <- function(f, theta, step) {
  if(f(theta) > f(theta + step)) {
    TRUE
  }
  FALSE
}

enqueue <- function(list, item) {
  list <- append(item, list)
  list
}

dequeue <- function(list) {
  first <- list[[1]]
  list <- list[-1]
  list(first, list)
}

bfs <- function(step_v, f, theta) {
  Q <- list(step_v)
  # Maybe set val to infinity?
  val <- step_v

  while(length(Q) != 0) {
    q_v <- dequeue(Q)
    val <- q_v[[1]]
    Q <- q_v[[2]]

    # Sufficient conditions for good step length    
    if(second_wolfe(val, f, theta) | reduces_obj(f, theta, val)) {
      break
    }
    
    # Since this is all based on Taylor, we do not stray too far from the original step
    next_up <- val * 1.5
    next_down <- val / 2
    Q <- enqueue(Q, next_up)
    Q <- enqueue(Q, next_down)
  }
  val
}

get_step <- function(B, f, theta) {
  f_theta <- f(theta)
  grad_theta0 <- get_gradient(f, theta, f_theta)
  step <- drop(-B %*% grad_theta0)
  step
  
  #Perform a bfs search on [theta, theta +step] to find a step length
  #step <- bfs(step, f, theta)
  
  #step
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
  iter <- 0

  while (iter < maxit) {
    iter <- iter + 1
    
    f_theta0 <- f(theta0, ...)
    grad_theta0 <- get_gradient(f_theta0)

    # Quasi Newton step from theta0 to theta1 
    delta <- drop(-B %*% grad_theta0)

    # New theta value 
    theta1 <- theta0 + delta
    print('new theta val')
    print(theta1)
    
    f_theta1 <- f(theta1, ...)
    grad_theta1 <- get_gradient(f_theta1)
    print('grad vec of f at new theta')
    print(grad_theta1)

    # Initial s and y vectors. Can this just be delta?
    step_k <- theta1 - theta0
    step_k_t <- t(step_k)
    y_k <- grad_theta1 - grad_theta0

    rho_inv_k <- drop(step_k_t %*% y_k)
    
    # rho_inv_k is a scalar
    rho_k <- 1/ rho_inv_k
    print('rho')
    print(rho_k)

    # Calculate next B matrix
    expression1 <- (diag(length(theta)) - rho_k*(step_k %*% t(y_k)))
    expression2 <- (diag(length(theta)) - rho_k*(y_k %*% step_k_t))
    expression3 <- rho_k * step_k %*% step_k_t
    # Is there some optimization to be done here by first calculating the rightmost matrix multiplication?
    B <- expression1 %*% B %*% expression2 + expression3
    print('new B')
    print(B)
    
    if(has_converged(f_theta0, grad_theta0, fscale, tol) == TRUE) {
      break
    }
    theta0 <- theta1
  }

  bfgs_res = list(
    f=f(theta0),
    theta=theta0,
    iter=iter,
    g=grad_theta0
  )
  bfgs_res
}

# Test function Simon provided (Used to check the Quasi method)
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ##

#mtrace(get_step)
bfgs(c(-1,2), rb, getg=TRUE, maxit = 40)
