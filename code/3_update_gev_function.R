# update function
if(!require(Deriv)){install.packages("Deriv")};library(Deriv)
if(!require(evd)){install.packages("evd")};library(evd)

#loss function
# logl=expression(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))
# func_loss = function()

#z_function
func_z <- function(dmatrix,mu,u,lam,rho){
  z = ifelse(abs(dmatrix %*% mu + u) > (lam/rho) ,
             (dmatrix %*% mu) + u - sign(u + (dmatrix %*% mu)) * (lam /rho), 0)
  return(z)
}


#u_function
func_u <- function(dmatrix,mu,z,u) {
  u <- u + (dmatrix %*% mu) - z
  return(u)
}

#dmatrix
mat_func = function(n) {
  m = matrix(0,n-1,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-1)) {
    m[i,i] = -1
    m[i,i+1] = 1
  }
  return(m)
}

## gevfunction 
gevreg = function(x, z, ctr_list)
{
  
  l2gev = function (tvec, dmatrix, rho, z_init, u_init)
  {
    # loc.vec = zz%*%tvec[1:(p+1)]
    loc.vec = tvec[1:n]
    # loglikelihood
    v1 = - sum(lgev(x, loc = loc.vec, scale = tvec[n+1], shape = tvec[n+2]))
    
    # lagrangian term
    v2 = (rho/2)*sum(((dmatrix %*% loc.vec) - z_init + u_init)^2)
    v = v1 + v2
    return(v)
  }
  
  lgev = function (x, loc = 0, scale = 1, shape = 0) 
  {
    if (min(scale) <= 0) 
      return( - 1e+6)
    if (length(shape) != 1) 
      stop("invalid shape")
    x <- (x - loc)/scale
    if (shape == 0) 
      d <- log(1/scale) - x - exp(-x)
    else {
      nn <- length(x)
      xx <- 1 + shape * x
      xxpos <- xx[xx > 0 | is.na(xx)]
      scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
      d <- numeric(nn)
      d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
        (1/shape + 1) * log(xxpos)
      d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    }
    return(d)
  }
  
  est = optim(tvec, l2gev, dmatrix = dmatrix,
              rho=init_rho, z_init = z_init, u_init = u_init, method = "BFGS",
              control = ctr_list)$par
  
  
  return(est)
}
## refit function

gevreg_m = function(tvec = tvec, x = x, zmat = zmat)
{
  p = ncol(zmat)
  l2gev_m = function (tvec = tvec, x = x, zmat = zmat)
  {
    # ns = length(xlist)
    # v1 = 0
    # for ( i in 1: length(xlist))
    # {
    loc.vec.reg = drop(zmat%*%tail(tvec, p))
    v = -sum(lgev(x, loc = loc.vec.reg, # v1이 필요없음.
                  scale = tvec[1], shape = tvec[2]))   # loss : +(-loglikelihood)
    # }
    # loglikelihood
    # regularization
    
    
    return(v)
  }
  
  lgev = function (x, loc = 0, scale = 1, shape = 0) 
  {
    if (min(scale) <= 0) 
      return( - 1e+6)
    if (length(shape) != 1) 
      stop("invalid shape")
    x <- (x - loc)/scale
    if (shape == 0) 
      d <- log(1/scale) - x - exp(-x)
    else {
      nn <- length(x)
      xx <- 1 + shape * x
      xxpos <- xx[xx > 0 | is.na(xx)]
      scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
      d <- numeric(nn)
      d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
        (1/shape + 1) * log(xxpos)
      d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    }
    return(d)
  }
  
  # tmp_loc = 0
  # for ( i in 1:ns)
  # {
  #   x = xlist[[i]]
  #   fit = fgev(x)
  #   tmp_loc = tmp_loc + fit$estimate[1]
  #   tvec[2*i] = fit$estimate[2]
  #   tvec[2*i+1] = fit$estimate[3]
  # }
  # tmp_loc = tmp_loc/ns
  # tvec[1] = tmp_loc
  
  fit = optim(tvec, l2gev_m, method = "BFGS", 
              x = x, zmat = zmat)
  return(fit) 
}
lgev = function (x, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    return( - 1e+6)
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x - exp(-x)
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -(1e+6)
  }
  return(d)
}