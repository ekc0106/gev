# rm(list = ls())

if(!require(tidyverse)){install.packages("tidyverse")};require(tidyverse)
if(!require(ggplot2)){install.packages("ggplot2")};require(ggplot2)


load("Pr_46.RData")
source("scad_gev_function.R")

st = c(114,138,211,247,285)

ctr_list = list()
ctr_list$maxit = 20
ctr_list$reltol = 1e-6

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# resolution in ADMM
epri = 1e-3
edul = 1e-3
# lamb = c(0.01,0.02,0.09,0.3,0.7,2,5)
# lamb = c(0,0.03,0.05,0.07,0.1,0.5,1)
# lamb = 0.04
# lamb = c(0.5,1)
# lamb= c(0,0.5)
# lamb = c(0,0.01,0.02,0.03,0.05,0.07,0.09,0.1,0.3)
a = 3.7 ; rho = 0.5
lamb <- c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15 , 1)
# lamb <- seq(0.01,0.2, by = 0.001)

# lamb = c(0.05,0.15,0.2,0.25)
# re_result2 = list()
# h = 1

results_list <- vector(mode = 'list', length = length(st))
results_list <- lapply(1:length(results_list), function(x) vector(mode = 'list',length = length(lamb)))
names(results_list) <- paste('st',st,sep = '_')
for(i in 1:length(results_list)){
  names(results_list[[i]]) <- paste('lamb', lamb, sep = '_')
}

for(i in 1:length(st)){
  for(j in 1:length(lamb)){
    # n_st <- table(Pr_46$stnlds)[i]
    results_list[[i]][[j]]$x_mu <- data.frame(x = rep(0,46), mu = rep(0,46)) %>% as_tibble()
    results_list[[i]][[j]]$tmp_z_u <- data.frame(tmp_z = rep(0,45), tmp_u = rep(0,45)) %>% as_tibble()
    results_list[[i]][[j]]$theta <- numeric(2)
    results_list[[i]][[j]]$df_aic_bic <- numeric(3) 
    results_list[[i]][[j]]$ggplot <- NA 
  }
}

for (j in 7:length(lamb)) {
  lam = lamb[j]
  cat("***lambda::",lam,"\n")
  
  
  for (i in 1:length(st)) {
    s_time = Sys.time()
    cat("**지역::",st[i],"\n")
    x = subset(Pr_46,stnlds == st[i])
    # x = subset(Pr_46,stnlds == 273)
    x = x$pr
    
    #init_value
    n = length(x)
    z = diag(n)
    true_beta = rep(mean(x),n)
    
    dmatrix = mat_func(n)
    init_rho = 0.5
    z_init = rep(0,n-1) ; u_init = rep(1,n-1)
    
    A = dmatrix
    B = -diag(n-1)
    mlefit <- fgev(x)
    start <- list()
    #start$scale <- sqrt(6 * var(x))/pi
    #start$loc <- mean(x) - 0.58 * start$scale 
    tvec = c(true_beta,0,0)
    tvec[1:n] = mlefit$estimate[1]
    tvec[n+1] = mlefit$estimate[2]
    tvec[n+2] = mlefit$estimate[3]
    
    #update
    # s_time = Sys.time()
    for (iter in 1:10000) {
      
      #mu, sigma, k optim
      old_tvec = gevreg(x = x,z = z,ctr_list = ctr_list)
      old_mu = old_tvec[1:n]
      
      # z update in ADMM
      tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = init_rho, a = a)
      
      AA = init_rho*t(A)%*%B
      sk = base::norm(drop(AA%*%(tmp_z -z_init)),"2") ## ?? z_init 이 맞나?
      rk = base::norm(drop(A%*%old_mu + B%*%tmp_z),"2")
      
      if (rk<epri & sk <edul) break
      
      # u update in ADMM
      tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
      
      u_init = tmp_u ; z_init = tmp_z
      tvec = old_tvec
      
      # cat("iter::",iter,"\n")
    }
    
    nonzero = length(tmp_z[which(tmp_z != 0)])
    df = nonzero + 3
    logl = dgev(x, loc = old_tvec[1:n], scale = old_tvec[n+1], shape = old_tvec[n+2], log = T)
    loss = -sum(logl)
    
    aic = 2*loss + 2*df
    bic = 2*loss + log(n)*df
    e_time = Sys.time()
    
    results_list[[i]][[j]]$x_mu$x <- x; results_list[[i]][[j]]$x_mu$mu <- tvec[1:n]
    results_list[[i]][[j]]$tmp_z_u$tmp_z <- tmp_z ; results_list[[i]][[j]]$tmp_z_u$tmp_u <- tmp_u
    results_list[[i]][[j]]$theta <- tvec[(n+1):(n+2)]
    results_list[[i]][[j]]$df_aic_bic <- c(df, aic, bic)
    
    gg <- ggplot(data = results_list[[i]][[j]]$x_mu, aes(x = 1973:2018, y = x)) + xlab('time') + ylab('precl') + geom_point(shape = 4)
    gg <- gg + geom_line(aes(x = 1973:2018, y = mu), col = 'red') + ggtitle(paste('st : ',st[i] , ', lambda : ', lamb[j])) +
      theme(plot.title = element_text(hjust = 0.5))
    results_list[[i]][[j]]$ggplot <- gg
    
    print(e_time- s_time)
    cat("primal error : ", rk, " dual error :", sk,'\n')
    cat("iter::",iter,"\n")
  }
}

{
# save(results_list, file = "results_list_1128.RData")
# ***lambda:: 0.05 
# **지역:: 114 
# Time difference of 17.93297 mins
# primal error :  4.552626  dual error : 1.471912 
# iter:: 10000 
# **지역:: 138 
# Time difference of 15.58346 mins
# primal error :  2.434133  dual error : 0.8677269 
# iter:: 10000 
# **지역:: 211 
# Time difference of 18.58405 mins
# primal error :  4.244641  dual error : 0.9718126 
# iter:: 10000 
# **지역:: 247 
# Time difference of 18.69302 mins
# primal error :  3.460455  dual error : 0.971187 
# iter:: 10000 
# **지역:: 285 
# Time difference of 17.54522 mins
# primal error :  3.687347  dual error : 1.170695 
# iter:: 10000 
# ***lambda:: 0.07 
# **지역:: 114 
# Time difference of 15.97342 mins
# primal error :  3.773169  dual error : 1.162166 
# iter:: 10000 
# **지역:: 138 
# Time difference of 4.058727 mins
# primal error :  0.0902459  dual error : 0 
# iter:: 10000 
# **지역:: 211 
# Time difference of 7.911472 mins
# primal error :  3.8874  dual error : 0.9623438 
# iter:: 10000 
# **지역:: 247 
# Time difference of 3.054789 mins
# primal error :  0.05547775  dual error : 0 
# iter:: 10000 
# **지역:: 285 
# Time difference of 8.890117 secs
# primal error :  0.0009997085  dual error : 0 
# iter:: 1162 
# ***lambda:: 0.1 
# **지역:: 114 
# Time difference of 3.209786 mins
# primal error :  0.3962766  dual error : 0 
# iter:: 10000 
# **지역:: 138 
# Time difference of 2.887511 mins
# primal error :  0.1898261  dual error : 0 
# iter:: 10000 
# **지역:: 211 
# Time difference of 2.672975 mins
# primal error :  0.004188779  dual error : 0 
# iter:: 10000 
# **지역:: 247 
# Time difference of 2.81543 mins
# primal error :  0.09025393  dual error : 0 
# iter:: 10000 
# **지역:: 285 
# Time difference of 1.347788 mins
# primal error :  0.0007845162  dual error : 0 
# iter:: 5917 
# ***lambda:: 0.15 
# **지역:: 114 
# Time difference of 2.900345 mins
# primal error :  0.06087825  dual error : 0 
# iter:: 10000 
# **지역:: 138 
# Time difference of 2.718092 mins
# primal error :  0.1109347  dual error : 0 
# iter:: 10000 
# **지역:: 211 
# Time difference of 2.829563 mins
# primal error :  0.4647491  dual error : 0 
# iter:: 10000 
# **지역:: 247 
# Time difference of 2.881819 mins
# primal error :  0.1371383  dual error : 0 
# iter:: 10000 
# **지역:: 285 
# Time difference of 18.14044 secs
# primal error :  0.0009748383  dual error : 0 
# iter:: 2288 
# ***lambda:: 1 
# **지역:: 114 
# Time difference of 2.734305 mins
# primal error :  0.267144  dual error : 0 
# iter:: 10000 
# **지역:: 138 
# Time difference of 2.675074 mins
# primal error :  0.08132487  dual error : 0 
# iter:: 10000 
# **지역:: 211 
# Time difference of 1.25 mins
# primal error :  0.0009994441  dual error : 0 
# iter:: 8219 
# **지역:: 247 
# Time difference of 2.665686 mins
# primal error :  0.02565204  dual error : 0 
# iter:: 10000 
# **지역:: 285 
# Time difference of 27.89967 secs
# primal error :  0.0009980172  dual error : 0 
# iter:: 2975 
}
# scad_result <- results_list
# save(scad_result, file = "scad_result.RData")
