rm(list = ls())

if(!require(tidyverse)){install.packages("tidyverse")};require(tidyverse)
if(!require(ggplot2)){install.packages("ggplot2")};require(ggplot2)


load("Pr_46.RData")
source("3_update_gev_function.R")

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

for (j in 1:length(lamb)) {
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
      tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = init_rho)
      
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
    cat("iter::",iter,"\n")
  }
}

# save(results_list, file = "results_list_1128.RData")


