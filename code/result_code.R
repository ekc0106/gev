## 결과
st = c(114,138,211,247,285)
lamb <- c(0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15 , 1)

load(file = 'C:/Users/kyucheol/Desktop/kma_data/results_list_1128.RData')
load(file = 'C:/Users/kyucheol/Desktop/kma_data/result_refit.RData')

abic_mat_fun <- function(result_list, aic_bic){
  # aic : 2, bic : 3
  mat <- matrix(0,length(st),length(lamb))
  
  for(i in 1:length(st)){
    for(j in 1:length(lamb)){
      mat[i,j] <- result_list[[i]][[j]]$df_aic_bic[aic_bic]
    }
  }
  colnames(mat) <- names(results_list[[1]]) ; rownames(mat) <- names(results_list)
  return(mat)
}
lam_min_idx_fun <- function(mat){
  lam_min_idx <- apply(mat, 1, which.min) %>% as.numeric
  return(lam_min_idx)
}
cp_fun <- function(result_list){
  change_point_mat1 <- matrix(0,length(st),length(lamb))
  
  
  for(i in 1:length(st)){
    for(j in 1:length(lamb)){
      change_point_mat1[i,j] <- result_list[[i]][[j]]$df_aic_bic[1]-3
    }
  }
  colnames(change_point_mat1) <- names(result_list[[1]]) ; rownames(change_point_mat1) <- names(result_list)
  return(change_point_mat1)
}
results_list_fin_fun <- function(result_list, lam_min_idx){
  results_list_fin <- vector(mode = 'list', length = length(st))
  names(results_list_fin) <- paste('st',st,sep = '_')
  
  for(idx in 1:length(results_list_fin)){
    results_list_fin[[idx]] <- result_list[[idx]][[lam_min_idx[idx]]]
    results_list_fin[[idx]]$lambda <- lamb[lam_min_idx[idx]]
  }
  return(results_list_fin)
}

## aic , bic
non_refit_aic <- abic_mat_fun(result_list = results_list, 2)
non_refit_bic <- abic_mat_fun(result_list = results_list, 3)
refit_aic <- abic_mat_fun(result_list = result_refit, 2)
refit_bic <- abic_mat_fun(result_list = result_refit, 3)

# lamb_min_idx
non_refit_lamb_a <- lam_min_idx_fun(non_refit_aic)
non_refit_lamb_b <- lam_min_idx_fun(non_refit_bic)
refit_lamb_a <- lam_min_idx_fun(refit_aic)
refit_lamb_b <- lam_min_idx_fun(refit_bic)

# change point
non_refit_cp <- cp_fun(results_list)
refit_cp <- cp_fun(result_refit)

# final_result

non_refit_final_aic <- results_list_fin_fun(result_list = results_list, lam_min_idx = non_refit_lamb_a)
non_refit_final_bic <- results_list_fin_fun(result_list = results_list, lam_min_idx = non_refit_lamb_b)
refit_final_aic <- results_list_fin_fun(result_list = result_refit, lam_min_idx = refit_lamb_a)
refit_final_bic <- results_list_fin_fun(result_list = result_refit, lam_min_idx = refit_lamb_b)

save(non_refit_final_aic, file = "non_refit_final_aic.RData")
save(non_refit_final_bic, file = "non_refit_final_bic.RData")
save(refit_final_aic, file = "refit_final_aic.RData")
save(refit_final_bic, file = "refit_final_bic.RData")
