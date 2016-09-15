source("./readArabdopsisData.R")
rstan_options(auto_write = TRUE)
options(mc.cores = 10)

arabi_data_std = arabi_data
mask_L = arabi_data_std$partner == "L"
arabi_data_std[ mask_L, traits] <- scale(arabi_data_std[ mask_L, traits])
arabi_data_std[!mask_L, traits] <- scale(arabi_data_std[!mask_L, traits])

input_data = list(N_s     = sum(!mask_L),
                  N_d     = sum( mask_L),
                  K       = length(traits),
                  n_RIL   = length(RIL_levels),
                  RIL_s   = arabi_data$RIL[!mask_L],
                  RIL_d   = arabi_data$RIL[ mask_L],
                  batch_s = arabi_data$block[!mask_L] - 1,
                  batch_d = arabi_data$block[ mask_L] - 1,
                  y_s     = as.matrix(arabi_data_std[!mask_L, traits]),
                  y_d     = as.matrix(arabi_data_std[ mask_L, traits])
                  )

fit = stan(file = "./VarianceDecomp.stan", data = input_data, iter = 400, chains = 10, 
           control = list(adapt_delta = 0.9))

summary(fit, pars = "G")
G = colMeans(extract(fit, pars = "G")[[1]])
R_s = colMeans(extract(fit, pars = "R_s")[[1]])
R_d = colMeans(extract(fit, pars = "R_d")[[1]])
diag(G)
