source("./readArabdopsisData.R")
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

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

fit = stan(file = "./VarianceDecomp.stan", data = input_data, iter = 200, chains = 3)
plot(fit, pars= "shrink")
plot(fit, pars= "shrink_env")
plot(fit, pars= "w")
plot(fit, pars= "w_env")
plot(fit, pars= c("w0", "wb_s", "wb_d", "wd"))

colMeans(rstan::extract(fit, pars = "shrink")[[1]]) > 0.4
colMeans(rstan::extract(fit, pars = "w")[[1]])

x = colMeans(rstan::extract(fit, pars = "shrink")[[1]]) > 0.4
x_env =  colMeans(rstan::extract(fit, pars = "shrink_env")[[1]]) > 0.4
x[] = as.numeric(x)
x_env[] = as.numeric(x_env)
plot(x * colMeans(rstan::extract(fit, pars = "w")[[1]]),
     x_env *  colMeans(rstan::extract(fit, pars = "w_env")[[1]]))
abline(0, 1)

# BLUPS

library(lme4)
library(evolqg)

fit_weight_L = lmer(weight ~ block + (1|RIL), data = arabi_data[mask_L,])
fit_weight_NONE = lmer(weight ~ block + (1|RIL), data = arabi_data[!mask_L,])
BLUPs_wht_L = ranef(fit_weight_L)[[1]]
BLUPs_wht_NONE = ranef(fit_weight_NONE)[[1]]
b_wht = bind_rows(data_frame(RIL = rownames(BLUPs_wht_L), weight = as.vector(scale(BLUPs_wht_L[[1]])), partner = "L"),
                  data_frame(RIL = rownames(BLUPs_wht_NONE), weight = as.vector(scale(BLUPs_wht_NONE[[1]])), partner = "NONE"))

fit_height_L = lmer(height ~ block + (1|RIL), data = arabi_data[mask_L,])
fit_height_NONE = lmer(height ~ block + (1|RIL), data = arabi_data[!mask_L,])
BLUPs_ht_L = ranef(fit_height_L)[[1]]
BLUPs_ht_NONE = ranef(fit_height_NONE)[[1]]
b_ht = bind_rows(data_frame(RIL = rownames(BLUPs_ht_L), height = as.vector(scale(BLUPs_ht_L[[1]])), partner = "L"),
                 data_frame(RIL = rownames(BLUPs_ht_NONE), height = as.vector(scale(BLUPs_ht_NONE[[1]])), partner = "NONE"))

fit_silique_L = lmer(silique ~ block + (1|RIL), data = arabi_data[mask_L,])
fit_silique_NONE = lmer(silique ~ block + (1|RIL), data = arabi_data[!mask_L,])
BLUPs_slq_L = ranef(fit_silique_L)[[1]]
BLUPs_slq_NONE = ranef(fit_silique_NONE)[[1]]
b_slq = bind_rows(data_frame(RIL = rownames(BLUPs_slq_L), silique = as.vector(scale(BLUPs_slq_L[[1]])), partner = "L"),
                  data_frame(RIL = rownames(BLUPs_slq_NONE), silique = as.vector(scale(BLUPs_slq_NONE[[1]])), partner = "NONE"))

blup_arabi_data = inner_join(b_slq, inner_join(b_wht, b_ht, by = c("RIL", "partner")), by = c("RIL", "partner"))

blup_input = list(N_s = nrow(blup_arabi_data %>% filter(partner == "NONE")),
                  N_d = nrow(blup_arabi_data %>% filter(partner == "L")),
                  J       = dim(markers)[2],
                  K       = length(traits),
                  miss_J  = miss_J,
                  row_mis = mis_pos[,1],
                  col_mis = mis_pos[,2],
                  n_RIL   = length(RIL_levels),
                  RIL_s   = as.numeric(filter(blup_arabi_data, partner == "NONE")$RIL),
                  RIL_d   = as.numeric(filter(blup_arabi_data, partner == "L")$RIL),
                  obs_mk  = obs_markers,
                  y_s     = as.matrix(filter(blup_arabi_data, partner == "NONE")[traits]),
                  y_d     = as.matrix(filter(blup_arabi_data, partner == "L")[traits]))

fit_blups = stan(file = "./SUR_HCp_BLUPS.stan", data = blup_input, iter = 200, chains = 3)

plot(fit_blups, pars= "shrink")
plot(fit_blups, pars= "shrink_env")
plot(fit_blups, pars= "w")
plot(fit_blups, pars= "w_env")


cut_off = 0.45
sig_s      = colMeans(rstan::extract(fit, pars       = "shrink")[[1]]) > cut_off
sig_d      = colMeans(rstan::extract(fit, pars       = "shrink_env")[[1]]) > cut_off
blup_sig_s = colMeans(rstan::extract(fit_blups, pars = "shrink")[[1]]) > cut_off
blup_sig_d = colMeans(rstan::extract(fit_blups, pars = "shrink_env")[[1]]) > cut_off
w_s        = colMeans(rstan::extract(fit, pars       = "w")[[1]])
w_d        = colMeans(rstan::extract(fit, pars       = "w_env")[[1]])
blup_w_s   = colMeans(rstan::extract(fit_blups, pars = "w")[[1]])
blup_w_d   = colMeans(rstan::extract(fit_blups, pars = "w_env")[[1]])
res = data.frame(S = apply(sig_s, 2, any),
                 blup_S = apply(blup_sig_s, 2, any),
                 D = apply(sig_d, 2, any),
                 blup_D = apply(blup_sig_d, 2, any),
                 w_S = apply(w_s, 2, Norm),
                 blup_w_S = apply(blup_w_s, 2, Norm),
                 w_D = apply(w_d, 2, Norm),
                 blup_w_D = apply(blup_w_d, 2, Norm))
plot(w_S ~ blup_w_S, res, pch = as.numeric(blup_sig_s)+1); abline(0, 1)
plot(w_D ~ blup_w_D, res); abline(0, 1)
plot(blup_w_S ~ blup_w_D, res); abline(0, 2)

colSums(res)
