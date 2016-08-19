source("./readArabdopsisData.R")
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

obs_markers = markers
miss_J = sum(is.na(obs_markers))
mis_pos = which(is.na(obs_markers), TRUE)
obs_markers[is.na(obs_markers)] = 2

arabi_data_std = arabi_data
mask_L = arabi_data_std$partner == "L"
arabi_data_std[ mask_L, traits] <- scale(arabi_data_std[ mask_L, traits])
arabi_data_std[!mask_L, traits] <- scale(arabi_data_std[!mask_L, traits])

    int<lower=0> N_s;      // individuals
    int<lower=0> N_d;      // individuals
    int<lower=1> K;      // traits
    int<lower=1> J;      // loci
    int<lower=1> miss_J; // missing loci
    int row_mis[miss_J];
    int col_mis[miss_J];
    int n_RIL;           // inbred lines
    int RIL_s[N_s];
    int batch_s[N_s];        // batch
    int RIL_d[N_d];
    int batch_d[N_d];        // batch
    vector[J] obs_mk[n_RIL];
    vector[K] y_s[N_s];
    vector[K] y_d[N_d];

    input_data = list(N_s     = sum(!mask_L),
                      N_d     = sum( mask_L),
                      J       = dim(markers)[2],
                      K       = length(traits),
                      miss_J  = miss_J,
                      row_mis = mis_pos[,1],
                      col_mis = mis_pos[,2],
                      n_RIL   = length(RIL_levels),
                      RIL_s   = arabi_data$RIL[!mask_L],
                      RIL_d   = arabi_data$RIL[ mask_L],
                      batch_s = arabi_data$block[!mask_L] - 1,
                      batch_d = arabi_data$block[ mask_L] - 1,
                      obs_mk  = obs_markers,
                      y_s     = as.matrix(arabi_data_std[!mask_L, traits]),
                      y_d     = as.matrix(arabi_data_std[ mask_L, traits])
                      )

fit = stan(file = "./SUR_HCp_byEnv.stan", data = input_data, iter = 100, chains = 1)
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
