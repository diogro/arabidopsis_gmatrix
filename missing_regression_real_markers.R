if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "dplyr", "tidyr", "readr", "lme4", "ggplot2", "cowplot", "MCMCglmm", "rstan", "mvtnorm")
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

raw_markers = read_csv("./data/BayxSha_2_Genotypes.csv")[-c(1:4),] %>% rename(RIL = order) %>% mutate(RIL = as.integer(RIL))
raw_markers[] = lapply(raw_markers, function(x) {
                           x[x == "A"] = 1;
                           x[x == "B"] = -1;
                           x[x == "C"] = 0;
                           x[x == "D"] = NA;
                       return(as.integer(x)) } )
names(raw_markers)[2:70] = paste0("G", 1:69)
markers = semi_join(raw_markers, arabi_data, id = "RIL")
arabi_data = semi_join(arabi_data, raw_markers, id = "RIL")
RIL_levels = unique(markers$RIL)
arabi_data$RIL = as.integer(factor(arabi_data$RIL, RIL_levels))
markers = as.matrix(select(markers, G1:G69))

obs_markers = markers
miss_J = sum(is.na(obs_markers))
mis_pos = which(is.na(obs_markers), TRUE)
obs_markers[is.na(obs_markers)] = 2
#obs_markers[is.na(obs_markers)] = 2

markers[is.na(markers)] = 1

sparseEffect = function(x, s = 2, p = c(0.1, 0.9)) ifelse(sample(c(TRUE, FALSE), x, T, p), rnorm(x, 1, s), 0)
N = dim(arabi_data)[1]
J = dim(markers)[2]
n_RIL = length(RIL_levels)
RIL   = arabi_data$RIL
add = matrix(sparseEffect(J*2), J, 2)

bRIL = matrix(rmvnorm(n_RIL, c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)), n_RIL, 2)
y = (markers)[RIL,] %*% add + bRIL[RIL,] + rnorm(N*2, 0, 1)

input_data = list(N = N, K = 2, J = J, miss_J = miss_J, row_mis = mis_pos[,1], col_mis = mis_pos[,2], n_RIL = n_RIL, RIL = RIL, y = y, obs_mk = obs_markers)
fit = stan(file = "./SUR_HCp.stan", data = input_data, iter = 200, chains = 3)
plot(fit, pars= "mis_mk")
plot(fit, pars= "beta_RIL")
plot(fit, pars= "shrink")
plot(fit, pars= "w")
plot(add, t(colMeans(rstan::extract(fit)$w)))
plot(bRIL, colMeans(rstan::extract(fit)$beta_RIL))
traceplot(fit, pars = "beta_RIL[1,1]", inc_warmup = T)
abline(0, 1)
