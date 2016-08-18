library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

markers = matrix(sample(c(0, 1), 100*10, replace = TRUE), 100, 10)
add = matrix(rnorm(10*2, 0, 1), 10, 2)
add[1,1] = 10

y = markers %*% add + rnorm(100*2, 0, 1)

obs_markers = markers
markers[3, 6]
obs_markers[1, 1] = 2
obs_markers[3, 6] = 2

model = "
data {
    int<lower=1> K;      // traits
    int<lower=1> J;      // loci
    int<lower=1> miss_J; // missing loci
    int<lower=0> N;      // individuals
    vector[J] obs_mk[N];
    vector[K] y[N];
}
parameters {
    # Effects matrix
    matrix[K,J] beta;

    vector<lower = 0, upper = 1>[miss_J] mis_mk;

    # Intercept
    vector[K] w0;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] L_sigma_R;
}
transformed parameters{
    vector[J] mk[N];

    for(j in 1:J){
        for(n in 1:N){
            if(obs_mk[n, j] != 2){
                mk[n, j] = obs_mk[n, j];
            }
        }
    }
    mk[1, 1] = mis_mk[1];
    mk[3, 6] = mis_mk[2];
}
model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_R;

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for (n in 1:N)
        mu[n] = w0 + beta * mk[n];

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    w0 ~ normal(0,5);

    L_Omega_R ~ lkj_corr_cholesky(2);
    L_sigma_R ~ cauchy(0, 2.5);

    mis_mk ~ uniform(0, 1);
}

"

input_data = list(K = 2, J = 10, miss_J = 2, N = 100, y = y, obs_mk = obs_markers)
fit = stan(model_code = model, data = input_data, iter = 2000, chains = 1)
plot(fit, pars= "mis_mk")
plot(add, t(colMeans(rstan::extract(fit)$beta)))
abline(0, 1)
