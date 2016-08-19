functions {
    // square root of a matrix (elementwise)
    matrix sqrt_mat(matrix x) {
        matrix[dims(x)[1], dims(x)[2]] res;
        for (n in 1:dims(x)[2]){
            for (m in 1:dims(x)[1]){
                res[m, n] = sqrt(x[m, n]);
            }
        }
        return res;
    }
    // square root of a vector (elementwise)
    vector sqrt_vec(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1]){
            res[m] = sqrt(x[m]);
        }
        return res;
    }
}
data {
    int<lower=0> N_s;      // individuals
    int<lower=0> N_d;      // individuals
    int<lower=1> K;      // traits
    int<lower=1> J;      // loci
    int<lower=1> miss_J; // missing loci
    int row_mis[miss_J];
    int col_mis[miss_J];
    int n_RIL;           // inbred lines
    int RIL_s[N_s];
    int RIL_d[N_d];
    vector[J] obs_mk[n_RIL];
    vector[K] y_s[N_s];
    vector[K] y_d[N_d];
}
parameters {
    # Horseshoe prior
    vector<lower=0>[K] r1_global;
    vector<lower=0>[K] r2_global;
    matrix<lower=0>[K, J] r1_local;
    matrix<lower=0>[K, J] r2_local;
    matrix<lower=0>[K, J] r1_localPlus;
    matrix<lower=0>[K, J] r2_localPlus;

    vector<lower=0>[K] r1_global_env;
    vector<lower=0>[K] r2_global_env;
    matrix<lower=0>[K, J] r1_local_env;
    matrix<lower=0>[K, J] r2_local_env;
    matrix<lower=0>[K, J] r1_localPlus_env;
    matrix<lower=0>[K, J] r2_localPlus_env;

    # Effects matrix
    matrix[K,J] beta;
    matrix[K,J] beta_env;

    vector<lower = -1, upper = 1>[miss_J] mis_mk;

    # Intercept
    vector[K] w0;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R_s;
    vector<lower=0>[K] r1_sigma_R_s;
    vector<lower=0>[K] r2_sigma_R_s;

    cholesky_factor_corr[K] L_Omega_R_d;
    vector<lower=0>[K] r1_sigma_R_d;
    vector<lower=0>[K] r2_sigma_R_d;

}
transformed parameters{
    vector[J] mk[n_RIL];

    // global and local variance parameters, and the input weights
    vector<lower=0>[K] L_sigma_R_s;
    vector<lower=0>[K] L_sigma_R_d;
    vector<lower=0>[K] tau;
    matrix<lower=0>[K, J] lambda;
    matrix<lower=0>[K, J] eta;
    matrix<lower=0>[K, J] etaLambda;
    matrix<lower=0>[K, J] sd_theta;
    matrix[K, J] w;
    vector<lower=0>[K] tau_env;
    matrix<lower=0>[K, J] lambda_env;
    matrix<lower=0>[K, J] eta_env;
    matrix<lower=0>[K, J] etaLambda_env;
    matrix<lower=0>[K, J] sd_theta_env;
    matrix[K, J] w_env;
    matrix[K,K] L_Sigma_R_s;
    matrix[K,K] L_Sigma_R_d;

    L_sigma_R_s = 2.5 * r1_sigma_R_s .* sqrt_vec(r2_sigma_R_s);
    L_Sigma_R_s = diag_pre_multiply(L_sigma_R_s, L_Omega_R_s);

    L_sigma_R_d = 2.5 * r1_sigma_R_d .* sqrt_vec(r2_sigma_R_d);
    L_Sigma_R_d = diag_pre_multiply(L_sigma_R_d, L_Omega_R_d);

    tau = r1_global .* sqrt_vec(r2_global) .* L_sigma_R_s;

    lambda = r1_local .* sqrt_mat(r2_local);

    eta = r1_localPlus .* sqrt_mat(r2_localPlus);

    etaLambda = lambda .* eta;
    for(j in 1:J){
        for(k in 1:K){
            sd_theta[k, j] = etaLambda[k, j] * tau[k];
        }
    }
    w = beta .* sd_theta;

    tau_env = r1_global_env .* sqrt_vec(r2_global_env) .* (L_sigma_R_d);

    lambda_env = r1_local_env .* sqrt_mat(r2_local_env);

    eta_env = r1_localPlus_env .* sqrt_mat(r2_localPlus_env);

    etaLambda_env = lambda_env .* eta_env;
    for(j in 1:J){
        for(k in 1:K){
            sd_theta_env[k, j] = etaLambda_env[k, j] * tau_env[k];
        }
    }
    w_env = beta_env .* sd_theta_env;

    for(j in 1:J){
        for(r in 1:n_RIL){
            if(obs_mk[r, j] != 2){
                mk[r, j] = obs_mk[r, j];
            }
        }
    }
    for(mj in 1:miss_J){
        mk[row_mis[mj], col_mis[mj]] = mis_mk[mj];
    }
}
model {
    vector[K] mu_s[N_s];
    vector[K] mu_d[N_d];

    for (n in 1:N_s)
        mu_s[n] = w0 + w * mk[RIL_s[n]];

    for (n in 1:N_d)
        mu_d[n] = w0 + w_env * mk[RIL_d[n]];

    y_s ~ multi_normal_cholesky(mu_s, L_Sigma_R_s);
    y_d ~ multi_normal_cholesky(mu_d, L_Sigma_R_d);

    w0 ~ normal(0,5);

    mis_mk ~ uniform(-1, 1);

    #// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
    to_vector(beta) ~ normal(0, 1);

    to_vector(r1_local) ~ normal(0.0, 1.0);
    to_vector(r2_local) ~ inv_gamma(0.5*3, 0.5*3);
    to_vector(r1_localPlus) ~ normal(0.0, 1.0);
    to_vector(r2_localPlus) ~ inv_gamma(0.5*3, 0.5*3);

    // half cauchy for tau
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    #// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
    to_vector(beta_env) ~ normal(0, 1);

    to_vector(r1_local_env) ~ normal(0.0, 1.0);
    to_vector(r2_local_env) ~ inv_gamma(0.5*3, 0.5*3);
    to_vector(r1_localPlus_env) ~ normal(0.0, 1.0);
    to_vector(r2_localPlus_env) ~ inv_gamma(0.5*3, 0.5*3);

    // half cauchy for tau
    r1_global_env ~ normal(0.0, 1.0);
    r2_global_env ~ inv_gamma(0.5, 0.5);

    L_Omega_R_s ~ lkj_corr_cholesky(2);
    r1_sigma_R_s ~ normal(0.0, 1.0);
    r2_sigma_R_s ~ inv_gamma(0.5, 0.5);
    L_Omega_R_d ~ lkj_corr_cholesky(2);
    r1_sigma_R_d ~ normal(0.0, 1.0);
    r2_sigma_R_d ~ inv_gamma(0.5, 0.5);
}
generated quantities {
    matrix[K, K] R_s;
    matrix[K, K] R_d;
    matrix[K, K] corrR_s;
    matrix[K, K] corrR_d;
    matrix[K, J] shrink;
    matrix[K, J] shrink_env;

    R_s = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R_s, L_Omega_R_s));
    R_d = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R_d, L_Omega_R_d));

    corrR_s = multiply_lower_tri_self_transpose(L_Omega_R_s);
    corrR_d = multiply_lower_tri_self_transpose(L_Omega_R_d);

    for(j in 1:J){
        for(k in 1:K){
            shrink[k, j] = 1 - 1/(1 + (etaLambda[k, j]^2));
        }
    }
    for(j in 1:J){
        for(k in 1:K){
            shrink_env[k, j] = 1 - 1/(1 + (etaLambda_env[k, j]^2));
        }
    }
}

