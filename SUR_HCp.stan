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
    int<lower=0> N;      // individuals
    int<lower=1> K;      // traits
    int<lower=1> J;      // loci
    int<lower=1> miss_J; // missing loci
    int row_mis[miss_J];
    int col_mis[miss_J];
    int n_RIL;           // inbred lines
    int RIL[N];
    vector[J] obs_mk[n_RIL];
    vector[K] y[N];
}
transformed data{
    vector[K] zeros;
    zeros = rep_vector(0.0, K);
}
parameters {
    # Horseshoe prior
    vector<lower=0>[K] r1_global;

    vector<lower=0>[K] r2_global;

    matrix<lower=0>[K, J] r1_local;
    matrix<lower=0>[K, J] r2_local;


    matrix<lower=0>[K, J] r1_localPlus;
    matrix<lower=0>[K, J] r2_localPlus;

    # Effects matrix
    matrix[K,J] beta;

    vector<lower = -1, upper = 1>[miss_J] mis_mk;

    # Intercept
    vector[K] w0;

    # RIL means
    vector[K] beta_RIL[n_RIL];

    # G matrix
    cholesky_factor_corr[K] L_Omega_G;
    vector<lower=0>[K] L_sigma_G;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] r1_sigma_R;
    vector<lower=0>[K] r2_sigma_R;
}
transformed parameters{
    vector[J] mk[n_RIL];

    // global and local variance parameters, and the input weights
    vector<lower=0>[K] L_sigma_R;
    vector<lower=0>[K] tau;
    matrix<lower=0>[K, J] lambda;
    matrix<lower=0>[K, J] eta;
    matrix<lower=0>[K, J] etaLambda;
    matrix<lower=0>[K, J] sd_theta;
    matrix[K, J] w;
    matrix[K,K] L_Sigma_R;

    L_sigma_R = 2.5 * r1_sigma_R .* sqrt_vec(r2_sigma_R);
    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    tau = r1_global .* sqrt_vec(r2_global) .* L_sigma_R;

    lambda = r1_local .* sqrt_mat(r2_local);

    eta = r1_localPlus .* sqrt_mat(r2_localPlus);

    etaLambda = lambda .* eta;
    for(j in 1:J){
        for(k in 1:K){
            sd_theta[k, j] = etaLambda[k, j] * tau[k];
        }
    }
    w = beta .* sd_theta;

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
    vector[K] mu[N];
    matrix[K,K] L_Sigma_G;

    L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);

    for (j in 1:n_RIL)
        beta_RIL[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    for (n in 1:N)
        mu[n] = w0 + w * mk[RIL[n]] + beta_RIL[RIL[n]];

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    w0 ~ normal(0,5);

    L_Omega_G ~ lkj_corr_cholesky(2);
    L_sigma_G ~ cauchy(0, 2.5);

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

    L_Omega_R ~ lkj_corr_cholesky(2);
    r1_sigma_R ~ normal(0.0, 1.0);
    r2_sigma_R ~ inv_gamma(0.5, 0.5);
}
generated quantities {
    matrix[K, K] G;
    matrix[K, K] R;
    corr_matrix[K] corrG;
    corr_matrix[K] corrR;
    matrix[K, J] shrink;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);

    for(j in 1:J){
        for(k in 1:K){
            shrink[k, j] = 1 - 1/(1 + (etaLambda[k, j]^2));
        }
    }
}

