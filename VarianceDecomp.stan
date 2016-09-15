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
    int n_RIL;           // inbred lines
    int RIL_s[N_s];
    int batch_s[N_s];        // batch
    int RIL_d[N_d];
    int batch_d[N_d];        // batch
    vector[K] y_s[N_s];
    vector[K] y_d[N_d];
}
transformed data{
    vector[2*K] zeros;
    zeros = rep_vector(0.0, 2*K);
}
parameters {
    # Intercept
    vector[K] w0;
    vector[K] wb_s;
    vector[K] wb_d;
    vector[K] wd;

    # RIL means
    vector[2*K] beta_RIL[n_RIL];

    # G matrix
    cholesky_factor_corr[2*K] L_Omega_G;
    vector<lower=0>[2*K] L_sigma_G;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R_s;
    vector<lower=0>[K] r1_sigma_R_s;
    vector<lower=0>[K] r2_sigma_R_s;

    cholesky_factor_corr[K] L_Omega_R_d;
    vector<lower=0>[K] r1_sigma_R_d;
    vector<lower=0>[K] r2_sigma_R_d;

}
transformed parameters{
    // global and local variance parameters, and the input weights
    vector<lower=0>[K] L_sigma_R_s;
    vector<lower=0>[K] L_sigma_R_d;
    matrix[K,K] L_Sigma_R_s;
    matrix[K,K] L_Sigma_R_d;

    L_sigma_R_s = 2.5 * r1_sigma_R_s .* sqrt_vec(r2_sigma_R_s);
    L_Sigma_R_s = diag_pre_multiply(L_sigma_R_s, L_Omega_R_s);

    L_sigma_R_d = 2.5 * r1_sigma_R_d .* sqrt_vec(r2_sigma_R_d);
    L_Sigma_R_d = diag_pre_multiply(L_sigma_R_d, L_Omega_R_d);

}
model {
    vector[K] mu_s[N_s];
    vector[K] mu_d[N_d];
    matrix[2*K,2*K] L_Sigma_G;

    L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);

    for (j in 1:n_RIL)
        beta_RIL[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    for (n in 1:N_s)
        mu_s[n] = w0 + batch_s[n] * wb_s +
        beta_RIL[RIL_s[n], 1:K];

    for (n in 1:N_d)
        mu_d[n] = w0 + batch_d[n] * wb_d + wd +
        beta_RIL[RIL_d[n], (K+1):(2*K)];

    y_s ~ multi_normal_cholesky(mu_s, L_Sigma_R_s);
    y_d ~ multi_normal_cholesky(mu_d, L_Sigma_R_d);

    w0 ~ normal(0,5);
    wd ~ normal(0,5);
    wb_s ~ normal(0,5);
    wb_d ~ normal(0,5);

    L_Omega_G ~ lkj_corr_cholesky(2);
    L_sigma_G ~ cauchy(0, 2.5);

    L_Omega_R_s ~ lkj_corr_cholesky(2);
    r1_sigma_R_s ~ normal(0.0, 1.0);
    r2_sigma_R_s ~ inv_gamma(0.5, 0.5);
    L_Omega_R_d ~ lkj_corr_cholesky(2);
    r1_sigma_R_d ~ normal(0.0, 1.0);
    r2_sigma_R_d ~ inv_gamma(0.5, 0.5);
}
generated quantities {
    matrix[2*K, 2*K] G;
    matrix[K, K] R_s;
    matrix[K, K] R_d;
    corr_matrix[2*K] corrG;
    corr_matrix[K] corrR_s;
    corr_matrix[K] corrR_d;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R_s = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R_s, L_Omega_R_s));
    R_d = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R_d, L_Omega_R_d));

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR_s = multiply_lower_tri_self_transpose(L_Omega_R_s);
    corrR_d = multiply_lower_tri_self_transpose(L_Omega_R_d);
}

