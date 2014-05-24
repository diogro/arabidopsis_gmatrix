data{
    int<lower=0> N;
    int<lower=0> branch[N];
    real partnerNONE[N];
    int RIL[N];
    int N_RIL;
}

transformed data{
    vector[2] zeros_RIL;
    for ( i in 1:2 ) zeros_RIL[i] <- 0;
}

parameters{
    real Intercept;
    real beta_partnerNONE;
    vector[2] vary_RIL[N_RIL];
    cov_matrix[2] Sigma_RIL;
}

model{
    real vary[N];
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    beta_partnerNONE ~ normal( 0 , 100 );
    // Varying effects
    for ( j in 1:N_RIL ) vary_RIL[j] ~ multi_normal( zeros_RIL , Sigma_RIL );
    // Fixed effects
    for ( i in 1:N ) {
        vary[i] <- vary_RIL[RIL[i],1]
                + vary_RIL[RIL[i],2] * partnerNONE[i];
        glm[i] <- vary[i] + Intercept
                + beta_partnerNONE * partnerNONE[i];
        glm[i] <- exp( glm[i] );
    }
    branch ~ poisson( glm );
}

generated quantities{
    real dev;
    real vary[N];
    real glm[N];
    dev <- 0;
    for ( i in 1:N ) {
        vary[i] <- vary_RIL[RIL[i],1]
                + vary_RIL[RIL[i],2] * partnerNONE[i];
        glm[i] <- vary[i] + Intercept
                + beta_partnerNONE * partnerNONE[i];
        dev <- dev + (-2) * poisson_log( branch[i] , exp(glm[i]) );
    }
}
 
