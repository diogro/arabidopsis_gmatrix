data{
    int N;
    real silique[N];
    real partnerNONE[N];
    int RIL[N];
    int bloob[N];
    int N_RIL;
    int N_block;
}

transformed data{
    vector[2] zeros_RIL;
    for ( i in 1:2 ) zeros_RIL[i] <- 0;
}

parameters{
    real Intercept;
    real<lower=0> sigma;
    vector[2] vary_RIL[N_RIL];
    cov_matrix[2] Sigma_RIL;
    real vary_block[N_block];
    real<lower=0> sigma_block;
}

model{
    real vary[N];
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    sigma_block ~ uniform( 0 , 100 );
    sigma ~ uniform( 0 , 100 );
    // Varying effects
    for ( j in 1:N_RIL ) vary_RIL[j] ~ multi_normal( zeros_RIL , Sigma_RIL );
    for ( j in 1:N_block ) vary_block[j] ~ normal( 0 , sigma_block );
    // Fixed effects
    for ( i in 1:N ) {
        vary[i] <- vary_RIL[RIL[i],1]
                + vary_RIL[RIL[i],2] * partnerNONE[i]
                + vary_block[bloob[i]];
        glm[i] <- vary[i] + Intercept;
    }
    silique ~ normal( glm , sigma );
}


