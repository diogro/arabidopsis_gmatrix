library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)

raw_arabi_data <- read.csv2("./raw_data.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT3)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

arabi_data$flower[is.na(arabi_data$flower)] <- 0
arabi_data = arabi_data[complete.cases(arabi_data),]
#arabi_data = arabi_data[arabi_data$flower > 0,]

arabi_data$weight  <- sqrt(arabi_data$weight)
arabi_data$silique  <- sqrt(arabi_data$silique)
arabi_data$iszero_silique = arabi_data$silique == 0
#arabi_data$silique[arabi_data$silique == 0] = NA
plot(silique~weight, arabi_data)
plot(silique~height, arabi_data)
table(arabi_data$flower)

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 5, scale = "free")

lme4_model = glmer(branch ~ partner + (partner|RIL),
                   family ="poisson", data = arabi_data)
summary(lme4_model)

traits = c('weight', 'height', 'silique', 'branch')
names_g = paste0(traits, rep(c('D', 'S'), each = 4))

branch_data <- list(N = dim(arabi_data)[1],
                    branch = arabi_data$branch,
                    partnerNONE = as.integer(as.factor(arabi_data$partner)),
                    RIL = as.integer(as.factor(arabi_data$RIL)),
                    N_RIL = length(unique(arabi_data$RIL)))
stan_model = stan(file = './branch.stan', data = branch_data, chain=1)
stan_samples = extract(stan_model, permuted = TRUE)


the_formula <- list(iszero_silique ~ partner + (partner|RIL) ,
                    silique        ~ partner + (partner|RIL))
silique_model = glmer2stan(the_formula, data=arabi_data,
                           family=list("binomial","gamma"),
                           sample = FALSE, calcDIC = FALSE)
write(silique_model$model, file = "silique.stan")


N_pi              = dim(arabi_data)[1]
N_mu              = sum(arabi_data$silique!=0)
iszero_silique_pi = as.integer(arabi_data$iszero_silique)
partnerNONE_pi    = as.integer(as.factor(arabi_data$partner))
RIL_pi            = as.integer(as.factor(arabi_data$RIL))
bin_total_pi      = rep(1,dim(arabi_data)[1])
silique_mu        = arabi_data$silique[!arabi_data$iszero_silique]
partnerNONE_mu    = as.integer(as.factor(arabi_data$partner))[!arabi_data$iszero_silique]
RIL_mu            = as.integer(as.factor(arabi_data$RIL))[!arabi_data$iszero_silique]
N_RIL             = length(unique(arabi_data$RIL))

silique_data <- list(N_pi              = dim(arabi_data)[1],
                     N_mu              = sum(arabi_data$silique!=0),
                     iszero_silique_pi = as.integer(arabi_data$iszero_silique),
                     partnerNONE_pi    = as.integer(as.factor(arabi_data$partner)),
                     RIL_pi            = as.integer(as.factor(arabi_data$RIL)),
                     bin_total_pi      = rep(1,dim(arabi_data)[1]),
                     silique_mu        = arabi_data$silique[!arabi_data$iszero_silique],
                     partnerNONE_mu    = as.integer(as.factor(arabi_data$partner))[!arabi_data$iszero_silique],
                     RIL_mu            = as.integer(as.factor(arabi_data$RIL))[!arabi_data$iszero_silique],
                     N_RIL             = length(unique(arabi_data$RIL)))
stan_model = stan(file = './silique.stan', data = silique_data, chain=1)
stan_samples = extract(stan_model, permuted = TRUE)
attach(stan_samples)
replicates = dim(vary_RIL)[1]
bin_sim = array(0, c(replicates, N_pi))
for(i in 1:replicates){
    vary_pi <- vary_RIL[i, RIL_pi,1] + vary_RIL[i, RIL_pi,2] * partnerNONE_pi
    glm_pi <- vary_pi + Intercept_pi[i] + beta_partnerNONE_pi[i] * partnerNONE_pi
    glm_pi <- inv.logit( glm_pi )
    for(j in 1:N_pi)
        bin_sim[i,j] = rbinom(1, 1, prob = glm_pi[j])
}
hist(apply(bin_sim, 1, sum)); abline(v = sum(arabi_data$iszero_silique), col = "red")
post_zeros = mean(apply(bin_sim, 1, sum))

silique_sim = array(0, c(replicates, N_mu))
for(i in 1:replicates){
    vary_mu <- vary_RIL[i, RIL_mu, 3] + vary_RIL[i, RIL_mu,4] * partnerNONE_mu
    glm_mu <- vary_mu + Intercept_mu[i] + beta_partnerNONE_mu[i] * partnerNONE_mu
    glm_mu <- exp( glm_mu )*theta_mu[i];
    for ( j in 1:N_mu ) {
        silique_sim[i,j] = rgamma(1, glm_mu[j] , theta_mu[i])
    }
}
par(mfrow=c(1, 2))
hist((apply(silique_sim, 2, mean))^2, xlim=c(0, 150))
hist(arabi_data$silique[!arabi_data$iszero_silique]^2, xlim=c(0, 150))
