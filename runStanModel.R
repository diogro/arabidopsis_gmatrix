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

arabi_data$weight  <- log(arabi_data$weight)
mask_0 = arabi_data$iszero_silique = arabi_data$silique == 0
arabi_data$silique[!mask_0]  <- log(arabi_data$silique[!mask_0])
#arabi_data$silique[arabi_data$silique == 0] = NA
plot(silique~weight, arabi_data)
plot(silique~height, arabi_data)
table(arabi_data$flower)

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 5, scale = "free")

lme4_model = glmer(branch ~ partner + (partner|RIL) + (1|block),
                   family ="poisson", data = arabi_data)
summary(lme4_model)

traits = c('weight', 'height', 'silique', 'branch')
names_g = paste0(traits, rep(c('D', 'S'), each = 4))

branch_model = glmer2stan(branch ~ partner + (partner|RIL) + (1|block),
                          data=arabi_data,
                          family="poisson"
                          sample = FALSE, calcDIC = FALSE)
write(silique_model$model, file = "branch.stan")

branch_data <- list(N = dim(arabi_data)[1],
                    branch = arabi_data$branch,
                    partnerNONE = as.integer(as.factor(arabi_data$partner)),
                    RIL = as.integer(as.factor(arabi_data$RIL)),
                    N_RIL = length(unique(arabi_data$RIL)))
branch_stan_model = stan(file = './branch.stan', data = branch_data, chain=1)
bm = extract(branch_stan_model, permuted = TRUE)


the_formula <- list(iszero_silique ~ partner + (partner|RIL) + block,
                    silique        ~ partner + (partner|RIL) + block)
silique_model = glmer2stan(the_formula, data=arabi_data,
                           family=list("binomial","gaussian"),
                           sample = FALSE, calcDIC = FALSE)
write(silique_model$model, file = "silique.stan")

N_pi              = dim(arabi_data)[1]
N_mu              = sum(!mask_0)
iszero_silique_pi = as.integer(mask_0)
partnerNONE_pi    = as.integer(as.factor(arabi_data$partner)) - 1
RIL_pi            = as.integer(as.factor(arabi_data$RIL))
block_pi          = as.integer(as.factor(arabi_data$block)) - 1
bin_total_pi      = rep(1,dim(arabi_data)[1])
silique_mu        = arabi_data$silique[!mask_0]
partnerNONE_mu    = partnerNONE_pi[!mask_0]
RIL_mu            = as.integer(as.factor(arabi_data$RIL))[!mask_0]
block_mu          = block_pi[!mask_0]
N_RIL             = length(unique(arabi_data$RIL))
silique_data <- list(N_pi              = N_pi,
                     N_mu              = N_mu,
                     iszero_silique_pi = iszero_silique_pi,
                     partnerNONE_pi    = partnerNONE_pi,
                     RIL_pi            = RIL_pi,
                     block_pi          = block_pi,
                     bin_total_pi      = bin_total_pi,
                     silique_mu        = silique_mu,
                     partnerNONE_mu    = partnerNONE_mu,
                     RIL_mu            = RIL_mu,
                     block_mu          = block_mu,
                     N_RIL             = N_RIL)
silique_stan_model = stan(file = './silique.stan', data = silique_data, chain=1)

sm = extract(silique_stan_model, permuted = TRUE)
replicates = dim(sm$vary_RIL)[1]
bin_sim = array(0, c(replicates, N_pi))
for(i in 1:replicates){
    vary_pi <- sm$vary_RIL[i, RIL_pi,1] +
               sm$vary_RIL[i, RIL_pi,2] * partnerNONE_pi
    glm_pi <- vary_pi +
              sm$Intercept_pi[i] +
              sm$beta_partnerNONE_pi[i] * partnerNONE_pi +
              sm$beta_block_pi[i] * block_pi
    glm_pi <- inv.logit( glm_pi )
    for(j in 1:N_pi)
        bin_sim[i,j] = rbinom(1, 1, prob = glm_pi[j])
}
hist(apply(bin_sim, 1, sum)); abline(v = sum(arabi_data$iszero_silique), col = "red")
post_zeros = mean(apply(bin_sim, 1, sum))
silique_sim = array(0, c(replicates, N_mu))
for(i in 1:replicates){
    vary_mu <- sm$vary_RIL[i, RIL_mu,3] +
               sm$vary_RIL[i, RIL_mu,4] * partnerNONE_mu
    glm_mu <- vary_mu +
              sm$Intercept_mu[i] +
              sm$beta_partnerNONE_mu[i] * partnerNONE_mu +
              sm$beta_block_mu[i] * block_mu
    for ( j in 1:N_mu )
        silique_sim[i,j] = rnorm(1, glm_mu[j], sm$sigma[i])
}
par(mfrow=c(2, 3))
hist((apply(silique_sim, 1, min))) ; abline(v =  min(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, mean))); abline(v = mean(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, max))) ; abline(v =  max(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, sum))) ; abline(v =  sum(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, sd)))  ; abline(v =   sd(((arabi_data$silique[!mask_0]))), col = "red")

partner = arabi_data$partner[!mask_0]
block = arabi_data$block[!mask_0]
form = "scale(data) ~ block + partner + (0 + partner|RIL_mu)"
extractHerit = function(data, lme_formula = form){
    x = (lmer(as.formula(lme_formula)))
    g_var = diag(VarCorr(x)[[1]])
    r_var = attributes(VarCorr(x))$sc^2
    g_var/(g_var + r_var)
}
silique_herit = apply(silique_sim, 1, extractHerit)
apply(silique_herit, 1, quantile, c(0.025, 0.5, 0.975))
apply(silique_herit, 1, mean)
