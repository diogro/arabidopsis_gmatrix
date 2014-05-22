library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)
library(glmer2stan)

raw_arabi_data <- read.csv2("./raw_data.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT3)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

arabi_data$flower[is.na(arabi_data$flower)] <- 0
arabi_data = arabi_data[complete.cases(arabi_data),]

arabi_data = arabi_data[arabi_data$flower > 0,]
arabi_data = arabi_data[arabi_data$height > 0,]
arabi_data$weight  <- scale(log(arabi_data$weight))
arabi_data$height  <- scale(as.numeric(arabi_data$height/10))
mask_0 = arabi_data$silique == 0
arabi_data$silique[!mask_0]  <- scale(log(arabi_data$silique[!mask_0]))
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

branch_model = glmer2stan(branch ~ 1 + (partner|RIL) + (1|block),
                          data=arabi_data,
                          family="poisson"
                          sample = FALSE, calcDIC = FALSE)
write(branch_model$model, file = "branch.stan")

branch_data <- list(N = dim(arabi_data)[1],
                    branch = arabi_data$branch,
                    partnerNONE = as.integer(as.factor(arabi_data$partner)),
                    RIL = as.integer(as.factor(arabi_data$RIL)),
                    N_RIL = length(unique(arabi_data$RIL)))
branch_stan_model = stan(file = './branch.stan', data = branch_data, chain=1)
bm = extract(branch_stan_model, permuted = TRUE)


the_formula <- list(iszero_silique ~ 1 + (partner|RIL) + (1|block),
                    silique        ~ 1 + (partner|RIL) + (1|block))
silique_model = glmer2stan(the_formula, data=arabi_data,
                           family=list("binomial","gaussian"),
                           sample = FALSE, calcDIC = FALSE)
write(silique_model$model, file = "silique.stan")

N_pi              = dim(arabi_data)[1]
N_mu              = sum(!mask_0)
iszero_silique_pi = as.integer(mask_0)
partnerNONE_pi    = as.integer(as.factor(arabi_data$partner)) - 1
RIL_pi            = as.integer(as.factor(arabi_data$RIL))
block_pi          = as.integer(as.factor(arabi_data$block))
bin_total_pi      = rep(1,dim(arabi_data)[1])
silique_mu        = arabi_data$silique[!mask_0]
partnerNONE_mu    = partnerNONE_pi[!mask_0]
RIL_mu            = as.integer(as.factor(arabi_data$RIL))[!mask_0]
block_mu          = block_pi[!mask_0]
N_RIL             = length(unique(arabi_data$RIL))
N_block           = length(unique(arabi_data$block))
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
                     N_RIL             = N_RIL,
                     N_block           = N_block)
silique_stan_model = stan(file = './silique.stan', data = silique_data, chain=1)

sm = extract(silique_stan_model, permuted = TRUE)
replicates = dim(sm$vary_RIL)[1]
bin_sim = array(0, c(replicates, N_pi))
for(i in 1:replicates){
    vary_pi <- sm$vary_RIL[i, RIL_pi,1] +
               sm$vary_RIL[i, RIL_pi,2] * partnerNONE_pi +
               sm$vary_block[i, block_pi, 1]
    glm_pi <- vary_pi + sm$Intercept_pi[i]
    glm_pi <- inv.logit( glm_pi )
    for(j in 1:N_pi)
        bin_sim[i,j] = rbinom(1, 1, prob = glm_pi[j])
}
hist(apply(bin_sim, 1, sum)); abline(v = sum(mask_0), col = "red")
post_zeros = mean(apply(bin_sim, 1, sum))
silique_sim = array(0, c(replicates, N_mu))
for(i in 1:replicates){
    vary_mu <- sm$vary_RIL[i, RIL_mu,3] +
               sm$vary_RIL[i, RIL_mu,4] * partnerNONE_mu +
               sm$vary_block[i, block_mu, 2]
    glm_mu <- vary_mu + sm$Intercept_mu[i]
    for ( j in 1:N_mu )
        silique_sim[i,j] = rnorm(1, glm_mu[j], sm$sigma[i])
}
par(mfrow=c(2, 3))
hist((apply(silique_sim, 1, min))) ; abline(v =  min(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, mean))); abline(v = mean(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, max))) ; abline(v =  max(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, sum))) ; abline(v =  sum(((arabi_data$silique[!mask_0]))), col = "red")
hist((apply(silique_sim, 1, sd)))  ; abline(v =   sd(((arabi_data$silique[!mask_0]))), col = "red")

extractHerit = function(x) diag(cov(cbind(x[,3], x[,3]+x[,2])))
silique_herit = t(apply(sm$vary_RIL_mu, 1, extractHerit)/rbind(sm$sigma, sm$sigma))
dimnames(silique_herit) = list(NULL, c("L", "NONE"))
rowMeans(silique_herit)
boxplot(silique_herit)

the_formula <- list(weight ~ 1 + (partner|RIL) + (1|block))
weight_model = glmer2stan(the_formula, data=arabi_data,
                           family="gaussian",
                           sample = FALSE, calcDIC = FALSE)
write(weight_model$model, file = "weight.stan")

N           = dim(arabi_data)[1]
weight      = arabi_data$weight[,1]
partnerNONE = as.integer(as.factor(arabi_data$partner)) - 1
RIL         = as.integer(as.factor(arabi_data$RIL))
block       = as.integer(as.factor(arabi_data$block))
N_RIL       = length(unique(arabi_data$RIL))
N_block     = length(unique(arabi_data$block))
weight_data <- list(N           = N,
                    weight      = weight,
                    partnerNONE = partnerNONE,
                    RIL         = RIL,
                    bloob       = block,
                    N_RIL       = N_RIL,
                    N_block     = N_block)
weight_stan_model = stan(file = './weight.stan', data = weight_data, chain=1)

wm = extract(weight_stan_model, permuted = TRUE)
replicates = dim(wm$vary_RIL)[1]
weight_sim = array(0, c(replicates, N))
for(i in 1:replicates){
    vary <- wm$vary_RIL[i, RIL,1] +
               wm$vary_RIL[i, RIL,2] * partnerNONE +
               wm$vary_bloob[i, block]
    glm <- vary + wm$Intercept[i]
    for ( j in 1:N )
        weight_sim[i,j] = rnorm(1, glm[j], wm$sigma[i])
}
par(mfrow=c(2, 3))
hist((apply(weight_sim, 1, min))) ; abline(v =  min(((arabi_data$weight))), col = "red")
hist((apply(weight_sim, 1, mean))); abline(v = mean(((arabi_data$weight))), col = "red")
hist((apply(weight_sim, 1, max))) ; abline(v =  max(((arabi_data$weight))), col = "red")
hist((apply(weight_sim, 1, sd)))  ; abline(v =   sd(((arabi_data$weight))), col = "red")
hist(weight_sim[1,]); hist(arabi_data$weight)

extractHerit = function(x) diag(cov(cbind(x[,1], x[,1]+x[,2])))
herit_weight = t(apply(wm$vary_RIL, 1, extractHerit)/rbind(wm$sigma, wm$sigma))
dimnames(herit_weight) = list(NULL, c("L", "NONE"))
colMeans(herit_weight)
boxplot(herit_weight)

weight_model = lmer(weight ~ 1 + (0 + partner|RIL) + (1|block),
                    data = arabi_data)
summary(weight_model)
varRIL = diag(VarCorr(weight_model)$RIL)
varRep = rep(VarCorr(weight_model)$block[1], 2)
varRes = rep(attributes(VarCorr(weight_model))$sc^2, 2)
h2 = varRIL/(varRIL + varRep + varRes)

the_formula <- list(height ~ 1 + (partner|RIL) + (1|block))
height_model = glmer2stan(the_formula, data=arabi_data,
                           family="gaussian",
                           sample = FALSE, calcDIC = FALSE)
write(height_model$model, file = "height.stan")

N           = dim(arabi_data)[1]
height      = arabi_data$height[,1]
partnerNONE = as.integer(as.factor(arabi_data$partner)) - 1
RIL         = as.integer(as.factor(arabi_data$RIL))
block       = as.integer(as.factor(arabi_data$block))
N_RIL       = length(unique(arabi_data$RIL))
N_block     = length(unique(arabi_data$block))
height_data <- list(N           = N,
                    height      = height,
                    partnerNONE = partnerNONE,
                    RIL         = RIL,
                    bloob       = block,
                    N_RIL       = N_RIL,
                    N_block     = N_block)
height_stan_model = stan(file = './height.stan', data = height_data, chain=1)
#print(height_stan_model)

hm = extract(height_stan_model, permuted = TRUE)
replicates = dim(hm$vary_RIL)[1]
height_sim = array(0, c(replicates, N))
for(i in 1:replicates){
    vary <- hm$vary_RIL[i, RIL,1] +
               hm$vary_RIL[i, RIL,2] * partnerNONE +
               hm$vary_block[i, block]
    glm <- vary + hm$Intercept[i]
    for ( j in 1:N )
        height_sim[i,j] = rnorm(1, glm[j], hm$sigma[i])
}
par(mfrow=c(2, 3))
hist((apply(height_sim, 1, min))) ; abline(v =  min(((arabi_data$height))), col = "red")
hist((apply(height_sim, 1, mean))); abline(v = mean(((arabi_data$height))), col = "red")
hist((apply(height_sim, 1, max))) ; abline(v =  max(((arabi_data$height))), col = "red")
hist((apply(height_sim, 1, sd)))  ; abline(v =   sd(((arabi_data$height))), col = "red")
hist(height_sim[1,]); hist(arabi_data$height)

extractHerit = function(x) diag(cov(cbind(x[,1], x[,1]+x[,2])))
herit_height = t(apply(hm$vary_RIL, 1, extractHerit)/rbind(hm$sigma, hm$sigma))
dimnames(herit_height) = list(NULL, c("L", "NONE"))
colMeans(herit_height)
boxplot(herit_height)

height_model = lmer(height ~ 1 + (0 + partner|RIL) + (1|block), data = arabi_data)
summary(height_model)
varRIL = diag(VarCorr(height_model)$RIL)
varRep = rep(VarCorr(height_model)$block[1], 2)
varRes = rep(attributes(VarCorr(height_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))
