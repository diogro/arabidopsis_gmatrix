library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)
library(glmer2stan)
library(MCMCglmm)

raw_arabi_data <- read.csv2("./raw_data.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT3)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

arabi_data$flower[is.na(arabi_data$flower)] <- 0
arabi_data = arabi_data[complete.cases(arabi_data),]

arabi_data = arabi_data[arabi_data$flower > 0,]
arabi_data = arabi_data[arabi_data$height > 0,]

#arabi_data$weight <- scale(log(arabi_data$weight))
#arabi_data$weight <- scale(sqrt(arabi_data$weight))
arabi_data$weight <- sqrt(arabi_data$weight)


mask_0 = arabi_data$silique == 0
arabi_data$silique[!mask_0]  <- scale(log(arabi_data$silique[!mask_0]))

plot(silique~weight, arabi_data)
plot(silique~height, arabi_data)
table(arabi_data$flower, arabi_data$partner)

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 5, scale = "free")

###################
## Silique
###################

N           = dim(arabi_data)[1]
silique     = arabi_data$silique
partnerNONE = as.integer(as.factor(arabi_data$partner)) - 1
RIL         = as.integer(as.factor(arabi_data$RIL))
block       = as.integer(as.factor(arabi_data$block))
N_RIL       = length(unique(arabi_data$RIL))
N_block     = length(unique(arabi_data$block))
silique_data<- list(N           = N,
                    silique     = silique,
                    partnerNONE = partnerNONE,
                    RIL         = RIL,
                    bloob       = block,
                    N_RIL       = N_RIL,
                    N_block     = N_block)
silique_stan_model = stan(file = './silique.stan', data = silique_data, chain=1)

sm = extract(silique_stan_model, permuted = TRUE)
replicates = dim(sm$vary_RIL)[1]
silique_sim = array(0, c(replicates, N))
for(i in 1:replicates){
    vary <- sm$vary_RIL[i, RIL,1] +
               sm$vary_RIL[i, RIL,2] * partnerNONE +
               sm$vary_block[i, block]
    glm <- vary + sm$Intercept[i]
    for ( j in 1:N )
        silique_sim[i,j] = rnorm(1, glm[j], sm$sigma[i])
}
par(mfrow=c(2, 3))
hist((apply(silique_sim, 1, min))) ; abline(v =  min(((arabi_data$silique))), col = "red")
hist((apply(silique_sim, 1, mean))); abline(v = mean(((arabi_data$silique))), col = "red")
hist((apply(silique_sim, 1, max))) ; abline(v =  max(((arabi_data$silique))), col = "red")
hist((apply(silique_sim, 1, sd)))  ; abline(v =   sd(((arabi_data$silique))), col = "red")
hist(silique_sim[1,]); hist(arabi_data$silique)

extractHerit = function(x) diag(cov(cbind(x[,1], x[,1]+x[,2])))
herit_silique = t(apply(sm$vary_RIL, 1, extractHerit)/rbind(sm$sigma, sm$sigma))
dimnames(herit_silique) = list(NULL, c("L", "NONE"))
colMeans(herit_silique)
boxplot(herit_silique)

silique_model = lmer(silique ~ 1 + (0 + partner|RIL) + (1|block),
                     data = arabi_data)
summary(silique_model)
varRIL = diag(VarCorr(silique_model)$RIL)
varRep = rep(VarCorr(silique_model)$block[1], 2)
varRes = rep(attributes(VarCorr(silique_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))

###################
## Weight
###################

the_formula <- list(weight ~ 1 + (partner|RIL) + (1|block))
weight_model = glmer2stan(the_formula, data=arabi_data,
                           family="gaussian",
                           sample = FALSE, calcDIC = FALSE)
write(weight_model$model, file = "weight.stan")

N           = dim(arabi_data)[1]
weight      = arabi_data$weight
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

weight_model = lmer(weight ~ 1 + (0 + partner|RIL) + (1|block), data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(weight_model)
varRIL = diag(VarCorr(weight_model)$RIL)
varRep = rep(VarCorr(weight_model)$block[1], 2)
varRes = rep(attributes(VarCorr(weight_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))


###################
## Height
###################

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

multi_model = lmer(value ~ variable + (0 + variable:partner|RIL) + (variable|block), data = m_arabi_data)
summary(multi_model)
VarCorr(multi_model)
