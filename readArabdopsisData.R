library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(MCMCglmm)
library(gridExtra)
library(gtools)

raw_arabi_data <- read.csv2("./raw_data.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch")

arabi_data$weight  <- sqrt(arabi_data$weight)
arabi_data$silique <- sqrt(arabi_data$silique)
arabi_data$branch  <- sqrt(arabi_data$branch)

arabi_data = arabi_data[complete.cases(arabi_data),]
arabi_data = arabi_data[arabi_data$height > 0,]

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_x")

mask_L = arabi_data$partner == "L"
traits = c('weight', 'height', 'silique', 'branch')

model_block_D = lmer(value ~ variable + (variable|block), data = filter(m_arabi_data, partner == 'L'))
model_block_S = lmer(value ~ variable + (variable|block), data = filter(m_arabi_data, partner == 'NONE'))
m_arabi_data_std = m_arabi_data
m_arabi_data_std$block <- NULL
m_arabi_data_std$value[m_arabi_data_std$partner == 'L'] = residuals(model_block_D)
m_arabi_data_std$value[m_arabi_data_std$partner == 'NONE'] = residuals(model_block_S)
arabi_data_std = dcast(m_arabi_data_std, partner+ID+RIL~variable)

arabi_data_std$weight  <- scale(arabi_data_std$weight)
arabi_data_std$silique <- scale(arabi_data_std$silique)
arabi_data_std$branch  <- scale(arabi_data_std$branch)
arabi_data_std$height  <- scale(arabi_data_std$height)

num_traits = 4
prior = list(R = list(R1 = list(V = diag(num_traits), n = 0.002),
                      R2 = list(V = diag(num_traits), n = 0.002)),
             G = list(G1 = list(V = diag(2*num_traits) * 0.02, n = 2*num_traits+1)))
arabi_model = MCMCglmm(cbind(weight, height, silique, branch) ~ trait:partner - 1,
                       random = ~us(trait:partner):RIL,
                       rcov   = ~us(trait:at.level(partner,    "L")):units +
                                 us(trait:at.level(partner, "NONE")):units,
                       family = rep("gaussian", num_traits),
                       prior = prior,
                       data = arabi_data_std)
Gs = array(arabi_model$VCV[,1:(4*num_traits*num_traits)], dim = c(1000, 2*num_traits, 2*num_traits))
G_mcmc = apply(Gs, 2:3, mean)
corr_G = cov2cor(G_mcmc)
names_g = paste0(traits, rep(c('D', 'S'), each= 4))
dimnames(G_mcmc)  = dimnames(corr_G) = list(names_g, names_g)
#corr_G[5:8, 1:4] = cov2cor(G_mcmc[5:8, 1:4])
#corr_G[1:4, 5:8] = cov2cor(G_mcmc[1:4, 5:8])
diag(G_mcmc)
G_mcmc_conf = apply(Gs, 2:3, quantile, c(0.025, 0.975))
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
herit = data.frame(trait = traits,
                   partner = rep(c('D', 'S'), each= 4),
                   value = diag(G_mcmc),
                   upper = diag(G_mcmc_conf[,,2]),
                   lower = diag(G_mcmc_conf[,,1]))

ggplot(herit, aes(trait, value, color = partner)) +
geom_point() + geom_errorbar(aes(ymin=lower, ymax = upper)) +
theme_classic() + labs(y = 'heritabilities', x = 'trait') + facet_wrap(~partner)
