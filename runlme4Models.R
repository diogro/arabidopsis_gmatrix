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

#raw_arabi_data <- read.csv2("./raw_data.csv")
raw_arabi_data <- read.csv("./data_clean_jason.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

arabi_data$flower[is.na(arabi_data$flower)] <- 0
arabi_data = arabi_data[complete.cases(arabi_data),]

arabi_data = arabi_data[arabi_data$flower > 0,]
arabi_data = arabi_data[arabi_data$height > 0,]

arabi_data$weight <- sqrt(arabi_data$weight)
arabi_data$silique  <- sqrt(arabi_data$silique)

mask_partner = arabi_data$partner == "L"
arabi_data$silique_std = arabi_data$silique
arabi_data$silique_std[ mask_partner] = scale(arabi_data$silique[ mask_partner])
arabi_data$silique_std[!mask_partner] = scale(arabi_data$silique[!mask_partner])
arabi_data$weight_std = arabi_data$weight
arabi_data$weight_std[ mask_partner] = scale(arabi_data$weight[ mask_partner])
arabi_data$weight_std[!mask_partner] = scale(arabi_data$weight[!mask_partner])
arabi_data$height_std = arabi_data$height
arabi_data$height_std[ mask_partner] = scale(arabi_data$height[ mask_partner])
arabi_data$height_std[!mask_partner] = scale(arabi_data$height[!mask_partner])

write.csv(arabi_data, "./data_clean.csv")

#mask_0 = arabi_data$silique == 0 | is.na(arabi_data$silique)
#arabi_data$silique[!mask_0]  <- scale(sqrt(arabi_data$silique[!mask_0]))

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

silique_model = lmer(silique ~ partner + (0 + partner|RIL),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(silique_model)
varRIL = diag(VarCorr(silique_model)$RIL)
varRep = rep(VarCorr(silique_model)$block[1], 2)
varRes = rep(attributes(VarCorr(silique_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRes))

silique_model = lmer(silique_std ~ partner + (0 + partner|RIL),
                     data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(silique_model)
varRIL = diag(VarCorr(silique_model)$RIL)
(h2 = varRIL)

silique_model_D = lmer(silique ~ 1 + (1|RIL),
                       data = filter(arabi_data, partner == "L"), na.action = 'na.omit')
summary(silique_model_D)
silique_model_S = lmer(silique ~ 1 + (1|RIL),
                       data = filter(arabi_data, partner == "NONE"), na.action = 'na.omit')
summary(silique_model_S)
varRIL = c("D" = VarCorr(silique_model_D)$RIL, "S" = VarCorr(silique_model_S)$RIL)
(h2 = varRIL)

###################
## Weight
###################

weight_model = lmer(weight ~ 1 + (0 + partner|RIL) + (1|block),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(weight_model)
varRIL = diag(VarCorr(weight_model)$RIL)
varRep = rep(VarCorr(weight_model)$block[1], 2)
varRes = rep(attributes(VarCorr(weight_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))

weight_model = lmer(weight_std ~ 1 + (0 + partner|RIL) + (1|block),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(weight_model)
varRIL = diag(VarCorr(weight_model)$RIL)
(h2 = varRIL)

###################
## Height
###################

height_model = lmer(height ~ partner + (0 + partner|RIL),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(height_model)
varRIL = diag(VarCorr(height_model)$RIL)
varRes = rep(attributes(VarCorr(height_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))

height_model = lmer(height_std ~ partner + (0 + partner|RIL),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(height_model)
varRIL = diag(VarCorr(height_model)$RIL)
varRes = rep(attributes(VarCorr(height_model))$sc^2, 2)
(h2 = varRIL)

m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight, height, silique), id.vars = c('partner', 'block', 'ID', 'RIL'))
m_arabi_data$trait_partner = paste0(m_arabi_data$partner, m_arabi_data$variable)
multi_model = lmer(value ~ trait_partner + (0 + trait_partner|RIL),
                   data = m_arabi_data, REML = FALSE, na.action = 'na.omit')
summary(multi_model)
VarCorr(multi_model)

m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight_std, height_std, silique_std), id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes( RIL, value, group = interaction(partner, variable, RIL), color = partner)) + geom_boxplot() + facet_wrap(~variable)
table(arabi_data$RIL, as.character(arabi_data$partner))
