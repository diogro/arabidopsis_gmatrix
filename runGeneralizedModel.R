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

arabi_data = arabi_data[complete.cases(arabi_data),]
arabi_data = arabi_data[arabi_data$height > 0,]

arabi_data$weight  <- sqrt(arabi_data$weight)
arabi_data$silique  <- scale(sqrt(arabi_data$silique))

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_x")

traits = c('weight', 'height', 'silique', 'branch')
names_g = paste0(traits, rep(c('D', 'S'), each = 4))

prior <- list(R = list(R1 = list(V = diag(1), nu = 0.002)),
              G = list(G1 = list(V = diag(2), nu = 0.002)))
model_b = MCMCglmm(branch ~ partner - 1,
                   random = ~idh(partner):RIL,
                   rcov   = ~units,
                   family = "ordinal",
                   verbose = FALSE,
                   saveX = TRUE,
                   prior = prior,
                   nitt = 30000, burnin = 20000, thin = 10,
                   pl = TRUE,
                   data = arabi_data)
summary(model_b)
posterior.mode(model_b$VCV)
herit_D_b =  model_b$VCV[, "L.RIL"]/(model_b$VCV[, "L.RIL"] + model_b$VCV[, "units"]+1)
herit_S_b =  model_b$VCV[, "NONE.RIL"]/(model_b$VCV[, "NONE.RIL"] + model_b$VCV[, "units"]+1)
HPDinterval(herit_S_b)
HPDinterval(herit_D_b)

prior <- list(R = list(R1 = list(V = diag(2), nu = 0.002)),
              G = list(G1 = list(V = diag(2), nu = 0.002)))
model_s = MCMCglmm(silique ~ 1,
                   random = ~us(partner):RIL,
                   rcov   = ~idh(partner):units,
                   family = "gaussian",
                   verbose = FALSE,
                   prior = prior,
                   #nitt = 1100000, burnin = 100000, thin = 500,
                   data = arabi_data)
summary(model_s)
posterior.mode(model_s$VCV)
herit_D_s =  model_s$VCV[, "L:L.RIL"]/(model_s$VCV[, "L:L.RIL"] + model_s$VCV[, "L.units"]+1)
herit_S_s =  model_s$VCV[, "NONE:NONE.RIL"]/(model_s$VCV[, "NONE:NONE.RIL"] + model_s$VCV[, "NONE.units"]+1)
posterior.mode(herit_S_s)
posterior.mode(herit_D_s)

plot(branch~heigth, arabi_data)
