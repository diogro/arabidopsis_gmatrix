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
#arabi_data = arabi_data[arabi_data$height > 0,]

arabi_data$weight  <- sqrt(arabi_data$weight)

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_x")

traits = c('weight', 'height', 'silique', 'branch')
names_g = paste0(traits, rep(c('D', 'S'), each = 4))

prior <- list(R = list(R1 = list(V = diag(1), nu = 0.002, fix = 1)),
              G = list(G1 = list(V = diag(1), nu = 2)))
model_b = MCMCglmm(branch ~ trait:partner - 1,
                   random = ~trait:RIL,
                   rcov   = ~trait:units,
                   family = "zapoisson",
                   verbose = FALSE,
                   prior = prior,
                   nitt = 30000, burnin = 20000, thin = 10,
                   pl = TRUE,
                   data = arabi_data)
summary(model_b)

prior <- list(R = list(R1 = list(V = diag(1), nu = 0.002, fix = 1)),
              G = list(G1 = list(V = diag(1), nu = 0.002)))
model_s = MCMCglmm(silique ~ trait:partner - 1,
                   random = ~trait:RIL,
                   rcov   = ~trait:units,
                   family = "hupoisson",
                   verbose = FALSE,
                   prior = prior,
                   saveX = TRUE,
                   nitt = 30000, burnin = 20000, thin = 10,
                   pl = TRUE,
                   data = arabi_data)
summary(model_s)
