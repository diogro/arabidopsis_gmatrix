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

prior <- list(R = list(R1 = list(V = diag(4), nu = 1, fix = 1)),
              G = list(G1 = list(V = diag(1), nu = 0.002),
                       G2 = list(V = diag(1), nu = 0.002)))
model = MCMCglmm(branch ~  1,
                 random = ~partner:RIL + partner:block,
                 rcov   = ~us(trait:partner):ID,
                 family = "zipoisson",
                 verbose = TRUE,
                 prior = prior,
                 data = arabi_data)
