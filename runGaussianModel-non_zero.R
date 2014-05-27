library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)
library(MCMCglmm)

raw_arabi_data <- read.csv2("./raw_data.csv")
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

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 5, scale = "free")

traits = c('weight', 'height', 'silique')
names_g = paste0(traits, rep(c('D', 'S'), each = 3))
padRmatrix <- function(x) {
    R = array(0, c(2*num_traits, 2*num_traits))
    R[1:num_traits, 1:num_traits] = x[,1:num_traits]
    R[(num_traits+1):(2*num_traits), (num_traits+1):(2*num_traits)] = x[,(num_traits+1):(2*num_traits)]
    dimnames(R) = list(names_g, names_g)
    return(R)
}

#############
# MCMCglmm with all non-zero traits
#############

arabi_data$partner = factor(arabi_data$partner)
num_traits = 3
prior = list(R = list(R1 = list(V = diag(num_traits), n = 0.002),
                      R2 = list(V = diag(num_traits), n = 0.002)),
             G = list(G1 = list(V = diag(2*num_traits) * 0.02, n = 2*num_traits+1),
                      G2 = list(V = diag(num_traits) * 0.02, n = num_traits+1)))
arabi_model = MCMCglmm(cbind(weight_std, height_std, silique_std) ~ partner:trait - 1,
                       random = ~us(trait:partner):RIL + us(trait):block,
                       rcov   = ~us(trait:at.level(partner, "L")):units +
                                 us(trait:at.level(partner, "NONE")):units,
                       family = rep("gaussian", num_traits),
                       verbose = TRUE,
                       nitt = 103000, burnin = 3000, thin = 10,
                       prior = prior,
                       data = arabi_data)
write.csv(arabi_data,"./data_clean.csv")

Gs = array(arabi_model$VCV[,grep("RIL", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, 2*num_traits, 2*num_traits))
Bs = array(arabi_model$VCV[,grep("block", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, num_traits, num_traits))
Rs = array(arabi_model$VCV[,grep("at.level", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, num_traits, 2*num_traits))
Rs = aaply(Rs, 1, padRmatrix)
G_mcmc = apply(Gs, 2:3, mean)
B_mcmc = apply(Bs, 2:3, mean)
R_mcmc = apply(Rs, 2:3, mean)
corr_G = cov2cor(G_mcmc)
corr_R = cov2cor(R_mcmc)
G_mcmc_conf = apply(Gs, 2:3, quantile, c(0.025, 0.975))
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
dimnames(G_mcmc) = dimnames(R_mcmc) = dimnames(G_mcmc_conf)[1:2] = dimnames(corr_R) = dimnames(corr_G) = list(names_g, names_g)
sim_strains = adply(1:1000, 1, function(index) mvtnorm::rmvnorm(1, arabi_model$Sol[index,], Gs[index,,]))
names(sim_strains) = gsub('trait', '', names(sim_strains))
names(sim_strains) = gsub('partner', '', names(sim_strains))
names(sim_strains) = gsub('L', 'D', names(sim_strains))
names(sim_strains) = gsub('NONE', 'S', names(sim_strains))
names(sim_strains) = gsub("([DS]):(.*)", "\\2\\1", names(sim_strains), perl=TRUE)
herit = summary(arabi_model)$Gcovariances[seq(1, 4*num_traits*num_traits, 2*num_traits+1),1:3]
herit <- data.frame(trait   = factor(rep(traits, 2), levels = traits),
                    partner = factor(rep(c('D', 'S'), each  = 3), levels = c('S', 'D')),
                    herit   = herit[,1],
                    lower   = herit[,2],
                    upper   = herit[,3], row.names = NULL)
summary(arabi_model)

herit_plot = ggplot(herit, aes(partner, herit)) +
geom_point() + geom_errorbar(aes(ymin=lower, ymax = upper)) +
theme_classic(base_size = 15) + labs(y = 'heritabilities', x = 'trait') + facet_wrap(~trait, scale="free_y")
ggsave("~/Desktop/heritabilities_arabi.png", herit_plot)

cast_phen = arabi_data
cast_phen$partner = as.character(levels(cast_phen$partner)[cast_phen$partner])
cast_phen[cast_phen == 'NONE'] <- 'S'
cast_phen[cast_phen == 'L']    <- 'D'
m.data = melt(cast_phen, id.vars= c('partner', 'RIL', 'ID'))
cast_phen = dcast(m.data, RIL~partner+variable, mean)
names(cast_phen) = gsub("([DS])_(.*)", "\\2\\1", names(cast_phen), perl=TRUE)

plot_function = function(trait1, trait2){
  out_plot = ggplot(sim_strains, aes_string(x = trait1, y = trait2, group = 1)) + geom_point(alpha = 0.3) + theme_classic(base_size = 20)
  out_plot = out_plot + geom_smooth(method="lm", color = 'black')
  out_plot = out_plot + geom_point(data = cast_phen, aes_string(x = trait1, y = trait2, group = 1), color = 'red', alpha = 0.3)
}
plots = list()
for(i in names(sim_strains)[-1]){
    plots[[i]] = list()
    for (j in names(sim_strains)[-1])
        if(i!=j) plots[[i]][[j]] = plot_function(i, j)
    names(plots[[i]]) = names(sim_strains)[names(sim_strains) != i][-1]
}
names(plots) = names(sim_strains)[-1]
names(plots[[4]])
tiff("~/Desktop/simulated_arabiopsis.tiff", heigh = 720, width = 1080)
grid.arrange(plots[[1]][[2]], plots[[1]][[4]], plots[[3]][[4]], plots[[2]][[3]], plots[[2]][[5]], plots[[4]][[5]], ncol = 3)
dev.off()
