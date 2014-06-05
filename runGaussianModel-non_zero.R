if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(reshape2)){install.packages("reshape2"); library(reshape2)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
#if(!require(lme4)){install.packages("lme4"); library(lme4)}
#if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(gridExtra)){install.packages("gridExtra"); library(gridExtra)}
#if(!require(gtools)){install.packages("gtools"); library(gtools)}
if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}

#################################
# Functions and definitions
#################################

#traits = c('weight', 'height', 'silique','branch')
traits = c('weight', 'height', 'silique')

num_traits = length(traits)
names_g = paste0(traits, rep(c('D', 'S'), each = num_traits))
padRmatrix <- function(x) {
    R = array(0, c(2*num_traits, 2*num_traits))
    R[1:num_traits, 1:num_traits] = x[,1:num_traits]
    R[(num_traits+1):(2*num_traits), (num_traits+1):(2*num_traits)] = x[,(num_traits+1):(2*num_traits)]
    dimnames(R) = list(names_g, names_g)
    return(R)
}
find_CI = function(x, prob = 0.95){
    n = length(x)
    xs = sort(x)
    nint = floor(prob*n)
    lowest_int = abs(xs[n] - xs[1])
    #print(lowest_int)
    for(i in 1:(n-nint)){
        current_int = abs(xs[i] - xs[i+nint])
        if(current_int <= lowest_int){
            lowest_int = current_int
            pos = i
        }
    }
    return(c(xs[pos], xs[pos+nint]))
}

#################################
# reading data
#################################

raw_arabi_data <- read.csv2("./data/raw_data.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT3)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

#####################################
# Removing missing and weird zeros
#####################################

arabi_data$flower[is.na(arabi_data$flower)] <- 0
arabi_data = arabi_data[complete.cases(arabi_data),]

arabi_data = arabi_data[arabi_data$flower > 0,]
arabi_data = arabi_data[arabi_data$height > 0,]

#####################################
# Factors and data transformations
#####################################

arabi_data$RIL   = as.factor(arabi_data$RIL)
arabi_data$block = as.factor(arabi_data$block)

arabi_data$weight  <- sqrt(arabi_data$weight)
arabi_data$silique <- sqrt(arabi_data$silique)
arabi_data$branch  <- sqrt(arabi_data$branch)

##########################################
# Scaling for heritabilitie estimations
##########################################

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
arabi_data$branch_std = arabi_data$branch
arabi_data$branch_std[ mask_partner] = scale(arabi_data$branch[ mask_partner])
arabi_data$branch_std[!mask_partner] = scale(arabi_data$branch[!mask_partner])

#####################
# Sanity check plot
#####################

m_arabi_data = melt(arabi_data, id.vars = c('partner', 'block', 'ID', 'RIL'))
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 5, scale = "free")

####################################
# MCMCglmm with all non-zero traits
###################################

arabi_data$partner = factor(arabi_data$partner)
prior = list(R = list(R1 = list(V = diag(num_traits), n = 0.002),
                      R2 = list(V = diag(num_traits), n = 0.002)),
             G = list(G1 = list(V = diag(2*num_traits) * 0.002, n = 2*num_traits + 1),
                      G2 = list(V = diag(num_traits) * 0.002, n = num_traits + 1)))
model_formula = paste0("cbind(",paste(paste0(traits, "_std"), collapse=','), ") ~ partner:trait - 1")
arabi_model = MCMCglmm(as.formula(model_formula),
                       random = ~us(trait:partner):RIL + us(trait):block,
                       rcov   = ~us(trait:at.level(partner, "L")):RIL:units +
                                 us(trait:at.level(partner, "NONE")):RIL:units,
                       family = rep("gaussian", num_traits),
                       verbose = TRUE,
                       nitt = 103000, burnin = 3000, thin = 10,
                       prior = prior,
                       data = arabi_data)
dimnames(arabi_model$Sol)[[2]] = gsub('trait', '', dimnames(arabi_model$Sol)[[2]])
dimnames(arabi_model$Sol)[[2]] = gsub('partner', '', dimnames(arabi_model$Sol)[[2]])
dimnames(arabi_model$Sol)[[2]] = gsub('L', 'D', dimnames(arabi_model$Sol)[[2]])
dimnames(arabi_model$Sol)[[2]] = gsub('NONE', 'S', dimnames(arabi_model$Sol)[[2]])
dimnames(arabi_model$Sol)[[2]] = gsub("([DS]):(.*)", "\\2\\1", dimnames(arabi_model$Sol)[[2]], perl=TRUE)

#####################################
# Extracting variance components
#####################################

Gs = array(arabi_model$VCV[,grep("RIL", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, 2*num_traits, 2*num_traits))
Bs = array(arabi_model$VCV[,grep("block", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, num_traits, num_traits))
Rs = array(arabi_model$VCV[,grep("at.level", dimnames(arabi_model$VCV)[[2]])], dim = c(10000, num_traits, 2*num_traits))
Rs = aaply(Rs, 1, padRmatrix)
corr_Gs = aaply(Gs, 1, cov2cor)
corr_Rs = aaply(Rs, 1, cov2cor)
G_mcmc = apply(Gs, 2:3, mean)
B_mcmc = apply(Bs, 2:3, mean)
R_mcmc = apply(Rs, 2:3, mean)
corr_G = apply(corr_Gs, 2:3, mean)
corr_R = apply(corr_Rs, 2:3, mean)
G_mcmc_conf = apply(Gs, 2:3, quantile, c(0.025, 0.975))
G_mcmc_conf = aperm(G_mcmc_conf, c(2, 3, 1))
containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
significant = !aaply(G_mcmc_conf, 1:2, containsZero)
dimnames(G_mcmc) = dimnames(R_mcmc) = dimnames(corr_Rs)[2:3] = dimnames(corr_Gs)[2:3] = dimnames(G_mcmc_conf)[1:2] = dimnames(corr_R) = dimnames(corr_G) = list(names_g, names_g)
summary(arabi_model)

######################
# Heritabilities
######################

herit = summary(arabi_model)$Gcovariances[seq(1, 4*num_traits*num_traits, 2*num_traits+1),1:3]
herit <- data.frame(trait   = factor(rep(traits, 2), levels = traits),
                    partner = factor(rep(c('D', 'S'), each  = num_traits), levels = c('S', 'D')),
                    herit   = herit[,1],
                    lower   = herit[,2],
                    upper   = herit[,3], row.names = NULL)
herit_plot = ggplot(herit, aes(partner, herit)) +
geom_point() + geom_errorbar(aes(ymin=lower, ymax = upper)) +
theme_classic(base_size = 15) + labs(y = 'heritabilities', x = 'trait') + facet_wrap(~trait, scale="free_y", nrow = 1)
#ggsave("~/Desktop/heritabilities_arabi.png", herit_plot)

######################
# Genetic Correlations
######################

gen_corrs = rbind(find_CI(corr_Gs[,2,1]),
                  find_CI(corr_Gs[,3,1]),
                  find_CI(corr_Gs[,3,2]),
                  find_CI(corr_Gs[,5,4]),
                  find_CI(corr_Gs[,6,4]),
                  find_CI(corr_Gs[,6,5]))
gen_corrs = data.frame(value = rowMeans(gen_corrs),
                       lower = gen_corrs[,1],
                       upper = gen_corrs[,2],
                       trait = c('weight_height', 'weight_silique', 'silique_height'),
                       partner = factor(rep(c('D', 'S'), each = num_traits), levels = c('S', 'D')))
gen_corrs_plot = ggplot(gen_corrs, aes(partner, value)) +
geom_point() + geom_errorbar(aes(ymin=lower, ymax = upper)) +
theme_classic(base_size = 15) + labs(y = 'genetic correlations', x = 'trait') + facet_wrap(~trait, scale="free_y", nrow = 1)
#ggsave("~/Desktop/genetic_correlations_arabi.png", gen_corrs_plot)

######################
# Poterior checks
######################

traits_std = paste0(paste0(traits, "_std"), rep(c('D', 'S'), each = num_traits))
n_RIL = length(unique(arabi_data$RIL))
n_block = length(unique(arabi_data$block))
N = length(arabi_data$RIL)
n_rep = 100
RIL_vector = as.numeric(arabi_data$RIL)
block_vector = as.numeric(arabi_data$block)
sim_array = array(0, dim = c(n_rep, N, 2*num_traits))
for(index in 1:n_rep){
    RIL_effect = array(0, dim = c(n_RIL, 2*num_traits))
    block_effect = array(0, dim = c(n_block, num_traits))
    for(RIL in unique(RIL_vector))
        RIL_effect[RIL,] = rmvnorm(1, sigma = Gs[index,,])
    for(block in unique(block_vector))
        block_effect[block,] = rmvnorm(1, sigma = Bs[index,,])
    for(i in 1:N){
        sim_array[index,i,] = arabi_model$Sol[index, traits_std] +
                              RIL_effect[RIL_vector[i],] +
                              rep(block_effect[block_vector[i]], 2) +
                              rmvnorm(1, sigma = Rs[index,,])
    }
}
dimnames(sim_array) = list(1:n_rep, 1:N, traits_std)

cast_phen = arabi_data
cast_phen$partner = as.character(levels(cast_phen$partner)[cast_phen$partner])
cast_phen[cast_phen == 'NONE'] <- 'S'
cast_phen[cast_phen == 'L']    <- 'D'
m.data = melt(cast_phen, id.vars= c('partner', 'RIL', 'ID', 'block'))
cast_phen = dcast(m.data, RIL~partner+variable, mean)
names(cast_phen) = gsub("([DS])_(.*)", "\\2\\1", names(cast_phen), perl=TRUE)

sim_strains = data.frame(sim_array[sample(1:n_rep, 1),,])
plot_function = function(trait1, trait2){
  out_plot = ggplot(sim_strains, aes_string(x = trait1, y = trait2, group = 1)) + geom_point(alpha = 0.3) + theme_classic(base_size = 20)
  out_plot = out_plot + geom_smooth(method="lm", color = 'black')
  out_plot = out_plot + geom_point(data = cast_phen, aes_string(x = trait1, y = trait2, group = 1), color = 'red', alpha = 0.3)
}
plots = list()
for(i in names(sim_strains)){
    plots[[i]] = list()
    for (j in names(sim_strains))
        if(i!=j) plots[[i]][[j]] = plot_function(i, j)
    names(plots[[i]]) = names(sim_strains)[names(sim_strains) != i]
}
names(plots) = names(sim_strains)
names(plots[[1]])
names(plots[[2]])
names(plots[[3]])
names(plots[[4]])
#png("~/Desktop/simulated_arabiopsis.png", heigh = 720, width = 1080)
grid.arrange(plots[['weight_stdD']][['height_stdD']], plots[['weight_stdD']][['silique_stdD']], plots[['height_stdD']][['silique_stdD']],
             plots[['weight_stdS']][['height_stdS']], plots[['weight_stdS']][['silique_stdS']], plots[['height_stdS']][['silique_stdS']], ncol = 3)
#dev.off()


obs_min = melt(sapply(cast_phen[traits_std], min, na.rm = T))
obs_min$variable = rownames(obs_min)
obs_max = melt(sapply(cast_phen[traits_std], max, na.rm = T))
obs_max$variable = rownames(obs_max)
sim_min = melt(adply(sim_array, 1, function(x) sapply(data.frame(x), min)))
sim_max = melt(adply(sim_array, 1, function(x) sapply(data.frame(x), max)))

ggplot(sim_min, aes(value)) + geom_histogram(binwidth = 0.1) + geom_vline(data = obs_min, aes(xintercept = value)) + facet_wrap(~variable)
ggplot(sim_max, aes(value)) + geom_histogram(binwidth = 0.1) + geom_vline(data = obs_max, aes(xintercept = value)) + facet_wrap(~variable)
