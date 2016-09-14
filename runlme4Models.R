if(!require(ggplot2))  { install.packages("ggplot2")  ;  library(ggplot2)  }
if(!require(reshape2)) { install.packages("reshape2") ;  library(reshape2) }
if(!require(dplyr))    { install.packages("dplyr")    ;  library(dplyr)    }
if(!require(lme4))     { install.packages("lme4")     ;  library(lme4)     }

#################################
# Functions and definitions
#################################

dir.create(file.path("./figures"), showWarnings = FALSE)

#traits = c('weight', 'height', 'silique','branch')
traits = c('weight', 'height', 'silique')

num_traits = length(traits)
names_g = paste0(traits, rep(c('D', 'S'), each = num_traits))

padRmatrix = function(x) {          # create square residual matrix
    R = array(0, c(2*num_traits, 2*num_traits))
    R[1:num_traits, 1:num_traits] = x[,1:num_traits]
    R[(num_traits+1):(2*num_traits), (num_traits+1):(2*num_traits)] = x[,(num_traits+1):(2*num_traits)]
    dimnames(R) = list(names_g, names_g)
    return(R)
}
find_CI = function(x, prob = 0.95){  # create credible intervals
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
checkStat = function(stat, title = ''){ # simple posterior checks
    obs_stat = melt(sapply(cast_phen[traits_std], stat, na.rm = T))
    obs_stat$variable = rownames(obs_stat)
    sim_stat = melt(adply(sim_array, 1, function(x) sapply(data.frame(x), stat)))
    ggplot(sim_stat, aes(value)) + geom_histogram() + geom_vline(data = obs_stat, aes(xintercept = value)) + facet_wrap(~variable) + ggtitle(title)
}

#################################
# reading data
#################################

raw_arabi_data = read.csv2("./data/raw_data.csv", as.is = T)
arabi_data = select(raw_arabi_data, ID, RIL, Block, Partner, HEIGHT, WEIGHT, SILIQUEN, NODEN, BOLT3)
names(arabi_data) = c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch", "flower")

#####################################
# Removing missing and weird zeros
#####################################

arabi_data$flower[is.na(arabi_data$flower)] = 0
arabi_data = arabi_data[complete.cases(arabi_data),]

arabi_data = arabi_data[arabi_data$flower > 0,]
arabi_data = arabi_data[arabi_data$height > 0,]

#####################################
# Factors and data transformations
#####################################

arabi_data$RIL     = as.factor(arabi_data$RIL)
arabi_data$block   = as.factor(arabi_data$block)

arabi_data$partner[arabi_data$partner == 'NONE'] = 'S'
arabi_data$partner[arabi_data$partner == 'L']    = 'D'
arabi_data$partner = as.factor(arabi_data$partner)

arabi_data$weight  = sqrt(arabi_data$weight)
arabi_data$silique = sqrt(arabi_data$silique)
arabi_data$branch  = sqrt(arabi_data$branch)

##########################################
# Scaling for heritability estimation
##########################################

mask_partner = arabi_data$partner == "D"
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

##########################
# Multivariate Model
##########################

#m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight_std, height_std, silique_std), id.vars = c('partner', 'block', 'ID', 'RIL'))
m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight, height, silique), id.vars = c('partner', 'block', 'ID', 'RIL'))
multi_model = lmer(value ~ variable*partner + variable*block + (0 + variable:partner|RIL),
                   data = m_arabi_data, REML = FALSE, na.action = 'na.omit', control=lmerControl(check.conv.singular="warning"))
summary(multi_model)
VarCorr(multi_model)
VarCorr(multi_model)$RIL


#Muito mau comportado e incoerente com os univariados. Suspeito que seja
#falta de replicação dentro de linhagens

m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight_std, height_std, silique_std), id.vars = c('partner', 'block', 'ID', 'RIL'))
#m_arabi_data = melt(select(arabi_data, ID, block, RIL, partner, weight, height, silique), id.vars = c('partner', 'block', 'ID', 'RIL'))
multi_model = lmer(value ~ variable + variable:partner + variable:block + (0 + variable:partner|RIL),
                   data = m_arabi_data, na.action = 'na.omit', control=lmerControl(check.conv.singular="warning"))
summary(multi_model)
VarCorr(multi_model)
VarCorr(multi_model)$RIL
diag(VarCorr(multi_model)$RIL)

##########################
# Univariate Models
##########################

###################
## Silique
###################

silique_model = lmer(silique ~ partner*block + (0 + partner|RIL),
                    data = arabi_data, na.action = 'na.omit')
summary(silique_model)
varRIL = diag(VarCorr(silique_model)$RIL)
varRes = rep(attributes(VarCorr(silique_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRes))

silique_model = lmer(silique_std ~ partner*block + (0 + partner|RIL),
                     data = arabi_data, na.action = 'na.omit')
summary(silique_model)
varRIL = diag(VarCorr(silique_model)$RIL)
varRes = rep(attributes(VarCorr(silique_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRes))

silique_model_D = lmer(silique_std ~ block + (1|RIL),
                       data = filter(arabi_data, partner == "D"))
summary(silique_model_D)
silique_model_S = lmer(silique_std ~ block + (1|RIL),
                       data = filter(arabi_data, partner == "S"))
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

# un-scaled per treatment
weight_model_D = lmer(weight ~ (1|RIL),
                       data = filter(arabi_data, partner == "D"))
summary(weight_model_D)
weight_model_S = lmer(weight ~ (1|RIL),
                       data = filter(arabi_data, partner == "S"))
summary(weight_model_S)
varRIL = c("D" = VarCorr(weight_model_D)$RIL, "S" = VarCorr(weight_model_S)$RIL)
(h2 = varRIL)

# scaled per treatment
weight_model_D = lmer(weight_std ~ (1|RIL),
                       data = filter(arabi_data, partner == "D"))
summary(weight_model_D)
weight_model_S = lmer(weight_std ~ (1|RIL),
                       data = filter(arabi_data, partner == "S"))
summary(weight_model_S)
varRIL = c("D" = VarCorr(weight_model_D)$RIL, "S" = VarCorr(weight_model_S)$RIL)
(h2 = varRIL)

###################
## Height
###################

height_model = lmer(height ~ partner + (0 + partner|RIL) + (1|block),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(height_model)
varRIL = diag(VarCorr(height_model)$RIL)
varRes = rep(attributes(VarCorr(height_model))$sc^2, 2)
(h2 = varRIL/(varRIL + varRep + varRes))

height_model = lmer(height_std ~ partner + (0 + partner|RIL) + (1|block),
                    data = arabi_data, REML = FALSE, na.action = 'na.omit')
summary(height_model)
varRIL = diag(VarCorr(height_model)$RIL)
varRes = rep(attributes(VarCorr(height_model))$sc^2, 2)
(h2 = varRIL)

