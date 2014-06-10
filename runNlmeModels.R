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

