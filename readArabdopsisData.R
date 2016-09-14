if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "dplyr", "tidyr", "readr", "lme4", "ggplot2", "cowplot", "MCMCglmm", "rstan")

raw_arabi_data <- read_csv("./data/phenotypes.csv")
arabi_data <- select(raw_arabi_data, ID, RIL, Block, envi, HEIGHT, WEIGHT, SILIQUEN, NODEN)
names(arabi_data) <- c("ID", "RIL", "block", "partner", "height", "weight", "silique", "branch")

arabi_data$height[arabi_data$height == 0] = NA
arabi_data$silique[is.na(arabi_data$height)] = NA

arabi_data = arabi_data[complete.cases(arabi_data),]
arabi_data = arabi_data[arabi_data$height > 0,]

arabi_data$silique = sqrt(arabi_data$silique)
arabi_data$branch = sqrt(arabi_data$branch)
arabi_data$weight = sqrt(arabi_data$weight)

m_arabi_data = gather(arabi_data, variable, value, height:branch)
ggplot(m_arabi_data, aes(x = value, color = partner)) +
geom_histogram() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap(~variable, ncol = 4, scale = "free_x")

traits = c("height", "weight", "silique")

raw_markers = read_csv("./data/BayxSha_2_Genotypes_imputed.csv", skip= 1)[-c(1:4),-2] %>% rename(RIL = order) %>% mutate(RIL = as.integer(RIL))
raw_markers[] = lapply(raw_markers, function(x) {
                           x[x == "A"] = 1;
                           x[x == "B"] = -1;
                           x[x == "C"] = 0;
                           x[x == "D"] = NA;
                       return(as.integer(x)) } )
names(raw_markers)[2:70] = paste0("G", 1:69)
markers = semi_join(raw_markers, arabi_data, id = "RIL")
arabi_data = semi_join(arabi_data, raw_markers, id = "RIL")
RIL_levels = unique(markers$RIL)
arabi_data$RIL = as.integer(factor(arabi_data$RIL, RIL_levels))
markers = as.matrix(select(markers, G1:G69))

obs_markers = markers
miss_J = sum(is.na(obs_markers))
mis_pos = which(is.na(obs_markers), TRUE)
obs_markers[is.na(obs_markers)] = 2
