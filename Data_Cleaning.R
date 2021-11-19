# DATA CLEANING -----------------------------------------------------------------

setwd('C:/Users/pietr/Desktop/Bayesian Statistics/Progetto')
data <- read.csv2('dati.csv',header = T,sep = ',')
data <- data.frame(data)
colnames(data)[3] <- 'sex' # 1 = Maschio, 2 = Femmina
colnames(data)[4] <- 'eziologia' # 1 = Trauma, 2 = Vascolare
colnames(data)[5] <- 'lateralità' # 1 = dx. 2 = sx
colnames(data)[6] <- 'data'
colnames(data)[13] <- 'SLSEP' # 1 = patologico, 2 = normale
colnames(data)[17] <- 'MLSEP' # 1 = patologico, 2 = presente