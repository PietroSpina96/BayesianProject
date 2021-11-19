#### clean ####
rm(list=ls())
cat("\014")
setwd('C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati')


#### IMPORTAZIONE DATI ####
data <- read.table('dati.csv', head=TRUE, sep=',', row.names = "id")
data <- data[,-c(5,9,21,22,23,24,25,26)]
data[,2] <- ifelse(data[,2] == 1, "M", "F")
data[,3] <- ifelse(data[,3] == 1, "TR", "VS")
data[,4] <- ifelse(data[,4] == 1, "DX", "SX")
data[,7] <- ifelse(data[,7] == 'FRMN', "SI", "NO")
data[,11] <- ifelse(data[,11] == 'B', 2, 1)
data[,8] <- as.factor(toupper(data[,8]))
data[,9] <- as.factor(toupper(data[,9]))
data[,12] <- as.factor(toupper(data[,12]))
data[,13] <- as.factor(toupper(data[,13]))
data[data[,12]=='NO ',12] <- c('NO')
data[data[,13]=='NO ',13] <- c('NO')
data[,12] <- factor(data[,12])
data[,13] <- factor(data[,13])  

colnames(data) <- c("AGE", "SEX", "EZIOLOGIA", "LATE", "MESIC", "MESIR", "FRMN",
                    "SLDX", "SLSX", "SLSEP", "SLSEPINC", "MLDX", "MLSX", 
                    "MLSEP", "MLSEP2", "GOSE", "LCF", "DRS")
attach(data)
options("scipen"=100, "digits"=4)

#View(data)

setwd('C:/Users/pietr/Desktop/Bayesian Statistics/Progetto/dati/potenziali_evocati')
res <- list()
latency <- c("SL", "ML")
late <- c("sx", "dx")
w <- 1
for(i in 1:26)
{
  for(j in 1:2)
  {
    for(k in 1:2)
    {
      fileimp <- paste(paste("(", i, ")/", i, sep=''), latency[j], 
                    paste(late[k], ".asc", sep=''), sep = ' ')
      
      restemp <- try(read.table(file=fileimp, skip = 13, dec= ','))
      if('try-error' %in% class(restemp)) next
      else
      {
        res[[w]] <- c(restemp[1:1600,])
        res[[w]]$sogg <- c(i)
        res[[w]]$latency <- latency[j]
        res[[w]]$late <- late[k]
        w <- w+1
      }   
    }
  }
}

#### GENERO MATRICI PER OGNI CONDIZIONE TESTATA AURICOLARE ####


importMatrix <- function(x, type, position)
{
  temp <- matrix(0, nrow=length(res[[1]][[1]]))
  ntmp <- c()
  j <- 1
  for(i in 1:length(x))
  {
    if(res[[i]]$latency == type[1] && res[[i]]$late == type[2])
    {
      if(position == 'au')
      {
        temp <- cbind(temp, as.vector(res[[i]][[1]]))
        ntmp <- c(ntmp, res[[i]]$sogg)
      }
      else
      {
        if(length(res[[i]]) == 11)
        {
          temp <- cbind(temp, as.vector(res[[i]][[3]]))
          ntmp <- c(ntmp, res[[i]]$sogg)
        }
        else
        {
          if(is.numeric(res[[i]][[4]]))
          {
            temp <- cbind(temp, as.vector(res[[i]][[4]]))
            ntmp <- c(ntmp, res[[i]]$sogg)
          }
        }
      }
    }
  }
  colnames(temp) <- as.character(c(0,ntmp))
  temp[,-1]
}


library(fda.usc)

f.data <- list()
f.data$ausxSL <- fdata(t(importMatrix(res, type = c('SL', 'sx'), position = 'au')))
f.data$audxSL <- fdata(t(importMatrix(res, type = c('SL', 'dx'), position = 'au')))
f.data$frsxSL <- fdata(t(importMatrix(res, type = c('SL', 'sx'), position = 'fr')))
f.data$frdxSL <- fdata(t(importMatrix(res, type = c('SL', 'dx'), position = 'fr')))
f.data$ausxML <- fdata(t(importMatrix(res, type = c('ML', 'sx'), position = 'au')))
f.data$audxML <- fdata(t(importMatrix(res, type = c('ML', 'dx'), position = 'au')))
f.data$frsxML <- fdata(t(importMatrix(res, type = c('ML', 'sx'), position = 'fr')))
f.data$frdxML <- fdata(t(importMatrix(res, type = c('ML', 'dx'), position = 'fr')))

plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250))
for(i in 1:26){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[i,], col = i, lwd=2)
}

#### PLOT CARINI ####

x11()
plot(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[1,], type = "l", ylim = c(-250, 250), xlab = 'time [ ms ]', ylab = expression(paste("Evoked Potential [ " * mu ~ "V ]")), lwd =3)
index = c(21,2,4,22)
cols = c("forestgreen","firebrick2","royalblue","darkorchid")
for(i in 1:4){
  lines(t(importMatrix(res, type = c('SL', 'sx'), position = 'au'))[index[i],], col = cols[i], lwd=3)
}
