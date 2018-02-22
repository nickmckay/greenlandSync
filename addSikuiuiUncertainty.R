library(geoChronR)
library(lipdR)
#add appropriate uncertainty into data..
NU= readLipd(here("data","Sikuiui.Thomas.2018.lpd"))


#Sikuiui first
#use select.data() to get the ageEnsemble and paleoData from each of the files
NU = mapAgeEnsembleToPaleoData(NU)
NU.ae = selectData(NU,varName = "ageEnsemble")
NU.depth = selectData(NU,varName = "depth")
NU.MAAT = selectData(NU,varName = "Sun2011MAAT")
NU.epsilon = selectData(NU, varName = "EpsilonC29-C23")


#Let's throw some uncertainty into MAAT

# ##UNCORRELATED NOISE!!!
# tempUnc = 2.5 #One sigma uncertainty
# nens = 1000
# NU.MAATEnsemble = NU.MAAT
# NU.MAATEnsemble$variableName = "NU.MAATEnsemble"
# NU.MAATEnsemble$values = matrix(NU.MAATEnsemble$values,ncol = nens,nrow = length(NU.MAAT$values))
# #lets add random normal uncertainty (uncorrelated)
# NU.MAATEnsemble$values = NU.MAATEnsemble$values + matrix(rnorm(1000,sd = tempUnc),ncol = nens,nrow = length(NU.MAAT$values))
# 
# ### END UNCORRELATED NOISE ###


##CORRELATED NOISE
tempUnc = 2.5 #One sigma uncertainty
nens = 1000
NU.MAATEnsemble = NU.MAAT
NU.MAATEnsemble$variableName = "NU.MAATEnsemble"
NU.MAATEnsemble$values = matrix(NU.MAATEnsemble$values,ncol = nens,nrow = length(NU.MAAT$values))
randos = matrix(NA,nens,nrow = length(NU.MAAT$values))

for(ne in 1:nens){
  randos[,ne] = arima.sim(list(order = c(1,0,0), ar = 0.7), n = length(NU.MAAT$values))
}

#detrend? Remove a linear trend from the uncertainty values
trendLines = predict(lm(randos ~ seq_len(length(NU.MAAT$values))),interval = prediction)
randos = randos - trendLines

#adjust for sd
randos = randos/sd(randos) * tempUnc

#take a look at the uncertainty and assess reasonableness
barplot(randos[,50])

#Add in corerelated noise
NU.MAATEnsemble$values = NU.MAATEnsemble$values + randos

### END UNCORRELATED NOISE ###

#add into dataset
NU$paleoData[[1]]$model[[1]] <- list()

NU$paleoData[[1]]$model[[1]]$ensembleTable[[1]] <- list()
NU$paleoData[[1]]$model[[1]]$ensembleTable[[1]]$MAAT <- NU.MAATEnsemble
NU$paleoData[[1]]$model[[1]]$ensembleTable[[1]]$age <- NU$paleoData[[1]]$measurementTable[[1]]$age
NU$paleoData[[1]]$model[[1]]$ensembleTable[[1]]$ageEnsemble <- NU$paleoData[[1]]$measurementTable[[1]]$ageEnsemble

NU$paleoData[[1]]$model[[1]]$method <-list()
NU$paleoData[[1]]$model[[1]]$method$description <- "Added autocorrelated noise to reconstructed temperature to mimic bias and uncertainty"
NU$paleoData[[1]]$model[[1]]$method$code <- "https://github.com/nickmckay/greenlandSync/blob/master/addSikuiuiUncertainty.R"

writeLipd(NU,path = here("data"))



