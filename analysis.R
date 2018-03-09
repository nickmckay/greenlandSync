#Load relevant libraries
library(lipdR)
library(geoChronR)
library(here)
library(changepoint)
library(tidyverse)
#and functions
source(here("syncFunctions.R"))


#Load in data
D = readLipd(here("data"))

#map age ensemble as needed
for(i in 1:length(D)){
  if(is.null(D[[i]]$paleoData[[1]]$measurementTable[[1]]$ageEnsemble)){
    #sediment?
    if(grepl("sed",D[[i]]$archiveType,ignore.case = T)){
      D[[i]] = mapAgeEnsembleToPaleoData(D[[i]],age.var = "ageEnsemble")
    }else{#ice core
      
      if(!is.null(D[[i]]$chronData[[1]]$model[[1]]$ensembleTable[[1]]$ageEnsemble)){
        D[[i]]$paleoData[[1]]$measurementTable[[1]]$ageEnsemble <- D[[i]]$chronData[[1]]$model[[1]]$ensembleTable[[1]]$ageEnsemble
      }else if(!is.null(D[[i]]$chronData[[1]]$model[[1]]$ensembleTable[[1]]$yearEnsemble)){
        ye = D[[i]]$chronData[[1]]$model[[1]]$ensembleTable[[1]]$yearEnsemble
        ye$variableName <- "ageEnsemble"
        ye$units <- "BP"
        ye$values <- 1950- ye$values
        D[[i]]$paleoData[[1]]$measurementTable[[1]]$ageEnsemble <- ye
      }else{
        print(paste("no ensemble data for",D[[i]]$dataSetName))
      }
    }
  }
}


TS = extractTs(D,whichtables = "all")


varNames = sapply(TS,"[[","paleoData_variableName")
sapply(TS,"[[","tableType")


vars2Analyze = c("Temp","temp","accumulation","d18O","seaIceCover","EpsilonC29-C23","Sun2011MAAT","C23d2H","C25d2H","C27d2H","C29d2H","AccRate","TTL CARB")
sTS = TS[varNames %in% vars2Analyze]
dsn = sapply(sTS,"[[","dataSetName")
#Data ready to analyze!

# Sea ice records:
#   CC04.Gibb.2015: sea ice in the labrador sea. 14C age model. Standard marine reservoir correction (Delta R is zero)
# CC70.Gibb.2015: sea ice in southern Baffin Bay. 14C age model. Standard marine reservoir correction (Delta R is zero)
# LS009.Ledu.2010: sea ice in Lancaster Sound.  14C age model. Delta R is 400. 14C only covers the early Holocene (so we’re fine for this study). I couldn’t figure out how to add extra “calibrated age” columns, which will eventually include the paleomagnetic tie points. I didn’t bother since these don’t affect the Early Holocene. Let me know if you’d like that info.
# It’d be great to compare each of the sea ice records with our aquatic leaf wax d2H, epsilon, and temperature records.
# 
# Ice core records:
#   GISP2.Kobashi.2017: temp, age. Maximum counting error for the Early Holocene is 2% (Vintner et al. 2006, JGR)
# GISP2.Stuiver.1995: d18O, age, depth. Would be nice to figure out how to put this on the same age model as Kobashi, although Kobashi doesn’t give us depth, so I’m not totally sure how to approach that…maybe there’s a version of this record in Iso2k on the GICC05 age model, which is what Kobashi used?
# NEEM.Rasmussen.2013: accumulation rate, age & depth (plus accumulation uncertainty)
# It’d be nice to compare the GISP2 temperature data to our temperature record.
# Also would be nice to compare NEEM accumulation to our epsilon record
# All three of these ice core records are on (or in the case of Stuiver 1995 should be on) the GICC05 scale. For this age model, maximum counting error for the Early Holocene is 2% (Vinther et al. 2006, JGR)

Q <- c(8,3,3,3,5,3,3,3,3,4,3,3,3,3,3, 3,3)
                   #meese.          #andersen are problematic. Andersen is too short and Meese is missing the ensemble data. 

pen = c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.01,  .1,.1,.1,.1,     .1,.1)
#cpt testing.

time.range <- c(7000,9500)
names(D)[8]

for(i in c(1:5,7,9:17)){
y = sTS[[i]]$paleoData_values
if(grepl("epsilon",sTS[[i]]$paleoData_variableName,ignore.case = T)){
 y[y < -20] <- NA
}
ae = sTS[[i]]$ageEnsemble
med.age <- apply(ae,1,median)
good = which(!is.na(y) & between(med.age,min(time.range),max(time.range)))
y <- y[good]
ae <- ae[good,]

plot(ae[,1],y,main = dsn[i],type="l")
cp <- cpt.mean(y,method = "SegNeigh",penalty="Manual",pen.value = pen[i],Q=Q[i])
plot(cp)
i
#get change point indices..
cps <- cp@cpts
if(length(cps)==3){#found two changepoints.
  cpts <- cps[1:2]
}else if(length(cps)==4){
  cpts <- cps[1:2]
}else if(length(cps)==5){
  cpts <- cps[1:4]
}else if(length(cps)==2){
  cpts <- c(cps[1],cps[1])
}else if(length(cps)==8){
  cpts <- cps[2:7]
}

all.cp = matrix(NA,nrow = ncol(ae), ncol = length(cpts) )
all.cp.med = c()
#average rows
for(nc in 1:length(cpts)){
all.cp[,nc] <- colMeans(ae[c(cpts[nc],cpts[nc]+1),])

all.cp.med[nc] = median(all.cp[,nc])
}

tsPlot <- plotTimeseriesEnsRibbons(ae,y,colorLow = "LightGray",colorHigh = "DarkGray") +
  scale_x_reverse(name = "Age (yr BP)",limits = time.range[order(-time.range)] ) + 
  ylab(paste0(sTS[[i]]$paleoData_variableName," (",sTS[[i]]$paleoData_units,")")) +
 
  ggtitle(paste0(dsn[i],"-",sTS[[i]]$paleoData_variableName))

#add pairs of lines
os = seq(1,length(cpts)-1,by=2)
tm = seq(2,length(cpts),by=2)

      tsPlot <-  tsPlot + geom_vline(aes(xintercept = all.cp.med[os]),color = "blue") + geom_vline(aes(xintercept = all.cp.med[tm]),color = "red")



      
      
cpHist <- ggplot()


for(nc in seq(1,length(cpts),by = 2)){
  
  cp1 <- data.frame(age = all.cp[,nc])
  cp2 <- data.frame(age = all.cp[,nc+1])
  
  
  cpHist <- cpHist+geom_density(data = cp1,aes(x = age, fill = "Termination"), alpha = 0.5, colour = NA) +
  scale_x_reverse(name = "Age (yr BP)",limits = time.range[order(-time.range)] )+
  theme_bw()+scale_fill_brewer("Changepoint",palette = "Set1")+
  ylab("Probability Density")


#if(cp1.med!=cp2.med){#if they're not the same
  cpHist <-  cpHist+geom_density(data = cp2, aes(x = age, fill = "Onset"), alpha = 0.5, colour = NA) +
    geom_vline(data = cp2, aes(xintercept = median(age)),color = "red")
#}
cpHist <-  cpHist+  
  geom_vline(data = cp1, aes(xintercept = median(age)),color = "blue")

if(i!=1){
  cpHist <-  cpHist+  theme(legend.position = c(0.1, 0.7))
}else{
  cpHist <-  cpHist+  theme(legend.position = c(0.8, 0.7))
}

}

p1 <- ggplot_gtable(ggplot_build(tsPlot))
p2 <- ggplot_gtable(ggplot_build(cpHist))

maxWidth = grid::unit.pmax(p1$widths[2:3], p2$widths[2:3])
p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
finalFig <- gridExtra::grid.arrange(p1, p2, heights = c(3, 2))

ggsave(finalFig,filename = here("figures",paste0(dsn[i],"-",sTS[[i]]$paleoData_variableName,"-Timeseries&changepoints.pdf")))

}



ensCoevalProb(cp1Age,cp2Age,thresh= 1000)



