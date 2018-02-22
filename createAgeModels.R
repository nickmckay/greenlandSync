#make Bacon/BAM age models
#Load relevant libraries
library(lipdR)
library(geoChronR)
library(here)
library(ggplot2)
library(tidyverse)
#Load in data
D = readLipd(here("data"))

#Loop through looking for models
for(i in 1:length(D)){
  L = D[[i]]
  if(is.null(L$chronData[[1]]$model)){#then create model...
    print(paste("Creating model for",L$dataSetName,"..."))
    
    #if there's a measurement table - BACON!
    if(!is.null(L$chronData[[1]]$measurementTable)  & grepl("sed",tolower(L$archiveType))){
      L = runBacon(L) #generate the bacon model
      
      #plot the model
      agePlot = plotChron(L)
      
      #add a layer that plots age from paleoData
      depth = selectData(L,varName = "depth")
      age = selectData(L,varName = "age")
      
      pd = tibble(age = age$values, depth = depth$values)
      pd <- filter(pd,!is.na(age))
      
      
      
      agePlot = agePlot + geom_line(data = pd, aes(y = depth, x = age,colour = "Original")) + 
        scale_color_manual("Model",values = c("green")) + xlab("Age (yr BP)") + scale_y_reverse("Depth (cm)")
      
      print(agePlot)
      ggsave(agePlot, filename = here("figures",paste0(L$dataSetName,"-newAgeModel.pdf")))
      
      answer = gtools::ask("Does this look good?")
    
      if(grepl("y",tolower(answer))){#save the output and file.
      writeLipd(L,path = here("data"))
      }else{
        stop("failed")
      }
    
    }else{#otherwise BAM!
      L=runBam(L,model = list(param = 0.02,ns = 1000,name = "poisson",resize = 0))
      writeLipd(L,path = here("data"))
    }
    
    
  }
  
}

