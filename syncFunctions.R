

#Define function to test for probability of age uncertain event synchroneity:
#Let's calculate the probability that two events are simultaneous within some threshold
ensCoevalProb = function(e1,e2,thresh=10,maxEns=1000){
  #input is two vectors of samples from a distribution
  e1 = c(na.omit(e1))#force to vector
  e2 = c(na.omit(e2))#force to vector
  if(is.na(maxEns)){
    maxEns1 = length(e1)
    maxEns2 = length(e2)
  }else{
    maxEns1 = min(c(maxEns,length(e1)))
    maxEns2 = min(c(maxEns,length(e2)))
  }
  #randomly select maxEns elements
  e1 = e1[sample.int(length(e1),size=maxEns1)]
  e2 = e2[sample.int(length(e2),size=maxEns2)]
  
  
  #force the two samples to be the same size be resampling the smaller 1.
  mediff = maxEns1-maxEns2
  
  if(mediff<0){
    #then e1 is shorter
    e1 = c(e1,e1[sample.int(length(e1),size=abs(mediff))])
  }else if(mediff>0){
    #then e2 is shorter
    e2 = c(e2,e2[sample.int(length(e2),size=abs(mediff))])
  }
  
  ediff = e1-e2
  outcome = matrix(0,nrow = length(ediff))
  for(e in 1:length(ediff)){
    if(abs(ediff[e])<=thresh){
      outcome[e]=1
    }
  }
  
  prob = mean(outcome)
  return(prob)
}



#