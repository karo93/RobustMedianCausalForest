myhodges_lehmann_estimator<- function(b,p)
{  sample_1 <- length(b)
   sample_2 <- length(p)
   k<-1
  if(sample_1==0){
    HL= median(p)
  } else if(sample_2==0){
    HL= median(b)
  } else{
    samplesize <- sample_1 * sample_2
    differences<- numeric(samplesize)
    for(i in 1:sample_1){
      for(j in 1:sample_2){
        differences[k] = b[i] - p[j]
        k<-k+1
      }
    }
    HL= median(differences)
  }
  return(HL)
 }
