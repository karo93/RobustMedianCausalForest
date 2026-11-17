my_HL <- function(y,wt,treatment){
  n <- length(y) 
  tmpsum <- sum(wt)
  splits <- numeric(n-1)
  for(i in 1:(n-1)) {
    tmp_1 <- 0
    y_i<-y[1:i]
    treat_i<-treatment[1:i]
    m_11 <-c(y_i[treat_i==1])
    m_10 <-c(y_i[treat_i==0])
    m_1 <- myhodges_lehmann_estimator(m_11,m_10)
    print(m_1)
    for(j in 1:i) {
      tmp_1<- tmp_1+sum(wt[j] * (abs(y[j] - m_1)))
    }
    tmp_2<-0
    y_j<-y[(i+1):n]
    treat_j<-treatment[(i+1):n]
    m_21 <-c(y_j[treat_j==1])
    m_20 <-c(y_j[treat_j==0])
    m_2 <- myhodges_lehmann_estimator(m_21,m_20)
    for(j in (i+1):n) {
      tmp_2 <- tmp_2+sum(wt[j] * (abs(y[j] - m_2)))
    }
    splits[i] <- (tmp_1 +tmp_2) / tmpsum
  }
  return(splits)
}

#for(i= 0;i < sample_1; i++){
 # for(j=0;j < sample_2;j++){
  #  differences[k] = b[i] - p[j];
   # k++;
#  }
#}
#HL= quick_select5(differences,samplesize);
