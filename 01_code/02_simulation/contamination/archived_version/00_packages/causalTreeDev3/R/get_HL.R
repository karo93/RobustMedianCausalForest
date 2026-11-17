get_HL <- function(y,wt,treatment){
  n <- length(y)
  
  getHL<-numeric(2*n-2)
  for(i in 1:(n-1)) {
    y_i<-y[1:i]
    treat_i<-treatment[1:i]
    m_11 <-c(y_i[treat_i==1])
    m_10 <-c(y_i[treat_i==0])
    getHL[i] <- myhodges_lehmann_estimator(m_11,m_10)
    y_j<-y[(i+1):n]
    treat_j<-treatment[(i+1):n]
    m_21 <-c(y_j[treat_j==1])
    m_20 <-c(y_j[treat_j==0])
    getHL[n-1+i] <- myhodges_lehmann_estimator(m_21,m_20)
  }
  return(getHL)
  #print(getHL)
}


#for(i= 0;i < sample_1; i++){
# for(j=0;j < sample_2;j++){
#  differences[k] = b[i] - p[j];
# k++;
#  }
#}
#HL= quick_select5(differences,samplesize);
