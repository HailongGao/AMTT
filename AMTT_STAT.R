##################################Calculate the value of the test statistic##################

midat<-function(X1,X2,t){
  # input: X1,X2@Sample
  #            t@Truncation parameter
  #
  #output: the value of the test statistic
  n1 = nrow(X1)
  n2 = nrow(X2)
  p<-ncol(X1)
  s1bar<-colMeans(X1)
  s2bar<-colMeans(X2)
  sigma<-((apply(X1, 2, var))*(n1-1)+(apply(X2, 2, var))*(n2-1))/(n1+n2)
  if(min(sigma)==0){
    stop("The individuals of two samples have equal observation for at least one OTU!")
  }
  Sta_i<-(((s1bar-s2bar)/sqrt(sigma))^2)
  Sta_sort<-sort(abs(Sta_i),decreasing = TRUE)
  sorted_index <- order(abs(Sta_i),decreasing = TRUE)
  sample1_rank<-X1[,sorted_index]
  sample2_rank<-X2[,sorted_index]
  mean1<-colMeans(sample1_rank)
  mean2<-colMeans(sample2_rank)
  sigma_rank<-((apply(sample1_rank, 2, var))*(n1-1)+(apply(sample2_rank, 2, var))*(n2-1))/(n1+n2)
  results<-rep(NA,length(t))
  for (i in 1:length(t)) {
    index<-p*t[i]
    index <- 1:index
    results[i] <- ((n1*n2)/(n1+n2))*(sum(((mean1[index] - mean2[index]) / sqrt(sigma_rank[index]))^2))
    
  }
  
  return(results=results)
}


########################calculate P value by parameterization method###################################

amidat_parameterization <- function(X1,X2,B,t) {
  # input: X1,X2@Sample
  #            B@The number of permutated datasets     
  #            t@Truncation parameter
  #
  #output: Combined P-values of  amidat
  
  n1 = nrow(X1)
  n2 = nrow(X2)
  p<-ncol(X1)
  Y<-rbind(X1,X2)
  results<-midat(X1,X2,t)
  
  t_h0_B = matrix(NA, nrow = B, ncol = length(t)) 
  for (b in 1:B){
    idx = sample.int(n1+n2, n1) 
    sample1_rank_b = Y[idx,] 
    sample2_rank_b = Y[-idx,] 
    t_h0_B[b,] = midat(sample1_rank_b,sample2_rank_b,t)
  }
  
  kafang_list <- list()
  abd_list <- list()
  
  for (i in 1:length(t)){
    
    K1<-sum(t_h0_B[,i])/B
    K2<-sum((t_h0_B[,i]-K1)^2)/B
    K3<-sum((t_h0_B[,i]-K1)^3)/B
    a<-K3/(4*K2)
    b<-K1-(2*K2^2)/K3
    d<-(8*K2^3)/K3^2
    kafang<-a*(rchisq(B, d))+b
    kafang_list[[i]] <- kafang
    abd_list[[i]] <- c(a, b, d)
  }
  Pi<-rep(NA,length(t))
  for(i in 1:length(t)){
  abd<-abd_list[[i]]
  P_i<-1-pchisq((results[i]-abd[2])/abd[1], abd[3])
  Pi[i]<-P_i}
  cauchy_statistic<-(sum(tan((0.5-Pi)*pi)))/length(t)
  P_final<-0.5-(atan(cauchy_statistic)/pi)
  return(list(results = results, t_h0_B = t_h0_B,kafang_list=kafang_list,Pi=Pi,P_final=P_final,abd_list=abd_list))
}

