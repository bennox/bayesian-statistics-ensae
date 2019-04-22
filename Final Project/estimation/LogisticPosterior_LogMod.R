# FONCTION DU MODELE 2
# Calcul du rapport des posteriors pour introduction dans l'algorithme Metroplis-Hastings
LogisticPosterior_LogMod<-function(Y, beta1, beta2, X){
  obj1<-exp(X%*%t(beta1))/(1+exp(X%*%t(beta1)))
  res11<-obj1^Y
  res12<-(1-obj1)^(1-Y)
  obj2<-exp(X%*%beta2)/(1+exp(X%*%beta2))
  res21<-obj2^Y
  res22<-(1-obj2)^(1-Y)
  
  return(prod(res12*res11/(res21*res22)))
}