# MODELE 2
#  Ratio des lois a posteriori pour le Metropolis-Hastings
RatioPosterior_LogMod<-function(Y,beta1,beta2,X){
  return(LogisticPosterior_LogMod(Y,beta1,beta2,X)*NormPosterior_LogMod(beta1,beta2 ))
}