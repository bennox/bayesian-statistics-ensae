# FONCTION MODELE 2
# Ratio des posteriors des parametres 
NormPosterior_LogMod<-function(beta1,beta2){
  return(prod(exp(-1/2*beta1^2)/exp(-1/2*beta2^2)))
}