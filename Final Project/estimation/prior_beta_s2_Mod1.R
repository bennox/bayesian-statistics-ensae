# Modele 1
#Loi a priori des parametres
prior_beta_s2_Mod1 = function (param){
  return(dmvnorm(param[1:4],mu0,sigma2*solve(S0))*dinvchisq(param[5],4))
}