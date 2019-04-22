# MODELE 1
# Loi a posteriori des parametres
postMod1 = function(param){
  # Partie gaussienne
  mu = solve(t(V)%*%V+solve(S0))%*%(t(V)%*%X1+S0%*%mu0) 
  S = sigma2*(t(V)%*%V+solve(S0))
  # Partie Inverse Gamma
  param1<-(4+n_pop)/2
  param2<-(1+t(X1)%*%X1+t(mu0)%*%(S0)%*%mu0-t(mu)%*%(t(V)%*%V+S0)%*%mu)/2
  return(dmvnorm(param[1:4],mu,S)*dinvgamma(param[5],param1,param2))
}
#post(c(mu0,1))