# FONCTION DU MODELE 2
# Algorithme Metropolis Hastings
MetropolisHastings_LogMod <- function(proposal,burning,iteration,Y,X){
  burn <- matrix(rep(0,length(proposal)*burning), ncol=length(proposal))
  iter <- matrix(rep(0,length(proposal)*iteration), ncol=length(proposal))
  
  burn[1,]<-proposal
  # BURNING
  for(i in 2:burning){
    betaT <- rmvnorm(1,proposal,diag(rep(1,4)))
    rho <- min(1,RatioPosterior_LogMod(Y,betaT,burn[i-1,],X)*dmvnorm(burn[i-1,],betaT,diag(rep(1,4)))/dmvnorm(betaT,burn[i-1,],diag(rep(1,4))))
    if(rbinom(1,1,rho)){
      burn[i,]<-betaT
    }else{
      burn[i,]<-burn[i-1,]
    }
    print(i)
  }
  
  # METROPOLIS HASTINGS
  iter[1,]<-burn[burning,]
  
  for(i in 2:iteration){
    betaT<-rmvnorm(1,proposal,diag(rep(1,4)))
    rho<-min(1,RatioPosterior_LogMod(Y,betaT,iter[i-1,],X)*dmvnorm(iter[i-1,],betaT,diag(rep(1,4)))/dmvnorm(betaT,iter[i-1,],diag(rep(1,4))))
    if(rbinom(1,1,rho)){
      iter[i,]<-betaT
    }else{
      iter[i,]<-iter[i-1,]
    }
    print(i)
  }
  return(list(burningMH=burn,iterationMH=iter))
}