# FONCTION DU MODELE 1 : 
## Log du produit de la posterior et des prior
log_l = function(param){ # param : vecteur des 5 parametres
  return(log(postMod1(param)*prior_beta_s2_Mod1(param)))
}
#log_l(mu0)
