# MODELE 2
# Vectorisation de la matrice d'adjacence
VectAdj<-function(Adjacence){
  n_pop = nrow(Adjacence)
  res<-matrix(rep(0,n_pop*(n_pop-1)/2), ncol=1)
  k<-1
  for(i in 1:(n_pop-1)){
    for(j in (i+1):n_pop){
      res[k]<-Adjacence[i,j]
      k<-k+1  
    }  
  }  
  return(res) 
}