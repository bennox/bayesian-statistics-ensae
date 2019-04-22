## FONCTION DU MODELE 2
# Transformation des donnees en matrices pour le modele
MatrixData<-function(distance,adjacence0,matF){
  n_pop = nrow(adjacence0)
  res<-matrix(rep(0,n_pop*(n_pop-1)/2*4),ncol=4)
  k<-1
  for(i in 1:(n_pop-1)){
    for(j in (i+1):n_pop){
      res[k,1]<-1
      res[k,2]<-distance[i,j]
      res[k,3]<-adjacence0[i,j]
      res[k,4]<-matF[i,j]
      k<-k+1
    }
  }
  return(res)
}