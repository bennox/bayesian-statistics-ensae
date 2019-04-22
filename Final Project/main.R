# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ENSAE - 3A - Statistique bayésienne
# Janvier 2017
# Tom Duchemin, Mehdi Miah, Benoit Robaglia
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# == Preambule ==============================================================

# Suppression de toutes les variables
rm(list=ls())

# Chargement des librairies
library(truncnorm)
library(igraph)
library(sna)
library(mvtnorm)
library(MHadaptive)
library(mvtnorm)
library(MCMCpack)
library(geoR)


# == Simulation a la periode t=1 ============================================

# Simulation du reseau a t=1 (matrice D0)
## la simulation est faite a partir de l'algorithme de Watts and Strogatz 
## - les parametres dim et size sont choisis tels que le nombre d'individus soit 
## proches de celui de l'article ; 
## - les parametres nei et p sont choisis tels que les statistiques du reseau
## soient proches de celles de l'article

D0 <- sample_smallworld(dim = 2, size = 23, nei = 2, p = 0.01)
adjD0 <- as_adj(D0, sparse=FALSE) # matrice d'adjacence

# Dans la simulation, le nombre d'individu n'est pas choisi par l'utilisateur
n_pop = nrow(adjD0) #nombre d'individus
cat(sprintf("La base simulée contient %.0f individus.", n_pop))

# Simulation des notes à la t=1 à partir d'une normale tronquée
X0 <- round(rtruncnorm(n_pop, a = 0, b = 4, mean = 2.6, sd = 0.8), 1)

# Distance entre chaque individus en termes de notes
distX <- matrix(sapply(X = seq(1,n_pop), Y = seq(1,n_pop),
                       function(X,Y) abs(X0[X]-X0[Y])), 
                ncol = n_pop)

# Nombre d'amis en commun à t=1
F0 <- matrix(rep(0,n_pop*n_pop),nrow = n_pop)

for(i in 1:n_pop){
  for(j in i:n_pop){
    res <- adjD0[i,]+adjD0[j,] 
    
    #res est un vecteur de taille n_pop contenant des 0, des 1 et des 2
    #un 2 signifie qu'il existe un individu k tel que i et j sont amis avec k (k peut être i ou j !)
    
    F0[i,j] <- ifelse(sum(res==2) - adjD0[i,j] >= 1, 1, 0) # nombre d'amis en commun 
  }
}

F0<-F0+t(F0)

# == Simulation à la période t=2 ============================================

# Parametres de l'article
alpha_0 = -2.56
alpha_x = -0.2
alpha_d = 2.52
alpha_f = 1.2

# Fonction d'utilite
U <- alpha_0 + alpha_x*distX + alpha_d*adjD0 + alpha_f*F0

# Probabilite de creer un lien sachant l'utilite
plien <- (exp(U)/(1+exp(U)))^2 

# Simulation du reseau a partir de la fonction d'utilite
adjD1 <- rgraph(n_pop, m = 1, tprob = plien, mode = "graph", diag = FALSE) 
D1 <- graph_from_adjacency_matrix(adjD1)

# Matrices utiles pour la generation des notes en t=2 (Y)
M <- sapply(X = seq(1,n_pop), function(X) sum(adjD1[X,]))

G <- matrix(rep(0,n_pop*n_pop), ncol = n_pop) # G : "matrice D row-normalise"
for(i in 1:n_pop){
  for(j in 1:n_pop){
    G[i,j]<-ifelse(M[i] == 0, 0, adjD1[i,j]/M[i])
  }
}

# Parametres de l'article
beta_0 = -0.13
beta_x = 0.74
beta_y_bar = 0.16
beta_x_bar = 0.11
sigma2 = 0.61

# Constuction des notes a la seconde periode
X1 <- solve(diag(n_pop)-beta_y_bar*G) %*% matrix(rep(beta_0,n_pop), ncol = 1) + 
  solve(diag(n_pop)-beta_y_bar*G) %*% (beta_x * X0 + beta_x_bar*G %*% X0) + 
  rnorm(n_pop,0,sigma2)


# == Statistiques sur les variables simulees ================================

# Affichage du réseau d'amitie a la periode t=1
plot(D0, edge.arrow.size = .01, vertex.label = NA, vertex.size = 5, edge.width = .01)

# Affichage du reseau d'amitie a la periode t=2
plot(D1, edge.arrow.size = .01, vertex.label = NA, vertex.size = 5)

# Nombre moyen d'amis à chaque période
sum(adjD0/(2*(n_pop-1))) # en t=1
sum(adjD1/(2*(n_pop-1))) # en t=2

# Dynamique de réseau
dynam<-matrix(rep(0,4),nrow=2)
colnames(dynam)<-c("Pas amis en 0","Amis en 0")
rownames(dynam)<-c("Pas amis en 1","Amis en 1")
adj<-adjD0+adjD1
for(i in 1:n_pop){
  dynam[1,1]<-length(adj[i,adj[i,]==0])+dynam[1,1]-1
  dynam[2,2]<-length(adj[i,adj[i,]==2])+dynam[2,2]
  dynam[1,2]<-length(adj[i,adj[i,]==1 & adjD0[i,]==1])+dynam[1,2]
  dynam[2,1]<-length(adj[i,adj[i,]==1 & adjD1[i,]==1])+dynam[2,1]
}
dynam<-dynam/2
print(dynam)

# Notes
mean(X0) # moyenne des notes en 1
mean(X1) # moyenne des notes en 2

############################################################################################
# == Modele 1 - Estimation par le modele lineaire ===========================================

# Matrice des donnees
tn = rep(1,n_pop)
V = matrix(c(tn,X0,G%*%X1,G%*%X0), nrow = n_pop, ncol =4)

# Valeur pour les priors des beta (on mettra 1 pour sigma2)
mu0 = c(-0.13,0.74,0.16,0.11)
S0 = diag(c(0.12,0.04,0.05,0.07))

# Fonctions necessaires pour l'algorithme de Metropolis-Hastings
source("./estimation/postMod1.R") # posterior du modele
source("./estimation/prior_beta_s2_Mod1.R") # prior des parametres
source("./estimation/log_l.R") # log du produit des deux fonctions precedentes


AlgoMod1<-Metro_Hastings(log_l,pars = c(mu0,1), iterations = 100000 , burn_in = 10000)
#Algorithme de metropolis hastings fourni par le package MHadaptive. Nous avons choisi un algorithme deja
#optimise pour pouvoir effectuer un grand nombre d'iterations. Un algorithme M-H a été codé pour les 2 modeles
#suivants.


## Convergence
par(mfrow=c(2,2))
plot(AlgoMod1$trace[,1],type="l",main="Convergence de beta_0")
plot(AlgoMod1$trace[,2],type="l",main="Convergence de beta_x")
plot(AlgoMod1$trace[,3],type="l",main="Convergence de beta_ybar")
plot(AlgoMod1$trace[,4],type="l",main="Convergence de beta_xbar")

# Moyenne et Ecart-type des estimations
# beta_ 0
mean(AlgoMod1$trace[,1]) 
sd(AlgoMod1$trace[,1])
# beta_x
mean(AlgoMod1$trace[,2])
sd(AlgoMod1$trace[,2])
# beta_ybar
mean(AlgoMod1$trace[,3])
sd(AlgoMod1$trace[,3])
# beta_xbar
mean(AlgoMod1$trace[,4])
sd(AlgoMod1$trace[,4])
# sigma2
mean(AlgoMod1$trace[,5])
sd(AlgoMod1$trace[,5])


############################################################################################
# == Modele 1  Estimation par le modele exogene ===========================================

source("./estimation/VectAdj.R")
source("./estimation/MatrixData.R")
source("./estimation/LogisticPosterior_LogMod.R")

# Fonctions nécessaires pour l'algorithme de Metropolis-Hastings
source("./estimation/NormPosterior_LogMod.R")
source("./estimation/RatioPosterior_LogMod.R")

# Algorithme de Metropolis Hastings
source("./estimation/MetropolisHastings_LogMod.R")

AlgoMod2 <- MetropolisHastings_LogMod(c(alpha_0, alpha_x, alpha_d, alpha_f),
                                  100, 5000, VectAdj(adjD1), MatrixData(distX,adjD0,F0))


# Convergence des parametres
mean1<-rep(AlgoMod2$iterationMH[1,1],length(Algo$iterationMH[,1]))
mean2<-rep(AlgoMod2$iterationMH[1,2],length(Algo$iterationMH[,2]))
mean3<-rep(AlgoMod2$iterationMH[1,3],length(Algo$iterationMH[,3]))
mean4<-rep(AlgoMod2$iterationMH[1,4],length(Algo$iterationMH[,4]))
for(i in 2:length(AlgoMod2$iterationMH[,1])){
  mean1[i]<-mean(AlgoMod2$iterationMH[1:i,1])
  mean2[i]<-mean(AlgoMod2$iterationMH[1:i,2])
  mean3[i]<-mean(AlgoMod2$iterationMH[1:i,3])
  mean4[i]<-mean(AlgoMod2$iterationMH[1:i,4])
}

plot(mean1,type="l",xlab="iterations",ylab="alpha_0",main="Convergence de alpha_0")
plot(mean2,type="l",xlab="iterations",ylab="alpha_x",main="Convergence de alpha_x")
plot(mean3,type="l",xlab="iterations",ylab="alpha_f",main="Convergence de alpha_f")
plot(mean4,type="l",xlab="iterations",ylab="alpha_d",main="Convergence de alpha_d")


# == Estimation par le modele endogene ======================================
