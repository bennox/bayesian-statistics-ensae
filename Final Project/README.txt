# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ENSAE - 3A - Statistique bay�sienne
# Janvier 2017
# Tom Duchemin, Mehdi Miah, Benoit Robaglia
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

## Ex�cution
1. lancer main.R
	a. Premi�rement et n�cessairement, lancer la section "Simulation des donn�es"
	b. Ensuite, lancer la section relative au mod�le d'int�r�t : mod�le 1, 2 ou 3.

## Organisation des fichiers
Le fichier main.R comprend la simulation des donn�es et les proc�dures d'estimation.
Les fonctions utilis�es dans ce main.R sont entr�s dans des codes s�par�s dans le dossier ./estimation (le fichier principal contient des commandes permettant d'importer ces fonctions).

## Pr�sentation des fonctions du dossier ./estimation
Les fonctions du dossier ./estimation sont les suivantes :
# Fonctions du mod�le 1
- post_Mod1.R : fonction qui renvoie la loi a posteriori des param�tres du mod�le
- prior_beta_s2_Mod1.R : fonction qui renvoie la loi a priori des param�tres
- log_l.R : fonction qui renvoie le logarithme du produit de la posterior et des prior des param�tres d'int�r�t (cette fonction a pour objectif d'�tre utilis�e dans un algorithme de Metropolis-Hastings)

# Fonctions du mod�le 2
- VectAdj.R : fonction qui vectorise une matrice d'adjacence
- MatrixData.R : fonction qui met les observations dans une matrix coh�rente avec la matrice d'adjacence vectoris�e
- NormPosterior_LogMod.R : fonction qui renvoie le ratio des fonctions proposal gaussiennes des param�tres de l'it�ration MH
- LogisticPosterior_LogMod.R : fonction qui renvoie le ratio des posteriors pour les param�tres de l'it�ration MH
- RatioPosterior_LogMod.R : fonction qui calcul le ratio des deux fonctions pr�c�dentes (pour calculer la probabilit� d'acceptation de l'it�ration MH)
- MetropolisHastings_LogMod.R : algorithme de Metropolis Hastings pour le mod�le 2

