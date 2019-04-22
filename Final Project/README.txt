# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ENSAE - 3A - Statistique bayésienne
# Janvier 2017
# Tom Duchemin, Mehdi Miah, Benoit Robaglia
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

## Exécution
1. lancer main.R
	a. Premièrement et nécessairement, lancer la section "Simulation des données"
	b. Ensuite, lancer la section relative au modèle d'intérêt : modèle 1, 2 ou 3.

## Organisation des fichiers
Le fichier main.R comprend la simulation des données et les procédures d'estimation.
Les fonctions utilisées dans ce main.R sont entrés dans des codes séparés dans le dossier ./estimation (le fichier principal contient des commandes permettant d'importer ces fonctions).

## Présentation des fonctions du dossier ./estimation
Les fonctions du dossier ./estimation sont les suivantes :
# Fonctions du modèle 1
- post_Mod1.R : fonction qui renvoie la loi a posteriori des paramètres du modèle
- prior_beta_s2_Mod1.R : fonction qui renvoie la loi a priori des paramètres
- log_l.R : fonction qui renvoie le logarithme du produit de la posterior et des prior des paramètres d'intérêt (cette fonction a pour objectif d'être utilisée dans un algorithme de Metropolis-Hastings)

# Fonctions du modèle 2
- VectAdj.R : fonction qui vectorise une matrice d'adjacence
- MatrixData.R : fonction qui met les observations dans une matrix cohérente avec la matrice d'adjacence vectorisée
- NormPosterior_LogMod.R : fonction qui renvoie le ratio des fonctions proposal gaussiennes des paramètres de l'itération MH
- LogisticPosterior_LogMod.R : fonction qui renvoie le ratio des posteriors pour les paramètres de l'itération MH
- RatioPosterior_LogMod.R : fonction qui calcul le ratio des deux fonctions précédentes (pour calculer la probabilité d'acceptation de l'itération MH)
- MetropolisHastings_LogMod.R : algorithme de Metropolis Hastings pour le modèle 2

