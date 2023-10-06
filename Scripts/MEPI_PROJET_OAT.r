# Packages
# install.packages("sensibility")
library(sensitivity)
library(ggplot2)
library(latex2exp)
library(numDeriv)
library(tidyverse)
library(latex2exp)
library(viridis)


# Fonction de base 
modAppli <- function(parametre){  
  
  # CONDITIONS DE SIMULATION
  temps = 2*365;  # Nombre de pas de temps en jours
  
  # Initialisation pour la sauvegarde de 4 sorties d'indicateurs epidemiologiques pour chaque jeu de parametres
  sorties <-
    matrix(0,
           nrow = nrow(parametre),
           ncol = 4,
           dimnames = list(c(paste0("scneario_", 1:nrow(parametre))),
                           c("tx_morbidite", "incidence_t730", "pic_infectieux", "prevalence_annee_1")))
  
  
  
  # Initialisation pour la sauvegarde des matrices des effectifs par scenarios
  sortie.MAT <- list()
  
  # Initialisation pour la sauvegarde des incidences journalieres
  sortie.nvinf <- list()
  
  # Boucle des scenarios, autant que de jeux de donnees
  for (i in 1:nrow(parametre)) { 
    
    # STRUCTURE & PARAMETRES DU MODELE
    # Parametres demographiques
    K = parametre[i,1];   # Capacite d'accueil du milieu
    sr = parametre[i,2];	# Sex-ratio
    m1 = parametre[i,3];	# Mortalite naturelle de la classe d'age 1
    m2 = parametre[i,4];	# Mortalite naturelle de la classe d'age 2
    m3 = parametre[i,5];	# Mortalite naturelle de la classe d'age 3
    f2 = parametre[i,6];	# Taux de fecondite de la classe d'age 2
    f3 = parametre[i,7];	# Taux de fecondite de la classe d'age 3
    portee = parametre[i,8];	# Nombre de descendants par reproduction
    t1 = parametre[i,9];	# Taux de maturation de la classe d'age 1
    t2 = parametre[i,10];	# Taux de maturation de la classe d'age 2
    
    # Parametres epidemiologiques
    trans = parametre[i,11]; # Taux de transmission de la maladie
    lat = parametre[i,12];	# Temps d'incubation
    rec = parametre[i,13];	# Temps de recuperation
    loss = parametre[i,14];	# Temps de perte d'immunite
    madd = parametre[i,15];	# Mortalite liee a la maladie
    
    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps), dimnames = list(c("J", "A1", "A2", "total"), c("S", "E", "I", "R")));  # Matrice des effectifs de type [classe,etat,temps]
    
    nvinf <- array(0, dim=c(temps));  # Vecteur de l'incidence journaliere
    
    # Condition initiales (La population est a sa structure d'equilibre, prealablement calculee)
    MAT[1,1,1] <- 27;  # 27 JS
    MAT[2,1,1] <- 23;  # 23 A1S
    MAT[3,1,1] <- 36;  # 3 A2S
    MAT[3,3,1] <- 1;  # 1 A2I
    
    # Effectifs totaux par etat de sante
    MAT[4,1,1] <- sum(MAT[1:3,1,1]);  #S
    MAT[4,2,1] <- sum(MAT[1:3,2,1]);  #E 
    MAT[4,3,1] <- sum(MAT[1:3,3,1]);  #I 
    MAT[4,4,1] <- sum(MAT[1:3,4,1]);  #R 
    
    # SIMULATIONS
    for (t in 1:(temps-1)) {  # Boucle temps
      
      # classe d'age J
      # RQ : les naissances sont non-contaminantes, les nouveaux nes etant dans l'etat S
      N <- sum(MAT[4,,t]);	# taille de la pop en t
      MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t] + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
      MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat) + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
      MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec) + lat*MAT[1,2,t]; 
      MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) + rec*MAT[1,3,t]; 
      
      # classe d'age A1
      MAT[2,1,t+1] <- MAT[1,1,t]*t1 + MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
      MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat) + trans*MAT[2,1,t]*MAT[4,3,t]/N;
      MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec) + lat*MAT[2,2,t];
      MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss) + rec*MAT[2,3,t];
      
      # classe d'age 3
      MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) + loss*MAT[3,4,t];
      MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)	+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
      MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec) + lat*MAT[3,2,t];
      MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss) + rec*MAT[3,3,t];
      
      # Calcule des effectifs totaux par etat de sante
      MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]);  #S
      MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]);  #E
      MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]);  #I
      MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);  #R
      
      # Calcule de l'incidence journaliere
      nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N
      
    }  # fin boucle temps
    
    
    
    # SORTIES PONCTUELLES
    # --- Taux de morbidite
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    
    # --- Incidence finale
    sortie2 <- nvinf[temps]
    
    # --- Pic infectieux (Nombre max d'infecte)
    sortie3 <- max(MAT[4,3,1:temps])
    
    # ---  Prevalence sur la premiere annee
    
    sortie4 <- sum(nvinf[1:365])
    
    # Integration des sorties ponctuelles a leur matrice de sortie
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
    # Integration des effectifs simules a leur liste de sortie
    # sortie.MAT[[i]] <- MAT
    # names(sortie.MAT)[i] <- paste0("scenario_", i)
    
    # # Integration des incidences journalieres a leur matrice de sortie
    # sortie.nvinf[[i]] <- nvinf
    # names(sortie.nvinf)[i] <- paste0("scenario_", i)
    
  }  # Fin de la boucle de scenario
  
  
  
  
  
  
  # Output de la fonction :
  # - Matrice des 4 Sorties ponctuelles
  # - Matrice des effectifs par scenario
  # - Vecteur des incidences journalieres
  
  return(list(indicateur_epidemio = sorties))
  
  
}  # Fin

#####  OAT STANDARD  -------------------------------------------



# On definit les gammes de variation par paramètres (10 valeurs par paramètres)
# --- On definit les gammes (_g) de variations

gamme_params <- 
  list(
    K = seq(from = 80, to = 120, length.out = 10),  # 100
    sr = seq(from = .4, to = .6, length.out = 10),  # .5
    m1 = seq(from = 0.0007, to = 0.003, length.out = 10),  # 0.0014 
    m2 = seq(from = .00015, to = 0.00045, length.out = 10),  # 0.00029
    m3 = seq(from = .000095, to = .0038, length.out = 10),  # 0.0019
    f2 = seq(from = .000095, to = .0038, length.out = 10),   # 0.0019
    f3 = seq(from = .0041, to = .01, length.out = 10),   # 0.0082
    portee = seq(from = 3, to = 7, length.out = 10),  # 5 
    t1 = seq(from = 1/340, to = 1/385, length.out = 10),  # 1 / 365
    t2 = seq(from = 1/340, to = 1/385, length.out = 10),  # 1 / 365
    
    trans =  seq(from = 0.2, to = 0.4, length.out = 10),  # 0.3
    lat =  seq(from = 1/8, to = 1/2, length.out = 10),  # 1 / 5
    rec =  seq(from = 1/25, to = 1/15, length.out = 10),  # 1 / 20
    loss = seq(from = 1/90, to = 1/110, length.out = 10),  # 1 / 100
    madd = seq(from = 0.0005, to = 0.0015, length.out = 10) # 0.001
  )



  
  names_para <- names(gamme_params) # Variable contenant les noms des paramètres
  parametres_initiaux <- PAR
  
# --- On cree une matrice avec en ligne les scenarii et en colonne les paramètres

mat_scenarios <-
  matrix(
    data = parametres_initiaux,
    nrow = 15 * 10,
    ncol = 15,
    dimnames = list(paste0("scenario", 1:(15 * 10)), names_para),
    byrow = T
  )


# --- Implementation des scenarios modifies (un paramètre modifie // original a chaque scenario)
ligne <- 1

for (i in 1:length(gamme_params)){
  # Pour chaque paramètre, implementer dans la matrice de scenario les 10 valeurs modifiees dans 10 scenarii
  mat_scenarios[ligne:(ligne+10-1), i] <- gamme_params[[i]][1:10]
  ligne <- ligne + 10
}


# On modélise chaque scénario avec la fonction modAppli
oat <- as.data.frame(modAppli(mat_scenarios)$indicateur_epidemio)

# On modélise avec les paramètres de référence le modèle m0
m0 <- modAppli(parametres_initiaux)$indicateur_epidemio



# ANALYSE OAT _______________________________
parametres_initiaux



# --- Fonction apply pour calculer l'indice de sensibilité par sortie du modèle par paramètre par scénario
list_as <- apply(
  X = oat,
  MARGIN = 1,
  FUN = function(x){abs(m0 - x)/(x + m0)}  # On pondère par x pour pouvoir comparer les paramètre entre eux
)

list_as <- t(list_as)  # on transpose le dataframe


# --- On associe chacun de ces ecart-type aux paramètres associées modifiées
# --- Puis on sélectionne la moyenne des 10 variances 

# Matrice qui va contenir la moyenne + SE des variances calculées
as_mean <-
  as.data.frame(matrix(
    data = 0,
    nrow = 15,
    ncol = 4,
    dimnames = list(names_para, colnames(m0))
  ))

as_se <-
  as.data.frame(matrix(
    data = 0,
    nrow = 15,
    ncol = 4,
    dimnames = list(names_para, colnames(m0))
  ))

ligne <- 1  # incrémentation ligne pour prendre les paramètres des 10 valeurs de sd pour chaques paramètres

# Moyenne + SE de chaque paramètre représentant l'impact de sa modification par rapport à m0 initial
for (i in 1:15){
  as_mean[i, 1:4] <- apply(list_as[ligne:(i*10),1:4], MARGIN = 2, FUN = mean)
  as_se[i, 1:4] <- apply(list_as[ligne:(i*10),1:4], MARGIN = 2, FUN = sd)
  
  ligne <- i*10 + 1
  
}  # fin boucle 'parametre'


# On ajoute une colonne "type de processus" et "parametre" pour la visualisation
as_mean["processus"] <- c(rep("demo", 10), rep("epidemio", 5))
as_se["processus"] <- c(rep("demo", 10), rep("epidemio", 5))
as_mean["params"] <- names_para
as_se["params"] <- names_para

# On reformate le tableau pour ggplot
as_mean <- pivot_longer(data = as_mean, cols = 1:4, names_to = "indicateur", values_to = "mean")
as_se<- pivot_longer(data = as_se, cols = 1:4, names_to = "indicateur", values_to = "se")


# On visualise les résultats de l'analyse de sensibilité
ggplot(data = as_mean, aes(x = factor(params, levels = names_para), y = mean, fill = indicateur)) +
  geom_bar(position = position_dodge(), stat = "identity", col = "black", lwd = 0.8, width = 0.6) + 
  geom_errorbar(aes(ymin = mean, ymax = mean + as_se$se, group = indicateur), 
                position = position_dodge(width = 0.6), width = 0.2, lwd = 0.8) +
  theme_classic() +
  xlab(NULL) +
  facet_wrap(~indicateur, ncol = 1, scales = "free_y", labeller = labeller(indicateur = c(incidence_t730 = "Incidence (t=730)", pic_infectieux = "Pic infectieux", prevalence_annee_1 = "Prévalence 1er année", tx_morbidite = "Taux de morbidité"))) + # Facetter par "indicateur"
  # theme(strip.text = element_blank() +  # Enlève labels des facettes
  ylab(TeX(paste0("Indice de sensibilité  ", "$(\\sqrt{sigma^2})$"))) +
  ylim(c(0, 0.1)) +
  scale_fill_manual(values=c("#ffffff", "#a8b5ae", "#587064", "#0e3123"), 
                    name=NULL,
                    breaks=c("incidence_t730", "pic_infectieux", "prevalence_annee_1", "tx_morbidite"),
                    labels=c("Incidence (t=730)", "Pic infectieux", "Prévalence 1er année", "Taux de morbidité")) +
  theme(legend.position = "none")

