# Packages
# install.packages("sensibility")
library(sensitivity)
library(ggplot2)
library(latex2exp)
library(numDeriv)
library(tidyverse)
library(latex2exp)
library(viridis)


# IMPORTATION FONCTION DE BASE ET DES VALEURS INITIALES
source("FONCTION_BASE.R")

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



# --- Variable contenant les noms des paramètres
names_para <- names(gamme_params)

  
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
oat <- modAppli(mat_scenarios)

# On modélise avec les paramètres de référence le modèle m0
m0 <- modAppli(parametres_initiaux)



# ANALYSE OAT STANDARD _______________________________

# --- Fonction apply pour calculer l'indice de sensibilité par sortie du modèle par paramètre par scénario
list_as <- apply(
  X = oat,
  MARGIN = 1,
  
  # Indice de sensibilite = (output_initial - output_modifie_i) / (output_initial + output_modifie_i)
  FUN = function(x){abs(m0 - x)/(x + m0)}  
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
    dimnames = list(names_para, c("tx_morbidite", "incidence_t730", "pic_infectieux", "prevalence_annee_1"))
  ))

as_se <-
  as.data.frame(matrix(
    data = 0,
    nrow = 15,
    ncol = 4,
    dimnames = list(names_para, c("tx_morbidite", "incidence_t730", "pic_infectieux", "prevalence_annee_1"))
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
ggplot(data = as_mean, aes(x = factor(params, levels = names_para), y = mean, fill = processus)) +
  geom_bar(position = position_dodge(), stat = "identity", col = "black", lwd = 0.8, width = 0.6) + 
  geom_errorbar(aes(ymin = mean, ymax = mean + as_se$se, group = indicateur), 
                position = position_dodge(width = 0.6), width = 0.2, lwd = 0.8) +
  theme_classic() +
  xlab(NULL) +
  facet_wrap(~indicateur, ncol = 2, nrow = 2, scales = "free_x", labeller = labeller(indicateur = c(incidence_t730 = "Incidence (t=730)", pic_infectieux = "Pic infectieux", prevalence_annee_1 = "Prévalence 1er année", tx_morbidite = "Taux de morbidité"))) + # Facetter par "indicateur"
  # theme(strip.text = element_blank() +  # Enlève labels des facettes
  ylab(TeX(paste0("Indice de sensibilité  ", "$(\\sqrt{sigma^2})$"))) +
  ylim(c(0, 0.1)) +
  scale_fill_manual(values=c("#ffffff","#CD0000"), 
                    name=NULL) +
  theme(legend.position = "none")

