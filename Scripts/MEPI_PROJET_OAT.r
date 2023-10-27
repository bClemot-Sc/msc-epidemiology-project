#####  OAT STANDARD  -------------------------------------------

# Packages
library(sensitivity)
library(ggplot2)
library(tidyverse)
library(latex2exp)  # Texte latex


# IMPORTATION FONCTION DE BASE ET DES VALEURS INITIALES
source("FONCTION_BASE.R")


# INITIALISATION ----------------------------------------------------------
# Gamme de variation de chaque paramètre 
gamme_params <- vector(mode = "list", length = length(ValNominale))
for (i in 1:length(ValNominale)) {
  gamme_params[[i]] <- seq(from = 0.75 * ValNominale[i], 
                           to = 1.25 * ValNominale[i], 
                           length.out = 10)
}



# Variable contenant les noms des paramètres
names_para <- names(ValNominale)

# Matrice avec en ligne les scenarii non-modifiés et en colonne les paramètres
mat_scenarios <-
  matrix(
    data = ValNominale,
    nrow = length(ValNominale) * 10,
    ncol = length(ValNominale),
    dimnames = list(paste0("scenario", 1:(length(ValNominale) * 10)), names_para),
    byrow = T
  )

# Implementation des scenarii modifiés 
ligne <- 1

for (i in 1:length(gamme_params)){
  mat_scenarios[ligne:(ligne+10-1), i] <- gamme_params[[i]][1:10]
  ligne <- ligne + 10
}



# MODELISATION ------------------------------------------------------------
# On modélise chaque scénario avec la fonction modAppli
oat <- modAppli(mat_scenarios)
dimnames(oat) <- list(c(1:150), c("tx_morbidite", "incidence_t730", "pic_infectieux", "prevalence_annee_1"))

# On modélise avec les paramètres de référence le modèle m0
m0 <- modAppli(matrix(ValNominale, nrow = 1))



# INDICE DE BAUER & HAMBY -------------------------------------------------
bauer <- as.data.frame(matrix(data = NA, nrow = length(ValNominale), ncol = 4))

# Calcul de l'indice (voir rapport pour détail)
k = 1
for(i in 1:4){
  for(j in seq(10, length(ValNominale)*10, 10)){
    
    bauer[k, i] <- (max(oat[((j-10)+1):j,i]) - min(oat[((j-10)+1):j,i])) / max(oat[((j-10)+1):j,i])
    k = k + 1
    
  }
  k = 1
}


# VISUALISATION -----------------------------------------------------------
# On renomme les lignes (paramètres) et les colonnes (sorties)
rownames(bauer) <- names_para 
colnames(bauer) <- c(
  "tx_morbidite",
  "incidence_t730",
  "pic_infectieux",
  "prevalence_annee_1"
)


# On ajoute une colonne "type de processus" et "parametre" pour la visualisation
bauer["processus"] <- c(rep("demo", 10), rep("epidemio", 5))
bauer["params"] <- names_para

# On reformate le tableau pour ggplot
bauer_gg <- pivot_longer(data = bauer, cols = 1:4, names_to = "indicateur", values_to = "indice")
bauer_gg$indicateur <-
  factor(
    x = bauer_gg$indicateur,
    levels = c(
      "tx_morbidite",
      "incidence_t730",
      "pic_infectieux",
      "prevalence_annee_1"
    )
  )

# On visualise les résultats de l'analyse de sensibilité
ggplot(data = bauer_gg, aes(x = factor(params, levels = names_para), y = indice, fill = processus)) +
  geom_bar(position = position_dodge(), stat = "identity", col = "black", lwd = 0.8, width = 0.6) + 
  theme_classic() +
  xlab(NULL) +
  facet_wrap(~indicateur, ncol = 2, nrow = 2, scales = "free",  labeller = labeller(indicateur = c(incidence_t730 = "Incidence (t=730)", pic_infectieux = "Pic infectieux", prevalence_annee_1 = "Prévalence 1ère année", tx_morbidite = "Taux de morbidité (t=730)"))) + # Facetter par "indicateur"
  # theme(strip.text = element_blank() +  # Enlève labels des facettes
  ylab("Indice de sensibilité") +
  ylim(c(0, max(bauer_gg$indice))) +
  scale_fill_manual(values=c("#ffffff","#CD0000"), name=NULL) +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "serif")) +
  scale_x_discrete(labels = c(
    TeX("$K$"),
    TeX("$\\psi$"),
    TeX("$m_L$"),
    TeX("$m_J$"),
    TeX("$m_A$"),
    TeX("$f_J$"),
    TeX("$f_A$"),
    TeX("$\\Theta$"),
    TeX("$\\tau_L$"),
    TeX("$tau_J$"),
    TeX("$\\beta$"),
    TeX("$\\sigma$"),
    TeX("$\\gamma$"),
    TeX("$\\lambda$"),
    TeX("$\\mu$")
  ))



# DISTRIBUTION DES SORTIES DU MODELES -------------------------------------
oat <- as.data.frame(oat, makes_name = T)
labels_sorties <- c("Taux morbidité (t=730)", "Incidence (t=730)", "Pic infectieux", "Prévalence 1ère année" )

x11(width = 14)
par(mfrow = c(1, 4))
line<-par(lwd=2)
for (i in 1:4) {
  hist(
    oat[, i],
    col = "#DCE1EA",
    freq = T,
    border = "black",
    xlab = NULL,
    ylab = NULL,
    main = labels_sorties[i],
    family = "serif",
    cex.lab = 2,
    cex.axis = 1.7,
    cex.main = 2.4
  )
  abline(v = m0[i], lty = 2, col = "red")
}


# EFFET PARAMETRES --------------------------------------------------------
# -- Normalisation entre 0 et 1 par sortie

# --- Noms paramètres
modalites <- c(
  "$K$",
  "$\\psi$",
  "$m_L$",
  "$m_J$",
  "$m_A$",
  "$f_J$",
  "$f_A$",
  "$\\Theta$",
  "$\\tau_L$",
  "$tau_J$",
  "$\\beta$",
  "$\\sigma$",
  "$\\gamma$",
  "$\\lambda$",
  "$\\mu$"
)

oat["params"] <- rep(modalites, each = 10)

# --- Modalités en expression LateX
oat$params <- factor(TeX(oat$params), levels = TeX(modalites))


oat["index"] <- rep(1:10, length(ValNominale))
oat_gg <- pivot_longer(data = oat, cols = 1:4, names_to = "sortie", values_to = "as")

# --- Réordonne les modalités
oat_gg$sortie <-
  factor(
    x = oat_gg$sortie,
    levels = c(
      "tx_morbidite",
      "incidence_t730",
      "pic_infectieux",
      "prevalence_annee_1"
    ),
    labels = c("Taux morbidité (t=730)", "Incidence (t=730)", "Pic infectieux", "Prévalence 1ère année" )
  )


for(i in levels(oat_gg$params)){
  print(ggplot(data = oat_gg[oat_gg$params == i,], mapping = aes(x = index, y = as)) +
          facet_wrap(~sortie, ncol = 2, nrow = 2, scales = "free_y") + # Facetter par "indicateur"
          geom_line() +
          theme_classic() +
          scale_x_continuous(breaks = 1:10) +
          xlab(TeX(i)) +
          ylab("Valeur de sortie")
  )
}

j = 1
for(i in levels(oat_gg$sortie)){
  print(
    ggplot(data = oat_gg[oat_gg$sortie == i,], mapping = aes(x = index, y = as)) +
      facet_wrap(~params, ncol = 3, nrow = 5, scales = "free_y",  
                 labeller = label_parsed) + # Facetter par "indicateur"
      geom_line() +
      theme_classic() +
      scale_x_continuous(breaks = 1:10) +
      scale_y_continuous(breaks = NULL) +  # Remove y-axis values
      xlab(TeX(i)) +
      ylab("Valeur de sortie") +
      xlab(levels(oat_gg$sortie)[j])
  )
  j = j + 1
}







# VISUALISATION THEORIQUE -------------------------------------------------
# Visualisation échantillonnage
par(mfrow=c(1,1))
plot(
  NULL,
  ylim = c(0, 11),
  xlim = c(44, 91),
  main = NULL,
  xlab = TeX("$\\phi$"),
  ylab = TeX("$\\beta$"),
  family = "serif")

points(y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), x = rep((44+91)/2, 10), pch = 19)
points(x = seq(45, 90, length.out = 10), y = rep((0+11)/2, 10), pch = 19)
