#### Analyse de senibilite avec la methode FAST et le package 'sensitivity'
#_____________________________________________________________________________

# RESSOURCES
# Script detaille de la fonction FAST : https://rdrr.io/cran/sensitivity/src/R/fast99.R
# Article utilisant la methode FAST : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7790364/


# PACKAGES
library(sensitivity)
library(ggplot2)
library(latex2exp)

# WORKING DIRECTORY
setwd("C:/Users/p_a_8/Documents/GitHub/MEPI-Projet/Scripts")  # PA

# IMPORTATION FONCTION DE BASE ET DES VALEURS INITIALES
source("FONCTION_BASE.R")


# ANALYSE FAST ------------------------------------------------------------


# INITIALISATION ----------------------------------------------------------
# --- Noms des parametres
names_para <- c("K","sr","m1","m2","m3","f2","f3","portee","t1","t2","trans","lat","rec","loss","madd")

# --- Liste des parametres de distribution pour chaque parametre
q.arg <-
  apply(cbind(0.75 * unname(ValNominale),
              1.25 * unname(ValNominale)),
        MARGIN = 1,
        function(x) {
          list(min = x[1], max = x[2])
        })

# ANALYSE FAST ------------------------------------------------------------
# FAST (100 iterations)
# --- Generation des valeurs de parametres pour les differents scenarii 
as_fast_100 <- sensitivity::fast99(model = NULL,
                    factors = names_para,
                    n = 100, 
                    q = "qunif",
                    q.arg = q.arg)



# --- On run le modele modAppli sur ces valeurs
sortie_100 <- modAppli(as_fast_100$X)
colnames(sortie_100) <- c("taux morbidité", "incidence t=730", "pic infectieux", "prévalene 1ère année")


# --- Calcul des variances et indices de Sobol pour chaque sortie pour chaque sortie du modele
fast_100_sortie1 <- as_fast_100
fast_100_sortie2 <- as_fast_100
fast_100_sortie3 <- as_fast_100
fast_100_sortie4 <- as_fast_100

# (tell() calcul à partir des scenarii les indices de Sobol et complète l'objet 'incomplet' x)
tell(x = fast_100_sortie1, y = sortie_100[,1])  # taux de morbidite
tell(x = fast_100_sortie2, y = sortie_100[,2])  # incidence t=730
tell(x = fast_100_sortie3, y = sortie_100[,3])  # pic infectieux
tell(x = fast_100_sortie4, y = sortie_100[,4])  # prevalence Aere annee



# _________________________________________________________


# FAST (10000 iterations)
# --- Generation des valeurs de parametres pour les differents scenarii 
as_fast_1000 <- sensitivity::fast99(model = NULL,
                                   factors = names_para,
                                   n = 1000, 
                                   q = "qunif",
                                   q.arg = q.arg)


# --- On run le modele modAppli sur ces valeurs
sortie_1000 <- modAppli(as_fast_1000$X)
colnames(sortie_1000) <- c("taux morbidité", "incidence t=730", "pic infectieux", "prévalene 1ère année")

# --- Calcul des variances et indices de Sobol pour chaque sortie pour chaque sortie du modele
fast_1000_sortie1 <- as_fast_1000
fast_1000_sortie2 <- as_fast_1000
fast_1000_sortie3 <- as_fast_1000
fast_1000_sortie4 <- as_fast_1000

tell(x = fast_1000_sortie1, y = sortie_1000[,1])  # taux de morbidite
tell(x = fast_1000_sortie2, y = sortie_1000[,2])  # incidence t=730
tell(x = fast_1000_sortie3, y = sortie_1000[,3])  # pic infectieux
tell(x = fast_1000_sortie4, y = sortie_1000[,4])  # prevalence Aere annee





# VISUALISATION -----------------------------------------------------------
# --- Version modifiee de la fonction plot du package fast99()
plot.fast99 <-
  function(x,
           ylim = c(0, 1),
           main = NULL,
           names.arg = NULL,
           ...) {
    S <- rbind(x$D1 / x$V, 1 - x$Dt / x$V - x$D1 / x$V)
    colnames(S) <- colnames(x$X)
    bar.col <- c("white", "grey")
    barplot(
      S,
      ylim = ylim,
      col = bar.col,
      main = main,
      names.arg = names.arg,
      family = "serif",
      cex.axis = 1.4, cex.main = 1.7
      
    )
    # legend("topright", c("main effect", "interactions"), fill = bar.col)
  }

# Noms modifiés des paramètres 
tex_label <- c(
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
  TeX("$\\beta_1$"),
  TeX("$\\sigma$"),
  TeX("$\\gamma$"),
  TeX("$\\lambda$"),
  TeX("$\\mu$")
)

# Indice de sensibilité par paramètre -----------------
x11(width = 14)  # Largeur de la fenêtre
par(mfrow = c(2, 4))
plot.fast99(fast_100_sortie1, main = "Taux de morbidité (t=730)", names.arg = tex_label, family = "serif")
plot.fast99(fast_100_sortie2, main = "Incidence (t=730)", names.arg = tex_label, family = "serif")
plot.fast99(fast_100_sortie3, main = "Pic infectieux", names.arg = tex_label, family = "serif")
plot.fast99(fast_100_sortie4, main = "Prevalence 1ère année", names.arg = tex_label, family = "serif")

plot.fast99(fast_1000_sortie1, main = "", names.arg = tex_label, family = "serif")
plot.fast99(fast_1000_sortie2, main = "", names.arg = tex_label, family = "serif")
plot.fast99(fast_1000_sortie3, main = "", names.arg = tex_label, family = "serif")
plot.fast99(fast_1000_sortie4, main = "", names.arg = tex_label, family = "serif")


# Echantillonnage --------------------

# as_fast_100_scale <- apply(as_fast_100$X, 2, FUN = function(x) scale(x, center = T, scale = T))
# as_fast_1000_scale <- apply(as_fast_1000$X, 2, FUN = function(x) scale(x, center = T, scale = T))
# 
# # plot(
#   NULL,
#   type = "l",
#   ylab = "Valeur du paramètre",
#   xlab = "Scénario",
#   ylim = c(min(as_fast_100_scale), max(as_fast_100_scale)),
#   xlim = c(0, nrow(as_fast_100_scale))
# )


# N = 100
sample_100 <- as_fast_100$X

par(mfrow = c(2, 2))
for (i in sample(1:15, size = 4)) {
  j = sample(1:15, size = 1)
  plot(
    sample_100[, i] ~ sample_100[, j],
    type = "p",
    ylab = paste0(colnames(sample_100)[i]),
    xlab = paste0(colnames(sample_100)[j]),
    main = paste0(colnames(sample_100)[i], " vs ", colnames(sample_100)[j])
  )
  
}

# N = 1000
sample_1000 <- as_fast_1000$X

par(mfrow = c(2, 2))
for (i in sample(1:15, size = 4)) {
  j = sample(1:15, size = 1)
  plot(
    sample_1000[, i] ~ sample_1000[, j],
    type = "p",
    ylab = paste0(colnames(sample_1000)[i]),
    xlab = paste0(colnames(sample_1000)[j]),
    main = paste0(colnames(sample_1000)[i], " vs ", colnames(sample_1000)[j])
  )
  
}

apply(X = sample_1000, MARGIN = 2, FUN = function(x) length(unique(x)))





# SORTIES -----------------------------------------------------------------
# N = 100
par(mfrow = c(2, 2))
for (i in 1:15) {
  for(j in 1:4){
    plot(
      sortie_100[,j] ~ sample_100[, i],
      type = "l",
      ylab = paste0(colnames(sortie_100)[j]),
      xlab = paste0(colnames(sample_100)[i]),
      main = paste0(colnames(sortie_100)[j], " ~ ", colnames(sample_100)[i])
    )
    
  }
}


# N = 1000
for (i in 1:15) {
  for (j in 1:4) {
    plot(
      sortie_1000[, j] ~ sample_1000[, i],
      type = "l",
      ylab = paste0(colnames(sortie_1000)[j]),
      xlab = paste0(colnames(sample_1000)[i]),
      main = paste0(colnames(sortie_1000)[j], " ~ ", colnames(sample_1000)[i])
    )
  }
}


# 
# # Visualisation théorique
# as_fast_100 <- sensitivity::fast99(model = NULL,
#                                    factors = names_para,
#                                    n = 100, 
#                                    q = "qunif",
#                                    q.arg = q.arg)
# 
# as_fast_500 <- sensitivity::fast99(model = NULL,
#                                    factors = names_para,
#                                    n = 500, 
#                                    q = "qunif",
#                                    q.arg = q.arg)
# 
# as_fast_1000 <- sensitivity::fast99(model = NULL,
#                                    factors = names_para,
#                                    n = 1000, 
#                                    q = "qunif",
#                                    q.arg = q.arg)
# 
# sample_100 <- as_fast_100$X
# sample_500 <- as_fast_500$X
# sample_1000 <- as_fast_1000$X
# 
# i = sample(1:15, size = 1)
# j = sample(1:15, size = 1)
# 
# par(mfrow = c(2, 3))
# 
# plot(
#     sample_100[, i] ~ sample_100[, j],
#     type = "p",
#     main = "",
#     xlab = "Paramètre i",
#     ylab = "Paramètre j", 
#     family = "serif",  cex = 0.2
#     )
#   
# 
# plot(
#   sample_500[, i] ~ sample_500[, j],
#   type = "p",
#   main = "",
#   xlab = "Paramètre i",
#   ylab = "Paramètre j", 
#   family = "serif",  cex = 0.2
# )
# 
# plot(
#   sample_1000[, i] ~ sample_1000[, j],
#   type = "p",
#   main = "",
#   xlab = "Paramètre i",
#   ylab = "Paramètre j", 
#   family = "serif",  cex = 0.2
# )
# 
# i = sample(1:15, size = 1)
# j = sample(1:15, size = 1)
# 
# 
# plot(
#   sample_100[, i] ~ sample_100[, j],
#   type = "p",
#   main = "",
#   xlab = "Paramètre k",
#   ylab = "Paramètre l", 
#   family = "serif",  cex = 0.2
# )
# 
# 
# plot(
#   sample_500[, i] ~ sample_500[, j],
#   type = "p",
#   main = "",
#   xlab = "Paramètre k",
#   ylab = "Paramètre l", 
#   family = "serif",  cex = 0.2
# )
# 
# plot(
#   sample_1000[, i] ~ sample_1000[, j],
#   type = "p",
#   main = "",
#   xlab = "Paramètre k",
#   ylab = "Paramètre l", 
#   family = "serif",  cex = 0.2
# )
