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

# ANALYSE FAST

# --- Noms des parametres
names_para <- c("K","sr","m1","m2","m3","f2","f3","portee","t1","t2","trans","lat","rec","loss","madd")

# --- Liste des parametres de distribution pour chaque parametre
q.arg <- list(
    list(min = 50, max = 150),
    list(min = 0.3, max = 0.7),
    list(min = 0.0007, max = 0.0028),
    list(min = 0.00015, max = 0.00045),
    list(min = 0.00094, max = 0.0038),
    list(min = 0.00094, max = 0.0038),
    list(min = 0.0041, max = 0.0164),
    list(min = 2, max = 8),
    list(min = 1/385, max = 1/345),
    list(min = 1/385, max = 1/345),
    list(min = 0.1, max = 0.5),
    list(min = 1/8, max = 1/2),
    list(min = 1/30, max = 1/10),
    list(min = 1/150, max = 1/50),
    list(min = 0.0005, max = 0.0015)
)



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
plot.fast99 <- function(x, ylim = c(0, 1), main = NULL, names.arg = NULL, ...) {
  S <- rbind(x$D1 / x$V, 1 - x$Dt / x$V - x$D1 / x$V)
  colnames(S) <- colnames(x$X)
  bar.col <- c("white","grey")
  barplot(S, ylim = ylim, col = bar.col, main = main, names.arg = names.arg)
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
par(mfrow = c(2, 2))
plot.fast99(fast_100_sortie1, main = "Taux de morbidite (n=100)", names.arg = tex_label)
plot.fast99(fast_1000_sortie1, main = "Taux de morbidite (n=1000)", names.arg = tex_label)

plot.fast99(fast_100_sortie2, main = "Incidence t=730 (n=100)", names.arg = tex_label)
plot.fast99(fast_1000_sortie2, main = "Incidence t=730 (n=1000)", names.arg = tex_label)

plot.fast99(fast_100_sortie3, main = "Pic infectieux (n=100)", names.arg = tex_label)
plot.fast99(fast_1000_sortie3, main = "Pic infectieux (n=1000)", names.arg = tex_label)

plot.fast99(fast_100_sortie4, main = "Prevalence 1ere annee (n=100)", names.arg = tex_label)
plot.fast99(fast_1000_sortie4, main = "Prevalence 1ere annee (n=1000)", names.arg = tex_label)



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

# N = 10000
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