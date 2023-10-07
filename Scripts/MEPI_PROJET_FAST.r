#### Analyse de senibilite avec la methode FAST et le package 'sensitivity'
#_____________________________________________________________________________

# RESSOURCES
# Script detaille de la fonction FAST : https://rdrr.io/cran/sensitivity/src/R/fast99.R
# Article utilisant la methode FAST : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7790364/


# PACKAGES
library(sensitivity)
library(ggplot2)
library(latex2exp)

# IMPORTATION FONCTION DE BASE ET DES VALEURS INITIALES
source("FONCTION_BASE.R")

# ANALYSE FAST

# --- Noms des parametres
names_para <- c("K","sr","m1","m2","m3","f2","f3","portee","t1","t2","trans","lat","rec","loss","madd")

# --- Liste des parametres de distribution pour chaque parametre
q.arg <- list(list(min = 80, max = 120),
              list(min = .4, max = .6),
              list(min = 0.0007, max = .003),
              list(min = .00015, max = .00045),
              list(min = .000095, max = .0038),
              list(min = .000095, max = .0038),
              list(min = .0041, max = .01),
              list(min = 3, max = 7),
              list(min = 1/385, max = 1/340),
              list(min = 1/385, max = 1/340),
              list(min = 0.2, max = .4),
              list(min = 1/8, max = 1/2),
              list(min = 1/25, max = 1/15),
              list(min = 1/110, max = 1/90),
              list(min = .0005, max = .0015)
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


# --- Visualisation
# Version modifiee de la fonction plot du package fast99()
plot.fast99 <- function(x, ylim = c(0, 1), main = NULL, ...) {
    S <- rbind(x$D1 / x$V, 1 - x$Dt / x$V - x$D1 / x$V)
    colnames(S) <- colnames(x$X)
    bar.col <- c("white","grey")
    barplot(S, ylim = ylim, col = bar.col, main = main)
    # legend("topright", c("main effect", "interactions"), fill = bar.col)
  }

par(mfrow = c(2, 2))
plot.fast99(fast_100_sortie1, main = "Taux de morbidite")
plot.fast99(fast_100_sortie2, main = "Incidence t=730")
plot.fast99(fast_100_sortie3, main = "Pic infectieux")
plot.fast99(fast_100_sortie4, main = "Prevalence 1ere annee")


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

# --- Calcul des variances et indices de Sobol pour chaque sortie pour chaque sortie du modele
fast_1000_sortie1 <- as_fast_1000
fast_1000_sortie2 <- as_fast_1000
fast_1000_sortie3 <- as_fast_1000
fast_1000_sortie4 <- as_fast_1000

tell(x = fast_1000_sortie1, y = sortie_1000[,1])  # taux de morbidite
tell(x = fast_1000_sortie2, y = sortie_1000[,2])  # incidence t=730
tell(x = fast_1000_sortie3, y = sortie_1000[,3])  # pic infectieux
tell(x = fast_1000_sortie4, y = sortie_1000[,4])  # prevalence Aere annee


# --- Visualisation
# Version modifiee de la fonction plot du package fast99()
plot.fast99 <- function(x, ylim = c(0, 1), main = NULL, ...) {
  S <- rbind(x$D1 / x$V, 1 - x$Dt / x$V - x$D1 / x$V)
  colnames(S) <- colnames(x$X)
  bar.col <- c("white","grey")
  barplot(S, ylim = ylim, col = bar.col, main = main)
  # legend("topright", c("main effect", "interactions"), fill = bar.col)
}

par(mfrow = c(2, 2))
plot.fast99(fast_1000_sortie1, main = "Taux de morbidite")
plot.fast99(fast_1000_sortie2, main = "Incidence t=730")
plot.fast99(fast_1000_sortie3, main = "Pic infectieux")
plot.fast99(fast_1000_sortie4, main = "Prevalence 1ere annee")
