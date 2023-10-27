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

# INITIALISATION 
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

# ANALYSE FAST
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


# FAST (1000 iterations)
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


# DISTRIBUTION DES SORTIES DU MODELES -------------------------------------
# Fast (n=100)
labels_sorties <- c("Taux morbidité (t=730)", "Incidence (t=730)", "Pic infectieux", "Prévalence 1ère année" )

x11(width = 14)
par(mfrow = c(1, 4))
line<-par(lwd=2)
for (i in 1:4) {
  hist(
    sortie_100[, i],
    col = "#DCE1EA",
    freq = T,
    border = "black",
    xlab = NULL,
    ylab = NULL,
    main = labels_sorties[i],
    family = "serif",
    cex.lab = 2,
    cex.axis = 2,
    cex.main = 2.4
  )
}

x11(width = 14)
par(mfrow = c(1, 4))
line<-par(lwd=2)
for (i in 1:4) {
  hist(
    sortie_1000[, i],
    col = "#DCE1EA",
    freq = T,
    border = "black",
    xlab = NULL,
    ylab = NULL,
    main = labels_sorties[i],
    family = "serif",
    cex.lab = 2,
    cex.axis = 2,
    cex.main = 2.4
  )
}















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


# ECHANTILLONNAGE --------------------
# Valeur des paramètres
sample_100 <- as_fast_100$X
sample_1000 <- as_fast_1000$X

tex_label <- c(
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

colnames(sample_100) <- tex_label
colnames(sample_1000) <- tex_label

par(mfrow = c(2, 3))

i = c(1, 6, 15)
j = c(8, 14, 11)

# 100
for(k in 1:3){
  plot(
    sample_100[, i[k]] ~ sample_100[, j[k]],
    type = "p",
    main = NULL,
    xlab = TeX(colnames(sample_100)[j[k]]),
    ylab = TeX(colnames(sample_100)[i[k]]),
    family = "serif",
    cex = 0.01,
    cex.lab = 1.5,
    cex.axis = 1.5
  )
}

title(
  outer = T,
  adj = 0.52,
  line = -3,
  main = "(a) 100",
  font.main = 2,
  cex.main = 1.5
)

# 1000
for(k in 1:3){
  plot(
    sample_1000[, i[k]] ~ sample_1000[, j[k]],
    type = "p",
    main = NULL,
    xlab = TeX(colnames(sample_1000)[j[k]]),
    ylab = TeX(colnames(sample_1000)[i[k]]),
    family = "serif",
    cex = 0.01,
    cex.lab = 1.5,
    cex.axis = 1.5
  )
}

title(
  outer = T,
  adj = 0.52,
  line = -21,
  main = "(b) 1000",
  font.main = 2,
  cex.main = 1.5
)



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



# VISUALISATION THEORIQUE
as_fast_500 <- sensitivity::fast99(model = NULL,
                                   factors = names_para,
                                   n = 500,
                                   q = "qunif",
                                   q.arg = q.arg)

sample_500 <- as_fast_500$X

colnames(sample_500) <- tex_label


i = sample(1:15, size = 1)
j = sample(1:15, size = 1)

par(mfrow = c(1, 1))
plot(
  sample_500[, i] ~ sample_500[, j],
  type = "p",
  main = NULL,
  xlab = "Paramètre i",
  ylab = "Paramètre j",
  family = "serif",  cex = 0.1, cex.lab = 1.5, cex.axis = 1.2
)
