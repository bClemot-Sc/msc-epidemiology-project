#### Analyse de senibilite avec la methode Morris et le package 'sensitivity'
#_____________________________________________________________________________

# PACKAGES
library(sensitivity)
library(ggplot2)
library(latex2exp)
library(gridExtra)


# IMPORTATION FONCTION DE BASE ET DES VALEURS INITIALES
source("FONCTION_BASE.R")


## Analyse de sensibilit? ? partir de la fonction morris()
AS.morris <- morris(
  model = modAppli,
  factors = c("K","sr","m1","m2","m3","f2","f3","portee","t1","t2","trans","lat","rec","loss","madd"),
  r = 100,
  design = list(
    type = 'oat',
    levels = 6,
    grid.jump = 3
  ),
  binf = c(80,.4,0.0007,.00015, .000095, .000095,.0041,3,1/385,1/385,0.2,1/8,1/25,1/110,0.0005),
  bsup = c(120, .6, .003, .00045, .0038, .0038, .01, 7, 1/340, 1/340, .4, 1/2, 1/15, 1/90, .0015)
)



## Exatraction des paranetres d'interet
# mu star
mu.star <- apply(abs(AS.morris$ee), 3, function(M){apply(M, 2, mean)})
# sigma
sigma <- apply(AS.morris$ee, 3, function(M){apply(M, 2, sd)})


# Dataframe
df.sorties <- data.frame(AS.morris$factors, mu.star, sigma, row.names = NULL)
colnames(df.sorties) <- c("Factors","mu.star1","mu.star2","mu.star3","mu.star4","sigma1","sigma2","sigma3","sigma4")

# Moyenne des sensibilites des 4 sorties
df.sorties$Moy.mu.star <- rowMeans(df.sorties[,2:5])
df.sorties$Moy.sigma <- rowMeans(df.sorties[,6:9])


# Premi?re sortie - tx morbidité
plot1 <- ggplot(data = df.sorties, aes(x = mu.star1, y = sigma1, label = Factors)) +
  geom_point() +
  geom_text(hjust = 1.5, vjust = 1) +
  xlab(TeX("$\\mu^*$")) +
  ylab(TeX("$\\sigma$")) +
  theme_minimal() +
  ggtitle("Taux de morbidité")

# Deuxieme sortie : incidence finale
plot2 <- ggplot(data = df.sorties, aes(x = mu.star2, y = sigma2, label = Factors)) +
  geom_point() +
  geom_text(hjust = 1.5, vjust = 1) +
  xlab(TeX("$\\mu^*$")) +
  ylab(TeX("$\\sigma$")) +
  theme_minimal() +
  ggtitle("Incidence t=730")

# Troisieme sortie : pic infectieux
plot3 <- ggplot(data = df.sorties, aes(x = mu.star3, y = sigma3, label = Factors)) +
  geom_point() +
  geom_text(hjust = 1.5, vjust = 1) +
  xlab(TeX("$\\mu^*$")) +
  ylab(TeX("$\\sigma$")) +
  theme_minimal() +
  ggtitle("Pic infectieux")

# Quatrieme sortie : prevalence premiere annee
plot4 <- ggplot(data = df.sorties, aes(x = mu.star4, y = sigma4, label = Factors)) +
  geom_point() +
  geom_text(hjust = 1.5, vjust = 1) +
  xlab(TeX("$\\mu^*$")) +
  ylab(TeX("$\\sigma$")) +
  theme_minimal() +
  ggtitle("Prévalence 1ère année")


grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

par(mfrow = c(1, 1))
# Moyenne des 4 sorties 
ggplot(data = df.sorties, aes(x = Moy.mu.star, y = Moy.sigma, label = Factors)) +
  geom_point() +
  geom_text(hjust = 1.5, vjust = 1) +
  xlab(TeX("$\\mu^*$")) +
  ylab(TeX("$\\sigma$")) +
  theme_minimal() +
  ggtitle("Sensibilité moyenne")

