#### Representation graphique du comportement du modele

### Importation de la fonction et des packages
source("MEPI_PROJET_FONCTION_MODIFIEE.r")
library(ggplot2)
library(reshape2)

### Creation d'un dataframe avec les effectifs pour chaque compartiment
## Recuperation des effectifs a partir de la fonction du modele
effectifs <- Sorties$n$scenario_1
## Creation d'un dataframe vide
df <- as.data.frame(matrix(NA, nrow = 2 * 365, ncol = 1 + 3 * 4 + 4 + 3))
colnames(df) <- c("temps","JS","JE","JI","JR","A1S","A1E","A1I","A1R","A2S","A2E","A2I","A2R","S","E","I","R","J","A1","A2")
df[,1] <- 1:(2*365)  # Ajout de la variable temps

## Obtentien des effectifs pour chaque classe d'age x etat de sante
for (i in df[,1]) {
  df[i,2] <- effectifs[1,1,i]  # Somme des JS
  df[i,3] <- effectifs[1,2,i]  # Somme des JE
  df[i,4] <- effectifs[1,3,i]  # Somme des JI
  df[i,5] <- effectifs[1,4,i]  # Somme des JR
  df[i,6] <- effectifs[2,1,i]  # Somme des A1S
  df[i,7] <- effectifs[2,2,i]  # Somme des A1E
  df[i,8] <- effectifs[2,3,i]  # Somme des A1I
  df[i,9] <- effectifs[2,4,i]  # Somme des A1R
  df[i,10] <- effectifs[3,1,i]  # Somme des A2S
  df[i,11] <- effectifs[3,2,i]  # Somme des A2E
  df[i,12] <- effectifs[3,3,i]  # Somme des A2I
  df[i,13] <- effectifs[3,4,i]  # Somme des A2R
}

## Calcule des sommes des effectifs par classe d'age
df$S <- df$JS + df$A1S + df$A2S
df$E <- df$JE + df$A1E + df$A2E
df$I <- df$JI + df$A1I + df$A2I
df$R <- df$JR + df$A1R + df$A2R

## Calcule des sommes des effectifs par etat de sante
df$J <- df$JS + df$JE + df$JI + df$JR
df$A1 <- df$A1S + df$A1E + df$A1I + df$A1R
df$A2 <- df$A2S + df$A2E + df$A2I + df$A2R

### Representation graphique ...
## ... de l'effectif total des etats de santé
# Reconfiguration du dataframe avec Reshape2
df.etats <- df[,c("temps","S","E","I","R")]
df.etats <- melt(as.data.frame(df.etats), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.etats) +
  geom_line(aes(x = temps, y = effectif, col = etat), linewidth = 0.7) +
  scale_color_manual(values = c("black","#CD69C9","#EE9A00","#66CDAA")) +
  labs(col = "Etat de sante") +
  ggtitle("Effectif total de chaque etat de sante") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs", breaks = c(0,20,40,60,80)) +
  theme_minimal()

## ... de l'effectif total des classes d'age
# Reconfiguration du dataframe avec Reshape2
df.classes <- df[,c("temps","J","A1","A2")]
df.classes <- melt(as.data.frame(df.classes), id.vars = "temps", variable.name = "classe", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.classes) +
  geom_line(aes(x = temps, y = effectif, col = classe), linewidth = 0.7) +
  scale_color_manual(values = c("royalblue2","orangered2","palegreen2")) +
  labs(col = "Classe d'age") +
  ggtitle("Effectif total de chaque classe d'age") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,40), breaks = c(0,10,20,30,40)) +
  theme_minimal()

## ... des effectifs de chaque etat pour la classe d'age J
# Reconfiguration du dataframe avec Reshape2
df.J <- df[,c("temps","JS","JE","JI","JR")]
df.J <- melt(as.data.frame(df.J), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.J) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 0.7) +
  scale_color_manual(values = c("royalblue1","royalblue2","royalblue3","royalblue4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(col = "Etat de sante", linetype = "Etat de sante") +
  ggtitle("Effectifs de chaque etat pour la classe d'age J") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_minimal()

## ... des effectifs de chaque etat pour la classe d'age A1
# Reconfiguration du dataframe avec Reshape2
df.A1 <- df[,c("temps","A1S","A1E","A1I","A1R")]
df.A1 <- melt(as.data.frame(df.A1), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.A1) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 0.7) +
  scale_color_manual(values = c("orangered1","orangered2","orangered3","orangered4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(col = "Etat de sante", linetype = "Etat de sante") +
  ggtitle("Effectifs de chaque etat pour la classe d'age A1") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_minimal()

## ... des effectifs de chaque etat pour la classe d'age A2
# Reconfiguration du dataframe avec Reshape2
df.A2 <- df[,c("temps","A2S","A1E","A2I","A2R")]
df.A2 <- melt(as.data.frame(df.A2), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.A2) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 0.7) +
  scale_color_manual(values = c("palegreen1","palegreen2","palegreen3","palegreen4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(col = "Etat de sante", linetype = "Etat de sante") +
  ggtitle("Effectifs de chaque etat pour la classe d'age A2") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_minimal()

## ... de l'incidence
# Recuperation des incidences a partir de la sortie du modele
incidences <- Sorties$incidence$scenario_1
# Creation du dataframe
df.incid <- data.frame("temps" = 1:(2*365), "incidence" = incidences)

# Representation graphique sur ggplot2
ggplot(data = df.incid) +
  geom_line(aes(x = temps, y = incidence), size = 0.7) +
  ggtitle("Evolution de l'incidence journaliere en fonction du temps") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Incidence") +
  theme_minimal()

## ... des effectifs des pathogènes libres
# Recuperation des incidences a partir de la sortie du modele
eff.pl <- Sorties$pathogene
# Creation du dataframe
df.pl <- data.frame("temps" = 1:(2*365), "pathogene" = eff.pl)

# Representation graphique sur ggplot2
ggplot(data = df.pl) +
  geom_line(aes(x = temps, y = pathogene), size = 0.7) +
  ggtitle("Effectif des pathogènes libres dans l'environnement") +
  scale_x_continuous(name = "Temps", breaks = c(0,150,300,450,600,750)) +
  scale_y_continuous(name = "Effectifs pathogenes libres") +
  theme_minimal()

