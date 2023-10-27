#### Representation graphique du comportement du modele

### Importation de la fonction et des packages
source("MEPI_PROJET_FONCTION.r")
library(ggplot2)
library(reshape2)
library(gridExtra)

### Creation d'un dataframe avec les effectifs pour chaque compartiment
## Recuperation des effectifs a partir de la fonction du modele
effectifs <- Sorties$n$scenario_1
## Creation d'un dataframe vide
df <- as.data.frame(matrix(NA, nrow = 1 * 365, ncol = 1 + 3 * 4 + 4 + 3))
colnames(df) <- c("temps","LS","LE","LI","LR","JS","JE","JI","JR","AS","AE","AI","AR","S","E","I","R","L","J","A")
df[,1] <- 1:(1*365)  # Ajout de la variable temps

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
df$S <- df$LS + df$JS + df$AS
df$E <- df$LE + df$JE + df$AE
df$I <- df$LI + df$JI + df$AI
df$R <- df$LR + df$JR + df$AR

## Calcule des sommes des effectifs par etat de sante
df$L <- df$LS + df$LE + df$LI + df$LR
df$J <- df$JS + df$JE + df$JI + df$JR
df$A <- df$AS + df$AE + df$AI + df$AR

### Representation graphique ...
## ... de l'effectif total des etats de sant?
# Reconfiguration du dataframe avec Reshape2
df.etats <- df[,c("temps","S","E","I","R")]
df.etats <- melt(as.data.frame(df.etats), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.etats) +
  geom_line(aes(x = temps, y = effectif, col = etat), linewidth = 1.2) +
  scale_color_manual(values = c("black","#CD69C9","#EE9A00","#66CDAA")) +
  labs(col = "Etat de santé") +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = "Effectifs", breaks = c(0,20,40,60,80)) +
  theme_classic(base_size = 20) +
  theme(text = element_text(family = "serif"))

## ... de l'effectif total des classes d'age
# Reconfiguration du dataframe avec Reshape2
df.classes <- df[,c("temps","L","J","A")]
df.classes <- melt(as.data.frame(df.classes), id.vars = "temps", variable.name = "classe", value.name = "effectif")

# Representation graphique avec ggplot2
ggplot(data = df.classes) +
  geom_line(aes(x = temps, y = effectif, col = classe), linewidth = 1.2) +
  scale_color_manual(values = c("royalblue2","orangered2","palegreen2")) +
  labs(col = "Classe d'âge") +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,40), breaks = c(0,10,20,30,40)) +
  theme_classic(base_size = 20) +
  theme(text = element_text(family = "serif"))

## ... des effectifs de chaque etat pour la classe d'age J
# Reconfiguration du dataframe avec Reshape2
df.J <- df[,c("temps","LS","LE","LI","LR")]
df.J <- melt(as.data.frame(df.J), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
plotJ <- ggplot(data = df.J) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 1.2) +
  scale_color_manual(values = c("royalblue1","royalblue2","royalblue3","royalblue4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(col = "Etat de santé", linetype = "Etat de santé") +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = "Effectifs", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  theme(text = element_text(family = "serif"))
plotJ

## ... des effectifs de chaque etat pour la classe d'age A1
# Reconfiguration du dataframe avec Reshape2
df.A1 <- df[,c("temps","JS","JE","JI","JR")]
df.A1 <- melt(as.data.frame(df.A1), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
plotA1 <- ggplot(data = df.A1) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 1.2) +
  scale_color_manual(values = c("orangered1","orangered2","orangered3","orangered4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(col = "Etat de santé", linetype = "Etat de santé") +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = " ", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  theme(text = element_text(family = "serif"))
plotA1

## ... des effectifs de chaque etat pour la classe d'age A2
# Reconfiguration du dataframe avec Reshape2
df.A2 <- df[,c("temps","AS","AE","AI","AR")]
df.A2 <- melt(as.data.frame(df.A2), id.vars = "temps", variable.name = "etat", value.name = "effectif")

# Representation graphique avec ggplot2
plotA2 <- ggplot(data = df.A2) +
  geom_line(aes(x = temps, y = effectif, col = etat, linetype = etat), linewidth = 1.2) +
  scale_color_manual(values = c("palegreen1","palegreen2","palegreen3","palegreen4")) +
  scale_linetype_manual(values = c("solid","dotted","dashed","twodash")) +
  labs(linetype = "Etat de santé") +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = "", limits = c(0,38), breaks = c(0,10,20,30)) +
  theme_classic(base_size = 20) +
  theme(text = element_text(family = "serif"))
plotA2

grid.arrange(plotJ, plotA1, plotA2, ncol = 3, nrow = 1)

## ... de l'incidence
# Recuperation des incidences a partir de la sortie du modele
incidences <- Sorties$incidence$scenario_1
# Creation du dataframe
df.incid <- data.frame("temps" = 1:(1*365), "incidence" = incidences[1:365])

# Representation graphique sur ggplot2
ggplot(data = df.incid) +
  geom_line(aes(x = temps, y = incidence), linewidth = 1.2) +
  scale_x_continuous(name = "Temps", breaks = c(0,100,200,300)) +
  scale_y_continuous(name = "Incidence") +
  theme_classic(base_size = 20) +
  theme(text = element_text(family = "serif"))
