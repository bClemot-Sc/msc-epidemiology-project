### Mod?le dont la sensibilit? doit ?tre analys?e dans le cadre du projet MODE-MPI 2023-2024

### Le mod?le est ici d?fini sous forme de fonction pour faciliter vos analyses de sensibilit? (AS)
### La fonction renvoie les sorties ponctuelles qui sont ? analyser dans l'AS

modAppli <- function(parametre){  

  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours) (ici on a donc deux annees, mais on peut directement le changer)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de param?tres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4) # Techniauement on a l'air de pouvoir mettre autant de set de parametres que l'on veut)

  # boucle des sc?narios de l'?chantillonnage de l'AS (car on peut faire entrer plusieurs jeux de parametres pour que le modele tourne plusieurs fois)
  for (i in 1:nrow(parametre)) { 

    # STRUCTURE & PARAMETRES DU MODELE

    # XX (On recupere ici les differents parametres)
    K = parametre[i,1];		# Capacite limite du milieu (Du fait de la structure 1-N/K qui apparait et de sa valeur)
    sr = parametre[i,2];	# Sex ratio
    m1 = parametre[i,3];	# Mortalite naturelle de la classe 1 
    m2 = parametre[i,4];	# Mortalite naturelle de la classe 2
    m3 = parametre[i,5];	# Mortalite naturelle de la classe 3
    f2 = parametre[i,6];	# Taux de fertilité de la classe 2
    f3 = parametre[i,7];	# Taux de fertilité de la classe 3
    portee = parametre[i,8];	# Nb d'enfant par reproduction
    t1 = parametre[i,9];	# Temps de grandir de la classe 1 a la classe 2, de base est de 1/365, donc de un an
    t2 = parametre[i,10];	# Temps de grandir de la classe 2 a la classe 3, de base est de 1/365 donc de un an

    # XX
    trans = parametre[i,11]; # Taux d'infection ? Qui induit le passage des individus sensibles d'etat 1 a des individus infectes d'etat 2
    lat = parametre[i,12];	# Taux d'incubation ? Qui induit le passage des individus d'etat 2 a 3
    rec = parametre[i,13];	# Taux d'infectes devenant immunises
    loss = parametre[i,14];	# Perte d'immunite des individus d'etat 4, etant potentiellement donc des individus immunises
    madd = parametre[i,15];	# Taux d'infectes mourrant dfe la maladie
    
    # L'etat 1 correspond aux individus sensibles ; l'etat 4 correspond aux individus immunises
    # L'etat 2 est probablement un etat de latence ; tandis aqe l'etat 3 est clairement l'etat infectieux

    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'?ge en ligne (derni?re ligne = pop tot), ?tat de sant? en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    # conditions initiales (la population est ? sa structure d'?quilibre, calcul?e par ailleurs)
    MAT[1,1,1] <- 27; # On a 27 petits en bonne sante
    MAT[2,1,1] <- 23; # On a 23 adultes en bonne sante
    MAT[3,1,1] <- 36; # On a 3 vieux en bonne sante
    MAT[3,3,1] <- 1;  # On a 1 vieux infectieux
    # effectifs par ?tat de sant?
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); MAT[4,2,1] <- sum(MAT[1:3,2,1]); MAT[4,3,1] <- sum(MAT[1:3,3,1]); MAT[4,4,1] <- sum(MAT[1:3,4,1]); # Fauiire les totaux sur la derniere ligne

    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
     # classe d'?ge xx
      # RQ : les naissances sont XX, les nouveaux n?s ?tant dans l'?tat XX
      N <- sum(MAT[4,,t]);	# taille de la pop en t
      # Entree et sorties pour la premiere classe d'age, avec tous les etats de sante
	MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N)     + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
	MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			    + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
	MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
	MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
     # classe d'?ge xx (la deuxieme ligne)
	MAT[2,1,t+1] <- MAT[1,1,t]*t1	    + MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N)   + loss*MAT[2,4,t];
	MAT[2,2,t+1] <- MAT[1,2,t]*t1	    + MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
	MAT[2,3,t+1] <- MAT[1,3,t]*t1	  + MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
	MAT[2,4,t+1] <- MAT[1,4,t]*t1	  + MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
     # classe d'?ge xx (lq troisieme ligne)
	MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
	MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
	MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
	MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
     # calcul des effectifs par ?tat de sant? (la quatrieme ligne donc)
	MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
	
	nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N

    }# fin boucle temps

    # sorties ponctuelles ? analyser
    # Les trois premiere ssont calculees avec la variable temps, donc comme on est en sortie de boucle, elles sont calculees pour l'etat a la fin du modele
    # La quatrieme quant a elle utilise nvinf qui stock quelque chose a chaue iteration de la boucle precedente
    # Taux de morbidité
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # Incidence finale
    sortie2 <- nvinf[temps]
    # Pic infectieux
    sortie3 <- max(MAT[4,3,1:temps])
    # Prévalence de la première année
    sortie4 <- sum(nvinf[1:365])
    
    # i etant l'iteration actuelle de scenario = de jeu de parametre, donc, de la meniere aue l'on peut avoir plusieurs jeux de parametres, on peut avoirt plusieurs  jeux de sorties
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
  }# fin boucle sc?narios AS
  return(list(sorties,MAT))
} # fin fonction du mod?le

# END

#### Application de la fonction
ValNominale = c(
  K = 100,
  sr = 0.5,
  m1 = 0.0014,
  m2 = 0.00029,
  m3 = 0.0019,
  f2 = 0.0019,
  f3 = 0.0082,
  portee = 5,
  t1 = 1 / 365,
  t2 = 1 / 365,
  trans = 0.3,
  lat = 1 / 5,
  rec = 1 / 20,
  loss = 1 / 100,
  madd = 0.001
)


PAR <- matrix(ValNominale, nrow = 1)
Sorties <- modAppli(PAR)

# Repr?sentations graphiques ----------------------------------------------

Effectifs <- Sorties[[2]]

## Evolution des effectifs des ?tats de sant?
# Creation d'un df pour les effectifs des ?tats de sant?
df <- matrix(NA, nrow = 2*365, ncol = 5)
colnames(df) <- c("temps","S","E","I","R")
df[,1] <- 1:(2*365)

# Ajout des sommes des etats de sante dans le df
for (i in df[,1]) {
  df[i,2] <- sum(Effectifs[1:3,1,i])
  df[i,3] <- sum(Effectifs[1:3,2,i])
  df[i,4] <- sum(Effectifs[1:3,3,i])
  df[i,5] <- sum(Effectifs[1:3,4,i])
}

plot(NULL, xlim=c(0,2*365), ylim=c(0,100),xlab = "Temps", ylab = "Effectifs")
lines(S ~ temps, data = df, col = "black")
lines(E ~ temps, data = df, col = "blue")
lines(I ~ temps, data = df, col = "red")
lines(R ~ temps, data = df, col = "green")
legend(650, 100, legend=c("S", "E","I","R"),col=c("black","blue", "red","green"), lty=1, cex=0.8)

# Nouveau dataframe avec less effectifs pour les classes d'ages
df2 <- matrix(NA, nrow = 2 * 365, ncol = 1 + 3 * 4)
colnames(df2) <- c("temps","JS","JE","JI","JR","A1S","A1E","A1I","A1R","A2S","A2E","A2I","A2R")
df2[,1] <- 1:(2*365)

for (i in df[,1]) {
  df2[i,2] <- Effectifs[1,1,i] # JS
  df2[i,3] <- Effectifs[1,2,i] # JE
  df2[i,4] <- Effectifs[1,3,i] # JI
  df2[i,5] <- Effectifs[1,4,i] # JR
  df2[i,6] <- Effectifs[2,1,i] # A1S
  df2[i,7] <- Effectifs[2,2,i] # A1E
  df2[i,8] <- Effectifs[2,3,i] # A1I
  df2[i,9] <- Effectifs[2,4,i] # A1R
  df2[i,10] <- Effectifs[3,1,i] # A2S
  df2[i,11] <- Effectifs[3,2,i] # A2E
  df2[i,12] <- Effectifs[3,3,i] # A2I
  df2[i,13] <- Effectifs[3,4,i] # A2R
}

plot(NULL, xlim=c(0,2*365), ylim=c(0,40),xlab = "Temps", ylab = "Effectifs")
lines(JS ~ temps, data = df2, col = "lightgrey")
lines(A1S ~ temps, data = df2, col = "darkgrey")
lines(A2S ~ temps, data = df2, col = "black")
lines(JE ~ temps, data = df2, col = "lightblue")
lines(A1E ~ temps, data = df2, col = "blue")
lines(A2E ~ temps, data = df2, col = "darkblue")
lines(JI ~ temps, data = df2, col = "pink")
lines(A1I ~ temps, data = df2, col = "red")
lines(A2I ~ temps, data = df2, col = "darkred")
lines(JR ~ temps, data = df2, col = "lightgreen")
lines(A1R ~ temps, data = df2, col = "green")
lines(A2R ~ temps, data = df2, col = "darkgreen")