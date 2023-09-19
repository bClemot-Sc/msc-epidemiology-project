#### Analyse de senibilite avec la methode Morris et le package 'sensitivity'

## Creation de la fonction du modele :
modAppli <- function(parametre){  
  
  # CONDITIONS DE SIMULATION
  temps = 2*365; # nb de pas de temps (en jours)
  # initialisation pour la sauvegarde de 4 sorties ponctuelles pour chaque jeu de paramètres
  sorties <- matrix(0, nrow=nrow(parametre), ncol=4)
  
  # boucle des scénarios de l'échantillonnage de l'AS
  for (i in 1:nrow(parametre)) { 
    
    # STRUCTURE & PARAMETRES DU MODELE
    
    # XX
    K = parametre[i,1];		# xx
    sr = parametre[i,2];	# xx
    m1 = parametre[i,3];	# xx
    m2 = parametre[i,4];	# xx
    m3 = parametre[i,5];	# xx
    f2 = parametre[i,6];	# xx
    f3 = parametre[i,7];	# xx
    portee = parametre[i,8];	# xx
    t1 = parametre[i,9];	# xx
    t2 = parametre[i,10];	# xx
    
    # XX
    trans = parametre[i,11]; # xx
    lat = parametre[i,12];	# xx
    rec = parametre[i,13];	# xx
    loss = parametre[i,14];	# xx
    madd = parametre[i,15];	# xx
    
    # INITIALISATION
    MAT <- array(0, dim=c(4,4,temps)); # nb indiv par classe d'âge en ligne (dernière ligne = pop tot), état de santé en colonne, pas de temps (dimension 3)
    nvinf <- array(0, dim=c(temps));
    # conditions initiales (la population est à sa structure d'équilibre, calculée par ailleurs)
    MAT[1,1,1] <- 27; # xx
    MAT[2,1,1] <- 23; # xx
    MAT[3,1,1] <- 36; # xx
    MAT[3,3,1] <- 1;  # xx
    # effectifs par état de santé
    MAT[4,1,1] <- sum(MAT[1:3,1,1]); MAT[4,2,1] <- sum(MAT[1:3,2,1]); MAT[4,3,1] <- sum(MAT[1:3,3,1]); MAT[4,4,1] <- sum(MAT[1:3,4,1]);
    
    # SIMULATIONS
    # boucle du temps
    for (t in 1:(temps-1)) { 
      # classe d'âge xx
      # RQ : les naissances sont XX, les nouveaux nés étant dans l'état XX
      N <- sum(MAT[4,,t]);	# taille de la pop en t
      MAT[1,1,t+1] <- MAT[1,1,t]*(1-m1-t1-trans*MAT[4,3,t]/N) + loss*MAT[1,4,t]      + max(0, sr*portee*(sum(MAT[2,,t])*f2 + sum(MAT[3,,t])*f3) * (1 - N/K)); 
      MAT[1,2,t+1] <- MAT[1,2,t]*(1-m1-t1-lat)			  + trans*MAT[1,1,t]*MAT[4,3,t]/N; 
      MAT[1,3,t+1] <- MAT[1,3,t]*(1-m1-madd-t1-rec)  		  + lat*MAT[1,2,t]; 
      MAT[1,4,t+1] <- MAT[1,4,t]*(1-m1-t1-loss) 		  + rec*MAT[1,3,t]; 
      # classe d'âge xx
      MAT[2,1,t+1] <- MAT[1,1,t]*t1	+ MAT[2,1,t]*(1-m2-t2-trans*MAT[4,3,t]/N) + loss*MAT[2,4,t];
      MAT[2,2,t+1] <- MAT[1,2,t]*t1	+ MAT[2,2,t]*(1-m2-t2-lat)			+ trans*MAT[2,1,t]*MAT[4,3,t]/N;
      MAT[2,3,t+1] <- MAT[1,3,t]*t1	+ MAT[2,3,t]*(1-m2-madd-t2-rec)		+ lat*MAT[2,2,t];
      MAT[2,4,t+1] <- MAT[1,4,t]*t1	+ MAT[2,4,t]*(1-m2-t2-loss)			+ rec*MAT[2,3,t];
      # classe d'âge xx
      MAT[3,1,t+1] <- MAT[2,1,t]*t2	+ MAT[3,1,t]*(1-m3-trans*MAT[4,3,t]/N) 	+ loss*MAT[3,4,t];
      MAT[3,2,t+1] <- MAT[2,2,t]*t2	+ MAT[3,2,t]*(1-m3-lat)				+ trans*MAT[3,1,t]*MAT[4,3,t]/N;
      MAT[3,3,t+1] <- MAT[2,3,t]*t2	+ MAT[3,3,t]*(1-m3-madd-rec)			+ lat*MAT[3,2,t];
      MAT[3,4,t+1] <- MAT[2,4,t]*t2	+ MAT[3,4,t]*(1-m3-loss)			+ rec*MAT[3,3,t];
      # calcul des effectifs par état de santé
      MAT[4,1,t+1] <- sum(MAT[1:3,1,t+1]); MAT[4,2,t+1] <- sum(MAT[1:3,2,t+1]); MAT[4,3,t+1] <- sum(MAT[1:3,3,t+1]); MAT[4,4,t+1] <- sum(MAT[1:3,4,t+1]);
      nvinf[t+1]   <- trans*MAT[4,1,t]*MAT[4,3,t]/N
      
    }# fin boucle temps
    
    # sorties ponctuelles à analyser
    # XX
    sortie1 <- (MAT[4,2,temps]+MAT[4,3,temps])/sum(MAT[4,,temps])
    # xx
    sortie2 <- nvinf[temps]
    # xx
    sortie3 <- max(MAT[4,3,1:temps])
    # xx
    sortie4 <- sum(nvinf[1:365])
    
    sorties[i,1] <- sortie1;
    sorties[i,2] <- sortie2;
    sorties[i,3] <- sortie3;
    sorties[i,4] <- sortie4;
    
  }# fin boucle scénarios AS
  return(sorties)
} # fin fonction du modèle

# END

## Test de la fonction du modele

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

PAR <- matrix(ValNominale, nrow=1)
modAppli(PAR)



## Chargement du package 
library(sensitivity)

## Analyse de sensibilité à partir de la fonction morris()
AS.morris <- morris(
  model = modAppli,
  factors = c("K","sr","m1","m2","m3","f2","f3","portee","t1","t2","trans","lat","rec","loss","madd"),
  r = 5,
  design = list(
    type = 'oat',
    levels = 6,
    grid.jump = 3
  ),
  binf = c(80,0.3,0.007,0.0002,0.001,0.001,0.007,3,1/350,1/350,0.2,1/3,1/15,1/80,0.005),
  bsup = c(120,0.8,0.0021,0.00038,0.0028,0.0028,0.0094,7,1/370,1/370,0.4,1/7,1/25,1/120,0.0015)
)

print(AS.morris)

## Representation du graphique de Morris
plot(AS.morris)


