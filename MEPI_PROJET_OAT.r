# Importation du main (fonction + parametres)
source("MEPI_PROJET_SCRIPT.r")


# Analyse de sensibilité : OAT  -------------------------------------------

# Packages
install.packages("sensibility")
library(sensitivity)

# On définit les gammes de variation par paramètres (10 valeurs par paramètres)
# --- On définit les gammes (_g) de variations
K_g = seq(from = 50, to = 150, length.out = 10)  # 100
sr_g = seq(from = .4, to = .6, length.out = 10)  # .5

m1_g = seq(from = 0, to = 0.01, length.out = 10)  # 0.0014 
m2_g = seq(from = 0, to = 0.002, length.out = 10)  # 0.00029
m3_g = seq(from = 0, to = 0.01, length.out = 10)  # 0.0019

f2_g = seq(from = 0, to = 0.02, length.out = 10)   # 0.0019
f3_g = seq(from = 0, to = 0.03, length.out = 10)   # 0.0082

portee_g = ceiling(seq(from = 1, to = 20, length.out = 10))  # 5 

t1_g = seq(from = 50, to = 400, length.out = 10)  # 1 / 365
t2_g = seq(from = 50, to = 400, length.out = 10)  # 1 / 365

trans_g =  seq(from = 0.05, to = 0.8, length.out = 10)  # 0.3
lat_g =  seq(from = 1, to = 30, length.out = 10)  # 1 / 5
rec_g =  seq(from = 1, to = 100, length.out = 10)  # 1 / 20
loss_g = seq(from = 10, to = 200, length.out = 10)  # 1 / 100
madd_g = seq(from = 0, to = 0.01, length.out = 10)  # 0.001

  
  # --- On crée une matrice avec en ligne les paramètres et en colonnes les valeurs possibles

