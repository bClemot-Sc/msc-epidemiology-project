# Importation du main (fonction + parametres)
source("MEPI_PROJET_FONCTION.r")


# Analyse de sensibilite : OAT  -------------------------------------------

# Packages
# install.packages("sensibility")
library(sensitivity)
library(ggplot2)
library(latex2exp)

# On definit les gammes de variation par paramètres (10 valeurs par paramètres)
# --- On definit les gammes (_g) de variations

gamme_params <- 
  list(
    K = seq(from = 50, to = 150, length.out = 10),  # 100
    sr = seq(from = .4, to = .6, length.out = 10),  # .5
    m1 = seq(from = 0, to = 0.01, length.out = 10),  # 0.0014 
    m2 = seq(from = 0, to = 0.002, length.out = 10),  # 0.00029
    m3 = seq(from = 0, to = 0.01, length.out = 10),  # 0.0019
    f2 = seq(from = 0, to = 0.02, length.out = 10),   # 0.0019
    f3 = seq(from = 0, to = 0.03, length.out = 10),   # 0.0082
    portee = ceiling(seq(from = 1, to = 20, length.out = 10)),  # 5 
    t1 = seq(from = 1/30, to = 1/400, length.out = 10),  # 1 / 365
    t2 = seq(from = 1/30, to = 1/400, length.out = 10),  # 1 / 365
    
    trans =  seq(from = 0.2, to = 0.4, length.out = 10),  # 0.3
    lat =  seq(from = 1/2, to = 1/30, length.out = 10),  # 1 / 5
    rec =  seq(from = 1/5, to = 1/35, length.out = 10),  # 1 / 20
    loss = seq(from = 1/10, to = 1/200, length.out = 10),  # 1 / 100
    madd = seq(from = 0, to = 0.01, length.out = 10) # 0.001
  )
  
# --- On cree une matrice avec en ligne les scenarii et en colonne les paramètres
parametres_initiaux <- PAR

mat_scenarios <-
  matrix(
    data = parametres_initiaux,
    nrow = 15 * 10,
    ncol = 15,
    dimnames = list(paste0("scenario", 1:(15 * 10)), names(gamme_params)),
    byrow = T
  )


# --- Implementation des scenarios modifies (un paramètre modifie // original a chaque scenario)
ligne <- 1


for (i in 1:length(gamme_params)){
  # Pour chaque paramètre, implementer dans la matrice de scenario les 10 valeurs modifiees dans 10 scenarii
  mat_scenarios[ligne:(ligne+10-1), i] <- gamme_params[[i]][1:10]
  ligne <- ligne + 10
}


# Création d'un algorithme d'évaluation de la sensibilité par paramètres
# L'idee est de comparer avec le modele initial m0
# On choisit d'étudier l'effet sur l'incidence

# --- Modèle m0
m0 <- modAppli(parametres_initiaux)
m0_incidence <- m0$incidence$scenario_1

# --- Modèles AS
m_as <- modAppli(mat_scenarios)

# --- on extrait la variance de m0 - modele_AS de chaque scenarios
# --- ... pour estimer la variabilité des écarts de prédictions
var_as <- vector(mode = "list", length = 15)
names(var_as) <- names(gamme_params)

list_as <- lapply(
  X = m_as$incidence,
  FUN = function(x)
    ifelse(sd(m0_incidence - x) != 0, sd(m0_incidence - x), 0)
)



# --- On associe chacun de ces ecart-type aux paramètres associées modifiées
# --- Puis on sélectionne la moyenne des 10 variances 

  # Matrice qui va contenir la moyenne + SE des variances calculées
as <-
  as.data.frame(matrix(
    data = 0,
    nrow = 15,
    ncol = 2,
    dimnames = list(names(gamme_params), c("mean_var", "se_var"))
  ))


ligne <- 1  # incrémentation ligne pour prendre les paramètres des 10 valeurs de sd pour chaques paramètres

  # Moyenne + SE de chaque paramètre représentant l'impact de sa modification par rapport à m0 initial
for (i in 1:15){
  as[i, "mean_var"] <- round(mean(unlist(list_as[ligne:(i*10)])), 3)
  as[i, "se_var"] <- round(sd(unlist(list_as[ligne:(i*10)]))/sqrt(10), 3)
  
  ligne <- i*10 + 1
  
}

  # On ajoute une colonne "type de processus" pour la visualisation
as["processus"] <- c(rep("demo", 10), rep("epidemio", 5))


# On visualise les résultats
ggplot(data = as, aes(x = factor(x = row.names(as), levels = row.names(as)),  # Pour conserver le nom + ordre des paramètres du modèle
                      y = mean_var, fill = processus)) +
  geom_bar(position = position_dodge(), stat = "identity", col = "black", lwd = 0.7, width = 0.6) + 
  geom_errorbar(aes(ymin = as$mean_var , ymax = as$mean_var + as$se_var), width = 0.2, lwd = 0.8) +
  theme_classic() +
  xlab(NULL) +
  ylab(TeX(paste0("Indice de sensibilité", "$(\\sqrt{sigma^2})$"))) +
  scale_fill_manual(values=c("#EEE5DE", "#8B8682"), 
                    name="Processus",
                    breaks=c("demo", "epidemio"),
                    labels=c("Démographique", "Epidémiologique")) 





