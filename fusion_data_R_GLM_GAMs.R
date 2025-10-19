# Charger les librairies
library(dplyr)
library(terra)


#1.load data files separately
tuna_data <- read.csv("D:/ecological_informatic/tuna_data_namzaf0122.csv")
#cyclone_data <- read.csv("D:/ecological_informatic/cyclone_parameters_2001_2022.csv")
#anticyclone_data <- read.csv("D:/ecological_informatic/anticyclone_parameters_2001_2022.csv")
cyclonic_data <- read.csv("D:/ecological_informatic/cyclonic_parameters_2001_2022.csv")
anticyclonic_data <- read.csv("D:/ecological_informatic/anticyclonic_parameters_2001_2022.csv")
eke_raster <- rast("D:/ecological_informatic/eke_2001_2022.nc")


# 2. Prepare the capture points
# Create a spatial object from the longitudes and latitudes
capture_points <- vect(tuna_data, geom = c("longitude", "latitude"), crs = "epsg:4326")

# 3.  Extract the values
# This is the key step. The `extract` function will retrieve the value from the raster
# at each position defined in `capture_points`.
# (Note: more complex logic is required to select the correct temporal layer)
eke_values <- terra::extract(eke_raster, capture_points)

# 5.Add the extracted values to the original data frame
tuna_data$EKE <- eke_values[,2] # Column 2 usually contains the values


# Now, tuna_data has a new column "EKE" and is ready for analysis.

# Merge the two data frames
# Use left_join() to keep all rows from tuna_data
# and add the matching columns from envirom_data.

# The join key is specified in the 'by' parameter
# Specify all columns that must match.

##merge tuna_data and cyclonic eddie_data
tuna_data1 <-left_join(
  tuna_data,
  cyclonic_data,
  by = c("Year","Month", "latitude", "longitude")
)
# check the result
head(tuna_data1)

#merging tuna_data1 and anticyclonic eddie_tata
tuna_data2 <-left_join(
  tuna_data1,
  anticyclonic_data,
  by = c("Year","Month", "latitude", "longitude")
)
head(tuna_data2)




##it is on 'tuna_data2' that you then apply the GLM and GAM models

library(MASS)
library(mgcv)
library(ggplot2)
library(RODBC)
library(doBy)
library(plotrix)
#---------------------------------------------------------------------
# ÉTAPE 0 : DATA PREPARATION
#---------------------------------------------------------------------

# --- Cleaning and preparation ---

# # It is crucial to treat year and month as factors
# so that the model considers them as separate categories
tuna_data2 <- tuna_data2 %>%
  mutate(
    Year = as.factor(Year),
    Month = as.factor(Month)
  )

#  Show first lines to check
head(tuna_data2)

# Simple method: get a summary of NAs by column
summary(tuna_data2)


write.csv(tuna_data2, "D:/ecological_informatic/tuna_data2.csv")

# More direct method: count NAs by column
sapply(tuna_data2, function(x) sum(is.na(x)))


# To be very precise, let's check the number of complete lines
# for the variables used in YOUR specific template
model_vars <- c("Catch_SWO", "Catch_ALB","Catch_YFT", "Catch_BET", "Year", "Month", "longitude", "latitude", "Effort")
lignes_completes <- complete.cases(tuna_data2[, model_vars])


sum(lignes_completes)

# Loading libraries and data...
# ... (your existing code)

# --- NEW STEP: DATA CLEANING ---

# List all variables needed for the ENTIRE analysis
# (GLM and GAM models)

all_model_vars <- c(
  "Year", "Month", "longitude", "latitude", "Effort",
   "Catch_SWO",  "Catch_ALB",  "Catch_YFT",  "Catch_BET",
  "EKE", "EA1","Radius1","U1", "c1", "SSH1", "RNL1", "EA2", "Radius2","U2", "c2", "SSH2", "RNL2" # Ajoutez ici toutes vos variables méso-échelle
)


# Create a new "clean" dataframe without any rows containing NAs
# in the important columns.
tuna_data_clean <- na.omit(tuna_data2[, all_model_vars])

# Check the dimensions
cat("Dimensions of the original dataset :", dim(tuna_data), "\n")
cat("Dimensions of the cleaned dataset :", dim(tuna_data_clean), "\n")



#---------------------------------------------------------------------
# # STEP 1: CPUE STANDARDIZATION WITH GLM (FOR EACH SPECIES)
#---------------------------------------------------------------------

# The objective is to model catch based on spatio-temporal factors.
# Effort is included via an offset, which allows the catch rate (CPUE) to be modeled directly.

# Note on the formula:
# - `Catch_SWO ~`: SWO catch is the response variable.
# - `Year + Month`: categorical effects of time.
# - `offset(log(Effort))`: This is the key point in CPUE standardization.
# offset, the model analyzes (Catch/Effort),
# i.e., the CPUE.

# --- 1.1 Standardization for Swordfish (SWO) ---
cat("--- GLM Launch for SWO ---\n")
glm_swo <- glm(
  Catch_SWO ~ Year + Month + longitude + latitude + offset(log(Effort)),
  family=gaussian(link="identity"), data = tuna_data_clean
)
# Show model summary
summary(glm_swo)

# Predict the standardized CPUE index and add it to the dataframe
# `type = "response"` gives us the prediction on the catch scale (not the log scale).
tuna_data_clean$glm_swo_pred <- predict(glm_swo, type = "response")


# --- 1.2 Standardization for Albacore (ALB) ---
cat("--- GLM Launch for ALB ---\n")
glm_alb <- glm(
  Catch_ALB ~ Year + Month + longitude + latitude + offset(log(Effort)),
  family=gaussian(link="identity"),data = tuna_data_clean
)
summary(glm_alb)
tuna_data_clean$glm_alb_pred <- predict(glm_alb, type = "response")

# --- 1.3 Standardization for Yellowfin Tuna (YFT) ---
cat("--- GLM Launch for YFT ---\n")
glm_yft <- glm(
  Catch_YFT ~ Year + Month + longitude + latitude + offset(log(Effort)),
  family=gaussian(link="identity"),data = tuna_data_clean
)

summary(glm_yft)

tuna_data_clean$glm_yft_pred <- predict(glm_yft, type = "response")

# --- 1.4 Standardization for Bigeye Tuna (BET) ---
cat("--- GLM Launch for BET ---\n")
glm_bet <- glm(
  Catch_BET ~ Year + Month + longitude + latitude + offset(log(Effort)),
  family=gaussian(link="identity"), data = tuna_data_clean
)
summary(glm_bet)
tuna_data_clean$glm_bet_pred <- predict(glm_bet, type = "response")


# Checking the dataframe with the new prediction columns
head(tuna_data_clean)
write.csv(tuna_data_clean, "D:/ecological_informatic/tuna_data_clean.csv")

#---------------------------------------------------------------------
# STEP 2: MODELING WITH GAM AND MESOSCALE VARIABLES
#---------------------------------------------------------------------
# Now, we use the catches predicted by the GLM (our standardized CPUE index)
# as the response variable in a new model (GAM) to see how they are
# influenced by environmental variables.

# --- Choosing the family for the GAM ---
# The response variable is now `glm_xxx_pred`. It is no longer a count
# but a continuous, positive value.
# - Gaussian: can work and guarantees positive predictions.



#construction of the covariate matrix, correlation, Spearman matrix
library(sf)
library(corrplot)
library(ggcorrplot)
library(Hmisc)
library(RColorBrewer)
library(ComplexHeatmap)

Ms_alb=tuna_data_clean[, c(10,11,12,13,14,15,16,17, 18,19,20,21,22,24)]
Ms_bet=tuna_data_clean[, c(10,11,12,13,14,15,16,17, 18,19,20,21,22,26)]
Ms_swo=tuna_data_clean[, c(10,11,12,13,14,15,16,17, 18,19,20,21,22,23)]
Ms_yft=tuna_data_clean[, c(10,11,12,13,14,15,16,17, 18,19,20,21,22,25)]

alb_M31 <- cor(Ms_alb, method = "spearman", use = "pairwise.complete.obs")
print(round(alb_M31,2))
bet_M31 <- cor(Ms_bet, method = "spearman", use = "pairwise.complete.obs")
print(round(bet_M31,2))
swo_M31 <- cor(Ms_swo, method = "spearman", use = "pairwise.complete.obs")
print(round(swo_M31,2))
yft_M31 <- cor(Ms_yft, method = "spearman", use = "pairwise.complete.obs")
print(round(yft_M31,2))

corrplot(alb_M31, method = "circle", type ="upper", order = "hclust",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         mar= c(0, 0, 1, 0))


resultats_corr <- rcorr(as.matrix(Ms_alb), type = "spearman")
# # Let's extract the matrices that interest us
matrice_r <- resultats_corr$r
matrice_p <- resultats_corr$P

# Optional: Display matrices in the console
 print("Correlation coefficient matrix for ALB (r):")
 print(round(matrice_r, 2))
 print("P-value matrix (P):")
 print(round(matrice_p, 4))


#  Creating the Visualization (Image Style)

# Define a color palette similar to the example (Red - White - Blue)
couleurs <- colorRampPalette(c("#0000AA", "white", "#CC0000"))(200)

corrplot(
  matrice_r,
  method = "color",       # Utilise des carrés de couleur
  type = "full",          # Affiche la matrice complète
  
  # --- Ajout du clustering hiérarchique ---
  order = "hclust",       # Réorganise les variables selon le clustering
  hclust.method = "complete", # Méthode de clustering (d'autres sont "ward.D", etc.)
  addrect = 2,            # Optionnel : Dessine des rectangles autour des groupes
  
  # --- Personnalisation des couleurs et du texte ---
  col = couleurs,             # Applique notre palette de couleurs personnalisée
  tl.col = "black",       # Couleur du nom des variables (x et y)
  tl.srt = 45,            # Rotation du nom des variables en bas
  
  # --- Ajout des coefficients sur le graphique ---
  addCoef.col = "black",  # Couleur des coefficients
  number.cex = 0.7,       # Taille du texte des coefficients
  
  # --- Gestion des p-values pour masquer les corrélations non-significatives ---
  p.mat = matrice_p,      # Fournit la matrice des p-values
  sig.level = 0.05,       # Seuil de significativité (vous pouvez le changer)
  insig = "blank",        # Laisse les cases vides si p > 0.05
  # Alternative : insig = "pch" pour mettre une croix
  
  # --- Ajout d'un titre et gestion des marges ---
  #title = "Matrice de corrélation de Spearman (avec clustering)",
  mar = c(0, 0, 1, 0)       # Marges pour le titre : c(bas, gauche, haut, droite)
)

resultats_corr1 <- rcorr(as.matrix(Ms_bet), type = "spearman")
# Extrayons les matrices qui nous intéressent
matrice_r1 <- resultats_corr1$r
matrice_p1 <- resultats_corr1$P

# Optionnel : Afficher les matrices dans la console
print("Correlation coefficient matrix for BET (r):")
print(round(matrice_r1, 2))
print("Matrice des p-values (P):")
print(round(matrice_p1, 4))

couleurs <- colorRampPalette(c("#0000AA", "white", "#CC0000"))(200)

corrplot(
  matrice_r1,
  method = "color",       # 
  type = "full",          # 
  
  # --- Ajout du clustering hiérarchique ---
  order = "hclust",       # 
  hclust.method = "complete", # 
  addrect = 2,            # 
  
  # --- Personnalisation des couleurs et du texte ---
  col = couleurs,             # Applique notre palette de couleurs personnalisée
  tl.col = "black",       # Couleur du nom des variables (x et y)
  tl.srt = 45,            # Rotation du nom des variables en bas
  
  # --- Ajout des coefficients sur le graphique ---
  addCoef.col = "black",  # Couleur des coefficients
  number.cex = 0.7,       # Taille du texte des coefficients
  
  # --- Gestion des p-values pour masquer les corrélations non-significatives ---
  p.mat = matrice_p1,      # Fournit la matrice des p-values
  sig.level = 0.05,       # Seuil de significativité (vous pouvez le changer)
  insig = "blank",        # Laisse les cases vides si p > 0.05
  # Alternative : insig = "pch" pour mettre une croix
  
  # --- Ajout d'un titre et gestion des marges ---
  #title = "Matrice de corrélation de Spearman (avec clustering)",
  mar = c(0, 0, 1, 0)       # Marges pour le titre : c(bas, gauche, haut, droite)
)

resultats_corr2 <- rcorr(as.matrix(Ms_swo), type = "spearman")
# Extrayons les matrices qui nous intéressent
matrice_r2 <- resultats_corr2$r
matrice_p2 <- resultats_corr2$P

# Optionnel : Afficher les matrices dans la console
print("Correlation coefficient matrix for SWO (r):")
print(round(matrice_r2, 2))
print("Matrice des p-values (P):")
print(round(matrice_p2, 4))

couleurs <- colorRampPalette(c("#0000AA", "white", "#CC0000"))(200)

corrplot(
  matrice_r2,
  method = "color",       # Utilise des carrés de couleur
  type = "full",          # Affiche la matrice complète
  
  # --- Ajout du clustering hiérarchique ---
  order = "hclust",       # Réorganise les variables selon le clustering
  hclust.method = "complete", # Méthode de clustering (d'autres sont "ward.D", etc.)
  addrect = 2,            # Optionnel : Dessine des rectangles autour des groupes
  
  # --- Personnalisation des couleurs et du texte ---
  col = couleurs,             # Applique notre palette de couleurs personnalisée
  tl.col = "black",       # Couleur du nom des variables (x et y)
  tl.srt = 45,            # Rotation du nom des variables en bas
  
  # --- Ajout des coefficients sur le graphique ---
  addCoef.col = "black",  # Couleur des coefficients
  number.cex = 0.7,       # Taille du texte des coefficients
  
  # --- Gestion des p-values pour masquer les corrélations non-significatives ---
  p.mat = matrice_p2,      # Fournit la matrice des p-values
  sig.level = 0.05,       # Seuil de significativité (vous pouvez le changer)
  insig = "blank",        # Laisse les cases vides si p > 0.05
  # Alternative : insig = "pch" pour mettre une croix
  
  # --- Ajout d'un titre et gestion des marges ---
  #title = "Matrice de corrélation de Spearman (avec clustering)",
  mar = c(0, 0, 1, 0)       # Marges pour le titre : c(bas, gauche, haut, droite)
)


resultats_corr3 <- rcorr(as.matrix(Ms_yft), type = "spearman")
# 
matrice_r3 <- resultats_corr3$r
matrice_p3 <- resultats_corr3$P

# 
print("Correlation coefficient matrix for YFT (r):")
print(round(matrice_r3, 2))
print("Matrice des p-values (P):")
print(round(matrice_p3, 4))

couleurs <- colorRampPalette(c("#0000AA", "white", "#CC0000"))(200)

corrplot(
  matrice_r3,
  method = "color",       # 
  type = "full",          # 
  
  # ---  ---
  order = "hclust",       # 
  hclust.method = "complete", # 
  addrect = 2,            # 
  
  # ---  ---
  col = couleurs,             # 
  tl.col = "black",       # 
  tl.srt = 45,            # 
  
  # ---  ---
  addCoef.col = "black",  # 
  number.cex = 0.7,       # 
  
  # ---  ---
  p.mat = matrice_p3,      # 
  sig.level = 0.05,       # 
  insig = "blank",        #
  # 
  
  # ---  ---
  
  mar = c(0, 0, 1, 0)       
)

#PEARSON CORRRELATION
errors_M31=cor(Ms_alb)
bet_M31=cor(Ms_bet)
swo_M31=cor(Ms_swo)
yft_M31=cor(Ms_yft)



corrplot.mixed(alb_M31, 'number')
corrplot.mixed(bet_M31, 'number')
corrplot.mixed(swo_M31, 'number')
corrplot.mixed(yft_M31, 'number')
corrplot(errors_M31, method ='color', title = "LBSPR_Non-Equilibrium")
corrplot.mixed(errors_M31, lower = 'shade', upper = 'pie', order = 'hclust')
corrplot.mixed(errors_M31, 'number')
corrplot.mixed(errors_M31, order = 'AOE')
corrplot.mixed(errors_M31, lower = 'shade', upper = 'pie', order = 'hclust')

corrplot(errors_M31, method = 'number')

library(usdm)

source("http://www.sthda.com/upload/rquery_cormat.r")
rquery.cormat(Ms31, type="flatten", graph=FALSE)
###VIF

vif(Ms31)
vif1=vifstep(Ms31, th=3)
vif2=vifstep(Ms31, th=5)


#chosen parameters("EKE", "EA1","Radius1", "SSH1", "RNL1", "EA2", "Radius2", "SSH2", "RNL2")



# --- 2.1 GAM Model for Swordfish (SWO) ---

# Note on the formula:
# - `s(EKE, k = 5)`: a smoothing function for Eddy kenecty energy (EKE).
# The model will find the shape of the nonlinear relationship between CPUE and EKE.
# `k` controls the flexibility of the curve (adjust if necessary)

cat("--- Launch of GAM for SWO ---\n")


gam_swo_meso1 <- gam(
  glm_swo_pred ~  s(EKE, k=5) + s(EA1, k = 5) + s(Radius1, k = 5) + s(SSH1, k=5) + s(RNL1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_swo_meso1)
plot(gam_swo_meso1, pages = 1, all.terms = TRUE,scheme=1,cex.axis = 1.8, cex.lab = 1.8,  shade.col = "lightblue")
gam.check(gam_swo_meso1)
AIC(gam_swo_meso1)

gam_swo_meso2 <- gam(
  glm_swo_pred ~ s(EKE, k=5) + s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k=5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_swo_meso2)
plot(gam_swo_meso2, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_swo_meso2)
AIC(gam_swo_meso2)

gam_swo_meso3 <- gam(
  glm_swo_pred ~ s(EKE, k = 5) + s(EA2, k = 5) + s(Radius2, k = 5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_swo_meso3)
plot(gam_swo_meso3, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_swo_meso3)
AIC(gam_swo_meso3)

gam_swo_meso4 <- gam(
  glm_swo_pred ~  s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k = 5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_swo_meso4)
plot(gam_swo_meso4, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_swo_meso4)
AIC(gam_swo_meso4)

AIC(gam_swo_meso1, gam_swo_meso2, gam_swo_meso3, gam_swo_meso4)

# --- 2.2 GAM Model for Albacore (ALB) ---
cat("--- Launch of GAM for SWO ALB ---\n")


summary(gam_alb_meso)
plot(gam_alb_meso, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_alb_meso)

gam_alb_meso1 <- gam(
  glm_alb_pred ~ s(EKE, k = 5) + s(EA1, k = 5) + s(Radius1, k = 5) + s(SSH1, k=5) + s(RNL1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_alb_meso1)
plot(gam_alb_meso1, pages = 1, all.terms = TRUE,scheme=1, cex.axis = 2, cex.lab = 2,  shade.col = "lightblue")
gam.check(gam_alb_meso1)

gam_alb_meso2 <- gam(
  glm_alb_pred ~ s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k=5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_alb_meso2)
plot(gam_alb_meso2, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_alb_meso2)


gam_alb_meso3 <- gam(
  glm_alb_pred ~ s(EKE, k = 5) + s(EA2, k = 5) + s(Radius2, k = 5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_alb_meso3)
plot(gam_alb_meso3, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_alb_meso3)

gam_alb_meso4 <- gam(
  glm_alb_pred ~ s(EKE, k = 5) + s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k=5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_alb_meso4)
plot(gam_alb_meso4, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_alb_meso4)


AIC(gam_alb_meso1, gam_alb_meso2, gam_alb_meso3, gam_alb_meso4)


# --- 2.3 GAM Model for Yellowfin Tuna (YFT) ---
cat("--- Launch of GAM for SWO YFT ---\n")

gam_yft_meso1 <- gam(
  glm_yft_pred ~ s(EKE, k = 5) +s(EA1, k = 5)+ s(SSH1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_yft_meso1)
plot(gam_yft_meso1, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_yft_meso1)
AIC(gam_yft_meso1)

gam_yft_meso2 <- gam(
  glm_yft_pred ~ s(EKE, k = 5) + s(SSH1, k=5) + s(RNL1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_yft_meso2)
plot(gam_yft_meso2, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_yft_meso2)
AIC(gam_yft_meso2)

gam_yft_meso3 <- gam(
  glm_yft_pred ~ s(EKE, k = 5) + s(EA2, k = 5) + s(SSH2, k = 5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_yft_meso3)
plot(gam_yft_meso3, pages = 1, all.terms = TRUE,scheme=1,cex.axis = 1.6, cex.lab = 1.6,  shade.col = "lightblue")
gam.check(gam_yft_meso3)
AIC(gam_yft_meso3)

gam_yft_meso4 <- gam(
  glm_yft_pred ~ s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k=5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_yft_meso4)
plot(gam_yft_meso4, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_yft_meso4)



AIC(gam_yft_meso1, gam_yft_meso2, gam_yft_meso3, gam_yft_meso4)


# --- 2.4 Launch of GAM for  Bigeye Tuna (BET) ---
cat("--- Lancement du GAM pour BET ---\n")
gam_bet_meso1 <- gam(
  glm_bet_pred ~ s(EKE, k = 5) + s(EA1, k = 5) + s(SSH1, k=5) + s(RNL1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_bet_meso1)
plot(gam_bet_meso1, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_bet_meso1)

gam_bet_meso2 <- gam(
  glm_bet_pred ~ s(EA1, k = 5) +  s(Radius1, k = 5) + s(SSH1, k=5) + s(RNL1, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_bet_meso2)
plot(gam_bet_meso2, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_bet_meso2)

gam_bet_meso3 <- gam(
  glm_bet_pred ~ s(EKE, k = 5) + s(EA2, k = 5) + s(Radius2, k = 5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_bet_meso3)
plot(gam_bet_meso3, pages = 1, all.terms = TRUE,scheme=1,cex.axis = 1.6, cex.lab = 1.6,  shade.col = "lightblue")
gam.check(gam_bet_meso3)

gam_bet_meso4 <- gam(
  glm_bet_pred ~ s(EA2, k = 5) + s(Radius2, k = 5) + s(SSH2, k=5) + s(RNL2, k=5),
  data = tuna_data_clean,
  family = gaussian(),
  method = "REML"
)
summary(gam_bet_meso4)
plot(gam_bet_meso4, pages = 1, all.terms = TRUE,scheme=1)
gam.check(gam_bet_meso4)

AIC(gam_bet_meso1, gam_bet_meso2, gam_bet_meso3, gam_bet_meso4)







