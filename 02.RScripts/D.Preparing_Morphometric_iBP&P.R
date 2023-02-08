
######################## SPECIES DELIMITATION ANALYSES ########################
################################## TUTORIAL ###################################


######################### JERONYMO DALAPICOLLA, 2021 ##########################
################################     iBPP    ##################################
######################## STEP 01: MORPHOLOGICAL DATA ##########################


##1. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/

##2. INPUTS FOR THIS STEP:
#A. MORPHOMETRIC DATA
#B. FUNCTIONS .R IN THE SAME FOLDER. GO TO THE GITHUB TO DOWNLOAD THEM

##3. INSTALL AND LOAD THE R PACKAGES:
#install the packages
install.packages(c("tcltk","pacman", "dplyr", "adegenet", "usdm"))
#load the packages
#to load multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load("tcltk", "dplyr", "adegenet", "usdm")

##4. CHOOSE A FOLDER FOR ANALYSIS. THE FILE MUST BE THERE! 
#Windows
choose.dir()
#Linux. You will need the tcltk package for this command:
setwd(tk_choose.dir())

##5. LOAD THE FUNCTIONS FOR THIS SCRIPTS. THE 'R' FILES MUST BE IN THE FOLDER CHOSE ON THE PREVIOUS STEP.
source("id.outliers.R")
source("contar.NA.R")
source("normalidade.R")

##6. LOAD THE FILE WITH MORPHOMETRICS DATA. IN MY FILE THERE ARE 39 COLUNMS, 1:17 CATEGORICAL VARIABLES AND 18:39 THE MORPHOMETRICS MEASUREMENTS.  
tab = read.csv(file="morphometric_data.csv")
#verify if the load of the file is ok.
head(tab)
tail(tab)
summary(tab)

##7. CHOOSE ONLY ADULTS FOR THE ANALYSIS. AGE_PATTON EQUALS TO 8,9, AND 10 FOR PROECHIMYS.
#copy the original data, for security proposals
tab_save = tab

#selecting the adults
tab = tab[tab$age_patton=="8" | tab$age_patton=="9" | tab$age_patton=="10", ]


##8. IDENTIFY THE OUTLIERS INSIDE THE MORPHOGROUPS/SPECIES/SPECIES GROUP AND IN ALL SAMPLES TOGETHER:
#all samples together
names(tab)
out = id.outliers(x=tab, quant=12:33, group=0, id="z", NUMBER=100, visual="biplot", res="MED", csv=T)


#########################################################
#########################################################
#verify if the samples indicated by the test are really outliers, if yes, replace them by "NA". Replace in the original file too. 
#########################################################
#########################################################



##9. LOAD THE FILE WITHOUT OUTLIERS!
tab = read.csv(file="morphometric_data_edited.csv")
#verify if the load of the file is ok.
head(tab)
tail(tab)
summary(tab)

##10. CHOOSE ONLY ADULTS FOR THE ANALYSIS. AGE_PATTON EQUALS TO 8,9, AND 10 FOR PROECHIMYS.
tab = tab[tab$age_patton=="8" | tab$age_patton=="9" | tab$age_patton=="10", ]


##11. REMOVING MISSING DATA
#original data
dim(tab[12:33]) ##314 inds and 22 variables

###Count NA
sum(is.na(tab[12:33])) #21 NA's
summary(tab[12:33])

#Variables
# MD = 6 * removed
# GSL = 3
# NL = 3
# RL = 3
# ZB = 1
# OccW = 1
# IFW = 1
# Pla = 1
# CIL = 1

#Individuals. check if the ind with more NA could be deleted (if you have 2 or + in its population)
inds_to_check = which(apply(tab[12:33], 1, function(x)sum(is.na(x))) != 0)
tab[inds_to_check,]

# 34 =1
# 42 =1
# 63 = 4 #removed
# 102 = 1 
# 126 = 1
# 199 = 1
# 200 = 1
# 211 = 1
# 224 = 3 # Removed
# 226 = 1
# 244 = 3 # Removed
# 252 = 1
# 264 = 1
# 294 = 1

tab = tab[-c(63,224,244),c(-33)] #311 inds and 21 variables
sum(is.na(tab[12:32])) #5
summary(tab[12:32])

#Variables with NA and check if the ind could be deleted (if you have 2 ou + by population)
# ZB = 1 #removed ind D1
# OccW = 1 #removed ind C2
# IFW = 1 #removed ind C2
# Pla = 1 # removed ind B1
# CIL = 1 # removed ind C1

#Individuals. check if the ind with more NA could be deleted (if you have 2 or + in its population)
inds_to_check = which(apply(tab[12:32], 1, function(x)sum(is.na(x))) != 0)
tab[inds_to_check,]

tab = tab[-c(101,198,210,224,291),]
sum(is.na(tab[12:32])) #0
summary(tab[12:32])




##12. TRANSFORM THE DATA IN LOG. DO NOT STANDARD THE DATA. iBP&P WILL DO IT LATTER.
#save a copy of raw data
data_raw = tab

#save a table with log transformed data
tab[12:32] = log(tab[12:32])

#test normality with qqplot
norm = normalidade (x=tab[12:32], quant=1:21, group=0)
names(tab)

## not normal by qqplots
tab_red = tab[,-c(25,26,27,32)]
names(tab_red)



##13. REMOVING ALLOMETRY WITH PCA WITH SIZE
#center = TRUE means centring by the mean, in this case is 0
#scannf	= FALSE, indicates the screeplot should not be displayed
#nf = 3, if scannf = FALSE, an integer must be indicated. It will be the number of kept axes in the PCA
pca_input = dudi.pca(tab_red[12:28], center = TRUE , scannf = FALSE, nf = 20)

pca_input$c1 #principal axes
pca_input$li #principal components

#Variable with higher contribution to PC1 (size) 
max(abs(pca_input$c1[,1]))
#0.2676335
#CIL

##11. REGRESSION TO REMOVE THE SIZE EFFECT:
reg = tab_red[12:28]
names(reg)
residuo1 = matrix(NA, 306, 17)


for(i in 1:length(reg)){
  residuo1[,i] = unlist(lm(reg[,i] ~ reg$CIL)[2])
}

dados2 = as.data.frame(residuo1)
names(dados2) = names(reg)
names(dados2)
data_log_pca = cbind(tab_red[1:11],dados2)


#Check CORRELATED VARIABLES TO REDUCE NUMBER OF :
vifcor(data_log_pca[12:28],th = 0.3) #11 vars ok
#max correlation ( PLa ~ BaL ):   0.2998007


vars_final = vifcor(data_log_pca[12:28],th = 0.3)
dados3 = data_log_pca[12:28][,vars_final@results$Variables]
data_final = cbind(tab_red[1:11],dados3)








##15.EXPORT THE RESULTS:
#save data:
write.table(data_final, "morpho_clean_log.csv", sep = ",", quote = F, row.names = F)

#create the morphological input by clades
names(data_final)
Individual = data_final$sample #ID for specimens
Species = data_final$clade #Code for putative species, should be the same for genetic data
input = cbind(Individual, Species, data_final[,12:22]) # columns with morphometric data

#save
write.table(input, file = "IBPP_MORPHO_CLADES.txt", sep="\t", quote = FALSE, row.names = FALSE)





#input for populations
Individual = data_final$sample #ID for specimens
Species = data_final$population #Code for putative species, should be the same for genetic data
input = cbind(Individual, Species, data_final[,12:22]) # columns with morphometric data

#save
write.table(input, file = "IBPP_MORPHO_POP.txt", sep="\t", quote = FALSE, row.names = FALSE)

#END;
