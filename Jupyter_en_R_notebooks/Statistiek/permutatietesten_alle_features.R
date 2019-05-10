getwd()
setwd("C:/Users/ninad/Documents/Thesis/Jupyter_en_R_notebooks/Statistiek")

features2 = read.csv(file = "Dataset_features_2", header = TRUE, fill = TRUE)

#install.packages("plyr")
library(plyr)
features2 = rename(features2, c("accumulated_euclidian_distance"="euclidian_distance"))

#install.packages("car")
library(car)

#set PBS as baseline for the model to compare with
features2_relevelled <- within(features2, treatment <- relevel(treatment, ref = "PBS"))

#write function to sample treatment within each patient
patients  =c(122, 129, 164, 714, 722, 751, 805, 858, 859)

  #shuffle treatments within patient
shuffle = function(dataset, number) {
  A = subset(dataset, patientnumber == number)
  newtreatment = sample(A[,'treatment'])
  perm_data = A
  perm_data$treatment = newtreatment
  perm_data
}
test = shuffle(features2, 722)

shuffle4 = function(dataset, number, treatment1, treatment2) {
  A = subset(dataset, patientnumber == number)
  B = subset(A, treatment == treatment1 | treatment == treatment2)
  C = subset(A, treatment != treatment1 & treatment != treatment2)
  newtreatment = sample(B[,'treatment'])
  B$treatment = newtreatment
  D = rbind(B, C)
  D
}
B = shuffle4(features2, 122, "CXCL1", "PBS")
#1.accumulated distance

#for loop in R is niet zo efficiënt want hij stored alles constant in het RAM geheugen en als de dataset te groot wordt dan switchet hij van 
#plaats, maar die plaatsen worden niet bijgehouden, eens het RAM geheugen dan vol zit, dan 
#run function over data, for number of times you want to do iteration
set.seed(1)
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}
b

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)

treatment = permutated_dataset$treatment
#post hoc test met multcomp library
if(!require(multcomp)){install.packages("multcomp")}
library(multcomp)
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)






##post hoc test
#paarsgewijs voor de treatments checken of er een verschil is (6 mogelijke combinaties : binomiaalcoëfficiënt (4   2))

#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (extreem)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> niet significant (p = 0.074)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##niet significant (p  = 0.339)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.096)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#significant (p = 0.0052)








#2 euclidian distance

#run function over data, for number of times you want to do iteration
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval          #0.00019996
#post hoc
library(multcomp)
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)








#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (extreem, p = 0.0003998401)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (p = 0.0003998401)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##heel significant (p = 0.0003998401)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$euclidian_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(euclidian_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#significant (p = 0.0003998401)










#3 average directness

#run function over data, for number of times you want to do iteration
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval

sum(FT)
sum(FT>=Ftreatment)

#post hoc test
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)





#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (extreem, p = 0.0003998401)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant? (p = 0.02958816)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##heel significant (p = 0.0003998401)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$average_directness ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(average_directness ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#significant (p = 0.0003998401)










#4mean_speed

#run function over data, for number of times you want to do iteration
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval

sum(FT)
sum(FT>=Ftreatment)

#post hoc test
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)






#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> niet significant (p = 0.6565374)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (p = 0.0003998401)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.05877649)

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##heel significant (p = 0.0003998401)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#significant (p = 0.00239904)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$mean_speed ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(mean_speed ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#significant (p = 0.009196321)












#5 stdev_velocity

#run function over data, for number of times you want to do iteration
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$stdev_velocity ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(stdev_velocity ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval

sum(FT)
sum(FT>=Ftreatment)

#post hoc test
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)
















#6.forward migration index x
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval

sum(FT)
sum(FT>=Ftreatment)

#post hoc test
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)



#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> heel significant (p = 0.0003998401)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> niet significant (p = 0.4578169)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.2279088)

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##heel significant (p = 0.0003998401)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_x ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_x ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.3154738)









#6.forward migration index y
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle(features2, patient)
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval

sum(FT)
sum(FT>=Ftreatment)

#post hoc test
model1.mcp<-glht(mod1,linfct=mcp(treatment="Tukey"))
summary(model1.mcp)

#1. CXCL1 - CXCL8
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "CXCL8")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> niet significant (p = 0.297481)


#2. CXCL1, PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

# -> significant (p = 0.0003998401)

# CXCL1 fMLP

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL1", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.05877649)

#CXCL8 PBS

set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "PBS")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

##heel significant (p = 0.0003998401)

#CXCL8 fMLP
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "CXCL8", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#heel significant (p = 0.0003998401)

#fMLP PBS
set.seed(1)
nreps <- 2500 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps) 
for (i in 2:nreps) {
  permutated_dataset = data.frame()
  for (patient in patients){
    newtreatment_dataset <- shuffle4(features2, patient, "PBS", "fMLP")
    permutated_dataset = rbind(permutated_dataset, newtreatment_dataset)
  }
  mod2 <- lm(permutated_dataset$forward_migration_index_y ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(forward_migration_index_y ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
ANOVA
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval
hist(FT, breaks = 100)
sum(FT>=Ftreatment)

#niet significant (p = 0.5621751)








#andere manier om het te proberen

newdf = data.frame()


patients  =c(122, 129, 164, 714, 722, 751, 805, 858, 859)


shuffle2 = function(dataset, number) {
  A = subset(dataset, patientnumber == number)
  newtreatment = sample(A[,'treatment'])
  perm_data = A
  perm_data$treatment = newtreatment
  perm_data
}
test = shuffle2(features2, 722)

typeof(test)

shuffle3=function(dataset){
  treatments = data.frame()
  for (patient in patients){
    treatments = rbind(treatments, shuffle2(dataset, patient))
  }
  treatments$treatment
}


set.seed(165)
B=10000
fOrig=anova(lm(accumulated_distance~treatment+patientnumber,data=features2_relevelled))$F[1]
fStar=sapply(X=1:B, FUN=function(b,y,groep)
{anova(lm(features2_relevelled$accumulated_distance~shuffle3(features2)+features2_relevelled$patientnumber))$F[1]})
fOrig

pval2=(sum(fStar>=fOrig)+1)/(B+1)
pval2

sum(fStar>=fOrig)


set.seed(150)
B=10000
fOrig=anova(lm(euclidian_distance~treatment+patientnumber,data=features2_relevelled))$F[1]
fStar=sapply(X=1:B, FUN=function(b,y,groep)
{anova(lm(features2_relevelled$euclidian_distance~shuffle3(features2)+features2_relevelled$patientnumber))$F[1]})
fOrig

pval2=(sum(fStar>=fOrig)+1)/(B+1)
pval2

sum(fStar>=fOrig)


set.seed(150)
B=5000
fOrig=anova(lm(mean_speed~treatment+patientnumber,data=features2_relevelled))$F[1]
fStar=sapply(X=1:B, FUN=function(b,y,groep)
{anova(lm(features2_relevelled$mean_speed~shuffle3(features2)+features2_relevelled$patientnumber))$F[1]})
fOrig

pval2=(sum(fStar>=fOrig)+1)/(B+1)
pval2

sum(fStar>=fOrig)

hist(fStar, breaks = 100)