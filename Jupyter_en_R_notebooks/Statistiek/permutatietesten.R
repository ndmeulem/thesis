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

#number of groups = 9 (9 patients)
length(features2[,"patientnumber"]==122)
H = subset(features2, patientnumber == 122)
length(H[,"treatment"]=="CXCL1")
G = combn(10, 5)
dim(G)

G = combn(9792, 4)
dim(G)

groepPerm <- sample(groep)    #hier moet dan treatment per patiënt nog komen
fStar=sapply(X=1:B, FUN=function(b,y,groep)
{anova(lm(y~groepPerm)$F[1])},y=log(dna$length),  #hier moet ~groepperm + patient nog komen
groep=dna$dose)

sample.df <- function(df) df[sample(df[,"treatment"])]

features2_permutated = sample.df(features2_relevelled[features2_relevelled$patientnumber==122,])
features2_permutated

test = sample(features2[,"treatment"])
features2_test = features2
features2_test$treatment=test

H = subset(features2, patientnumber == 122)
length(H[,"treatment"])
perm_122 = sample(H[,'treatment'])
features2_perm_122 =  H
features2_perm_122$treatment = perm_122

I = subset(features2, patientnumber == 129)
length(I[,'treatment'])
perm_129 = sample(I[,'treatment'])
features2_perm_129 =  I
features2_perm_129$treatment = perm_129

J = subset(features2, patientnumber == 164)
perm_164 = sample(J[,'treatment'])
features2_perm_164 =  J
features2_perm_164$treatment = perm_164

K = subset(features2, patientnumber == 722)
perm_722 = sample(K[,'treatment'])
features2_perm_722 =  K
features2_perm_722$treatment = perm_722

L = subset(features2, patientnumber == 714)
perm_714 = sample(L[,'treatment'])
features2_perm_714 =  L
features2_perm_714$treatment = perm_714

M = subset(features2, patientnumber == 751)
perm_751 = sample(M[,'treatment'])
features2_perm_751 =  M
features2_perm_751$treatment = perm_751

N = subset(features2, patientnumber == 805)
perm_805 = sample(N[,'treatment'])
features2_perm_805 = N
features2_perm_805$treatment = perm_805

O = subset(features2, patientnumber == 858)
perm_858 = sample(O[,'treatment'])
features2_perm_858 =  O
features2_perm_858$treatment = perm_858

P = subset(features2, patientnumber == 859)
perm_859 = sample(P[,'treatment'])
features2_perm_859 =  P
features2_perm_859$treatment = perm_859

features2_permutated = rbind(features2_perm_122, features2_perm_129, features2_perm_164, features2_perm_714,
                             features2_perm_722, features2_perm_751, features2_perm_805, features2_perm_858, features2_perm_859)

set.seed(250)
B=10000
f = lm(accumulated_distance~treatment+patientnumber,data=features2_relevelled)
fOrig = aov(f)
summary(fOrig)
fStar=sapply(X=1:B, FUN=function(b,y,treatment)
{anova(lm(y~features2_permutated$treatment+features2_permutated$patientnumber))$F[1]},y=(features2_relevelled$accumulated_distance))
fStar
summary(fStar)

pval2=(sum(fStar>=12.19)+1)/(B+1)
pval2

condition = features2_perm_122[,'treatment'] == "CXCL1"
test = features2_perm_122[condition,]
hist(test$accumulated_distance)

write.csv(features2_permutated, file = "permutated_dataframe_features2")

condition = features2[,"treatment"] == "CXCL1"
test = features2[condition,]
hist(test$accumulated_distance)


#maak een for loop om per patientnummer een gerandomiseerde data subset te maken (maw, kort regel 36 tot 83 in)
#doe dan anova voor die ene permutatie zonder de apply en B en zo
#probeer dan die for loop om te zetten in een apply functie en doe dan 10000 iteraties om random samples te maken (voeg B terug toe)


set.seed(165)
B=10000
fOrig = anova(lm(accumulated_distance~treatment+patientnumber,data=features2_relevelled))$F[1]
fStar=sapply(X=1:B, FUN=function(b,y,groep) {anova(lm(y~sample(groep)))$F[1]},
             y=features2_relevelled$accumulated_distance,groep=features2_relevelled$treatment) 
pval=(sum(fStar>=fOrig)+1)/(B+1)
pval


#Follow tutorial from Vermont University: Manly's approach -- Unrestricted Permutation of Observations

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]

print( "Resampling as in Manly with unrestricted sampling of observations. ")
# Now start resampling
nreps <- 5000 
FT <- numeric(nreps)    #Set up space to store F values as calculated.
FP <- numeric(nreps)  
FT[1] <- Ftreatment        # The first F of our 5000 
FP[1] <- Fpatientnumber
for (i in 2:nreps) {
  newtreatment <- sample(features2_relevelled$treatment, 9792)
  mod2 <- lm(features2_relevelled$accumulated_distance ~ newtreatment + features2_relevelled$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}
probT <- length(FT[FT >= Ftreatment + .Machine$double.eps ^0.5])/nreps
probP <- length(FP[FP >= Fpatientnumber+ .Machine$double.eps ^0.5])/nreps       
### The addition of "+ .Machine$double.eps" is an aid against two numbers that differ only by
### floating point computer calculations at the extreme.

cat(" The probability value for treatment is ", probT, "\n") 
cat(" The probability value for patientnumber is ", probP, "\n")
lijst = list(FT >= Ftreatment + .Machine$double.eps ^0.5)
is.element(TRUE, lijst)









#resampling nested designs
treat <- as.factor(rep(1:36, each = 10))
gender <- as.factor(rep(1:9, each = 50))
n <- 10
N <- length(dv)
g <- nlevels(gender)
tt <- nlevels(ther)

# F values from Nested.R
mod1 <- summary(aov(dv ~ gender*ther))
Fthwingen1 <- mod1[[1]]$"F value"[2]
Fgen1 <- mod1[[1]]$"Mean Sq"[1]/mod1[[1]]$"Mean Sq"[2]
Fthwingen <- numeric(nreps)
Fgen <- numeric(nreps)
# First we will compute the Gender effect
#  Randomly shuffle Therapists across Gender
#    Make a matrix containing therapists on columns
dvmat <- matrix(dv, nrow = tt, byrow = TRUE)

for (i in 1:nreps) {
  worms <- dvmat[sample(1:10,10),]   #worms???
  newdv <- as.vector(t(worms))
  a <- summary(aov(newdv~gender*ther))
  # This gives us the appropriate terms because it removes
  # the 1 df that would go to the interaction from th, giving ther(gender)
  # But if you enter the terms in the other order it won't work.
  MSgen <- a[[1]]$"Mean Sq"[1]
  MSth <- a[[1]]$"Mean Sq"[2]
  MSerrorwin <-  a[[1]]$"Mean Sq"[3]
  Fgen[i] <- MSgen/MSth
  
  S1 <- dv[gender == 1]
  S2 <- dv[gender == 2]
  temp1 <- sample(S1, length(S1), replace = FALSE)
  temp2 <- sample(S2, length(S2), replace = FALSE)
  newdv2 <- c(temp1, temp2)
  mod6<- summary(aov(newdv2 ~ gender * ther))
  Fthwingen[i] <- mod6[[1]]$"F value"[2]
  #Fgen[i] <- MSgen/  mod6[[1]]$"Mean Sq"[2]
}

probgen <- length(Fgen[Fgen >=  Fgen1])/nreps
probth <- length(Fthwingen[Fthwingen >= Fthwingen1])/nreps
par(mfrow = c(2,2))
hist(Fthwingen, breaks = 50)
hist(Fgen, breaks = 50)
# This works, but the question is what is the denominator to test gender.
# It should be MS(ther(gender)), but that is calculated with a different randomization.
# I have used MSther(gender) from first analysis.




















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
  mod2 <- lm(permutated_dataset$accumulated_distance ~ permutated_dataset$treatment + permutated_dataset$patientnumber)
  b <- summary(aov(mod2))
  FT[i] <- b[[1]]$"F value"[1]
  FP[i] <- b[[1]]$"F value"[2]
}

# Standard Anova on these data
mod1 <- lm(accumulated_distance ~ treatment + patientnumber, data = features2_relevelled)
ANOVA <- summary(aov(mod1))
cat( " The standard ANOVA for these data follows ","\n")
Ftreatment <-  ANOVA[[1]]$"F value"[1]   # Saving F values for future use
Fpatientnumber <-  ANOVA[[1]]$"F value"[2]


pval=(sum(FT>=Ftreatment)+1)/(nreps+1)
pval






probT <- length(FT[FT >= Ftreatment + .Machine$double.eps ^0.5])/nreps
probP <- length(FP[FP >= Fpatientnumber+ .Machine$double.eps ^0.5])/nreps       
### The addition of "+ .Machine$double.eps" is an aid against two numbers that differ only by
### floating point computer calculations at the extreme.

cat(" The probability value for treatment is ", probT, "\n") 
cat(" The probability value for patientnumber is ", probP, "\n")







set.seed(165)
B=10000
fOrig = anova(lm(accumulated_distance~treatment+patientnumber,data=features2_relevelled))$F[1]
fStar=sapply(X=1:B, FUN=function(b,y,treatment) {anova(lm(y~sample(groep)))$F[1]},
             y=features2_relevelled$accumulated_distance,groep=features2_relevelled$treatment) 
pval=(sum(fStar>=fOrig)+1)/(B+1)
pval

test = sample(c('test', 'random', 'stress', 'anxiety'))

test = double()