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