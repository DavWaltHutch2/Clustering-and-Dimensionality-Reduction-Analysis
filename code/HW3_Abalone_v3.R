#############################################################################################
###################################   ABALONE ANALYSIS ######################################
#############################################################################################


##############################################
###########  SET-UP ENVIRONMENT ##############
##############################################

##CLEAR WORK ENVIRONMENT
rm(list = ls())


##LOAD PACKAGES
##library("EMCluster")
##library("mclust")
##library('plot3D')
##library("RandPro")
##library("e1071")
##library("ica")
##library("fastICA")
##library("caret")
##library("matlib")
##library("cluster")
##library("fpc")
##library("clValid")
##library("infotheo")
##library("corrplot")
if(!require("EMCluster")) install.packages("EMCluster"); library("EMCluster")
if(!require("mclust")) install.packages("mclust"); library("mclust")
if(!require("plot3D")) install.packages("plot3D"); library("plot3D")
if(!require("RandPro")) install.packages("RandPro"); library("RandPro")
if(!require("e1071")) install.packages("e1071"); library("e1071")
if(!require("ica")) install.packages("ica"); library("ica")
if(!require("fastICA")) install.packages("fastICA"); library("fastICA")
if(!require("caret")) install.packages("caret"); library("caret")
if(!require("matlib")) install.packages("matlib"); library("matlib")
if(!require("cluster")) install.packages("cluster"); library("cluster")
if(!require("fpc")) install.packages("fpc"); library("fpc")
if(!require("clValid")) install.packages("clValid"); library("clValid")
if(!require("infotheo")) install.packages("infotheo"); library("infotheo")
if(!require("corrplot")) install.packages("corrplot"); library("corrplot")



##GET FILE
file <- "./data/Abalone_Dataset.csv"
data <- read.csv(file, header = FALSE)
colnames(data) <- c("Sex","Length","Diameter","Height","Whole Weight", 
                    "Shucked Weight", "Viscera Weight", "Shell Weight","Rings")

##ENCODE DATA
features <- dummyVars(" ~ Sex", data = data)
features.encoded <- data.frame(predict(features, newdata = data))
data.encoded <- cbind(data[,!(colnames(data) %in% c("Sex"))],features.encoded)
data <- data.encoded
str(data)

##SPLIT DATA
data.x <- data[,!(colnames(data) %in% c("Rings"))]
data.y <- data$Rings
data.org  <- data


##RECLASSIFY RINGS - SMALL, MODERATE, LARGE
data.z <- data.y
data.z <- as.integer(data.z)
barplot(table(data.z))
data.z[data.z %in% c(1:7)] <- 1
data.z[data.z %in% c(8:11)] <- 2
data.z[data.z %in% c(12:29)] <- 3
barplot(table(data.z))
data.z <- as.factor(data.z)
temp <- data.y
data.y <- data.z
data.z <- temp
table(data.y)
table(data.z)

#########################################################
############### GLOBAL FUNCTIONS ########################
#########################################################
homogeneity_score <- function(labels_true, labels_pred)
{
  C <- labels_true
  K <- labels_pred
  
  return( 1 - condentropy(C,K)/entropy(C))
}

completeness_score <- function(labels_true, labels_pred)
{
  C <- labels_true
  K <- labels_pred
  
  return( 1 - condentropy(K,C)/entropy(K))
}

v_measure <- function(labels_true, labels_pred)
{
  
  retVal <- 2 * homogeneity_score(labels_true, labels_pred)*completeness_score(labels_true, labels_pred)/
    (homogeneity_score(labels_true, labels_pred) + completeness_score(labels_true, labels_pred))
  return(retVal)
}

##########################################################
###############  PART 1: CLUSTERING  #####################
##########################################################

## K-MEANS
##########################################################

##SET SEED
set.seed(76952)

##SCALE DATA
data.input <- as.data.frame(scale(data.x))

##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.input)^2
for (i in 2:15)
{
  model.kmeans <- kmeans(data.input,centers=i,nstart = 10, iter.max = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]
  
plot(2:15, wss, type="b", xlab="Number of Clusters",ylab="Within Group SSE", main = "K-Means Elbow Graph") 
plot(2:15, sil, type="b", xlab = "Number of Clusters", ylab = "Average Silhouette Score", main = "K-Means Silhouette Score")


##CREATE MODEL
centers <- 7 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.input,centers = centers, nstart = 1, iter.max = 20)
end_time <- Sys.time()
clusters.kmeans <- model.kmeans$cluster
total_time <- end_time - start_time

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))
h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)

##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time, sep = ""))


##PLOT 3D - KMEANS
colvar <- as.numeric(clusters.kmeans)
scatter3D(data.input[,7], data.input[,5],data.input[,6], colvar = colvar, theta = 0, phi = 20, 
          main = "Plot of K-Means Clusters")

##PLOT 3D - ORIGINAL
colvar <- as.numeric(data.y)
scatter3D(data.input[,7], data.input[,5],data.input[,6], colvar = colvar, theta = 0, phi = 20, 
          main = "Plot of Original Clusters")









## Expectation Maximazation
##########################################################

##SET SEED
set.seed(367747)

##EM CLUSTERING
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
start_time <- Sys.time()
model.em <- Mclust(data.x, modelNames = models, G= 1:30)
end_time <-Sys.time()
total_time <- end_time - start_time
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT


##VALIDATE RESULTS
clusters.em  <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))
h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)

print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time, sep = ""))


##PLOT 3D - EXPECTATION MAXIMIZATION
colvar <- as.numeric(clusters.em)
scatter3D(data.input[,7], data.input[,5],data.input[,6], col = c("dark blue","dark red"), theta = 0, phi = 20, 
          main = "Plot of Expectation Maximization Clusters")

##PLOT 3D - ORIGINAL
colvar <- as.numeric(data.y)
scatter3D(data.input[,7], data.input[,5],data.input[,6], colvar = colvar, theta = 0, phi = 20, 
          main = "Plot of Original Clusters")



##########################################################
##########  PART 2: DIMENSIONALITY REDUCTION  ############
##########################################################


## PRINCIPAL COMPONENT ANALYSIS (PCA)
##########################################################

##SET SEED
set.seed(387647)

##PEROFRM PCA
model.pca <- prcomp(data.x, center = TRUE, scale = TRUE)

summary(model.pca)

##PLOT CUMMULATIVE IMPORTANCE
cummulative.importance <- summary(model.pca)$importance[3,]
plot(cummulative.importance, type = "b", xlab = "Principal Components", ylab = "Cumulative Importance",
     main = "Cumulative Importance")
 
##PLOT EIGEN VALUES
barplot(model.pca$sdev^2, xlab = "Principal Component", ylab = "Eigenvalue", main = "Eigenvalue Distribution")

##PLOT 3D
colvar <- as.numeric(data.y)
col <- as.character(data.y)
col[col == "1"] <- "blue"
col[col == "2"] <- "green"
col[col == "3"] <- "red"

scatter3D(model.pca$x[,1], model.pca$x[,2], model.pca$x[,3], colvar = colvar, theta = 0, phi = 20, 
          main = "Plot of Best Principle Components (PC1, PC2, PC3)")

plot(model.pca$x[,1], model.pca$x[,2], col = col)

##GET DATA
data.pca <- model.pca$x[,1:2]
str(data.pca)



## INDEPENDENT COMPONENT ANALYSIS (ICA)
##########################################################

##EXPLORATORY ANALYSIS
corr <- cor(data.x)
corrplot(corr, method = "circle", tl.cex = .75)


##SET SEED
set.seed(4766)

##CHOOSE NUMBER OF COMPONENTS
accuracy.all <- NULL
for(i in 2:10)
{
  print(paste("Components: ", i))
  ica.results <- fastICA(data.x, n.comp = i)
  S <- as.data.frame(ica.results$S)
  S <- cbind(S,"diagnosis" = data.y)
  
  model.svm <- svm(diagnosis~., data = S, kernel = "linear")
  pred <- predict(model.svm, S)
  cf <- confusionMatrix(pred, data.y)
  accuracy.all <- c(accuracy.all, round(cf$overall["Accuracy"]*100, 2))
  
}
plot(accuracy.all, type = "b", xlab = "Total Independent Components", ylab ="SVM Accuracy", 
     main = "SVM Accuracy vs Total Independent Components")

##PLOT KURTOSIS DISTRIBUTION
n.comp <- 4
ica.model <- fastICA(data.x, n.comp = n.comp)
S <- ica.model$S

kurtosis.all = NULL
for(i in 1:n.comp)
{
  print(kurtosis(S[,i]))
  kurtosis.all <- c(kurtosis.all, abs(kurtosis(S[,i])))
}

order <- order(kurtosis.all, decreasing = TRUE)
barplot(kurtosis.all[order], names.arg = order, ylab = "Excess Kurtosis", xlab = "Independent Components", 
        main = "Ordered Independent Components")


##PLOT KURTOSIS
colvar <- as.numeric(data.y)
col <- as.character(data.y)
col[col == "1"] <- "blue"
col[col == "2"] <- "green"
col[col == "3"] <- "red"

plot(S[,3], col = col, xlab = "Index", ylab = "IC 3", main = "High Kurtosis Independent Component (IC3)")

scatter3D(S[,2], S[,1], S[,4], colvar = colvar, theta = 0, phi = 20, 
          main = "Low Kurtosis Independent Components (IC1, IC2, IC4)") ##HIGH KURTOSIS

plot(S[,1], col = col, xlab = "Index", ylab = "Independent Component 1", main = "Independent Component 1")




##GET DATA 
data.ica <- ica.model$S
str(data.ica)

data.ica.all <- data.ica
str(data.ica.all)
data.ica.high <- data.ica[,3]
str(data.ica.high)
data.ica.low <- data.ica[, c(2,1,4)]
str(data.ica.low)


## RANDOMIZED PROJECTIONS
##########################################################
##SET SEED
set.seed(47676446)

accuracy.all <- NULL
rows <- dim(data.x)[2]
cols <- dim(data.x)[2]
trials <- 25
for( i in 1:cols)
{
  print(paste("Columns: ", i))
  accuracy.avg <- NULL
  for(j in 1:trials)
  {
    proj.matrix <- form_matrix(rows,i,FALSE)
    data.transformed <- as.matrix(data.x) %*% proj.matrix
    data.transformed <- as.data.frame(data.transformed)
    data.transformed <- cbind(data.transformed, "diagnosis" = data.y)
    
    model.svm <- svm(diagnosis~., data = data.transformed, kernel = "linear")
    pred <- predict(model.svm, data.transformed)
    cf <- confusionMatrix(pred, data.y)
    accuracy.avg <- c(accuracy.avg, round(cf$overall["Accuracy"]*100, 2))
  }

  accuracy.all <- c(accuracy.all, mean(accuracy.avg))
}
plot(accuracy.all, type = "b", xlab = "Total Random Projections",ylab = "SVM Accuracy", 
     main = "SVM Accuracy vs Total Random Projections")


##MAP TO LOWER DIMENSION
rows <- dim(data.x)[2]
cols <- 6
proj.matrix <- form_matrix(rows,cols,FALSE, projection = "gaussian")
data.transformed <- as.matrix(data.x) %*% proj.matrix
data.transformed <- as.data.frame(data.transformed)
##data.transformed <- cbind(data.transformed, "diagnosis" = data.y)


##PLOT 3D
colvar <- as.numeric(data.y)
col <- as.character(data.y)
col[col == "1"] <- "blue"
col[col == "2"] <- "green"
col[col == "3"] <- "red"
scatter3D(data.transformed[,1], data.transformed[,2], data.transformed[,3], colvar = colvar, 
          theta = 90, phi = 20, main = "Plot Over Random Projections") 
plot(data.transformed[,1], data.transformed[,2], col = col)

##GET DATA
data.random <- data.transformed
str(data.random)



## RFE
##########################################################

##SET SEED
set.seed(6564)

##DEFINE IMPORTANT FEATURES
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results <- rfe(data.x, data.y, rfeControl=control)
print(results)
plot(results, type=c("g", "o"), main = "Recursive Feature Elimination")

##GET DATA
predictors <- predictors(results)

data.rfe <- data[,predictors]
str(data.rfe)




##########################################################
####  PART 3: DIMENSIONALITY REDUCTION AND CLUSTERING  ###
##########################################################

##DATA
##########################################################
#data.pca
#data.ica
#data.random
#data.rfe


## K-MEANS
##########################################################

##data.pca
#############################
data.reduced <- data.pca
str(data.reduced)

##SET SEED
set.seed(56347)

##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.reduced)^2
max.clusters <- 20
for (i in 2:max.clusters)
{
  model.kmeans <- kmeans(data.reduced,centers=i,nstart = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]

plot(2:max.clusters, wss, type="b", xlab="Number of Clusters",ylab="Within Group SSE", main = "Elbow Graph") 
plot(2:max.clusters, sil, type="b", xlab = "Number of Clusters", ylab = "Average Silhouette Score", main = "Silhouette Score")


##CREATE MODEL
centers <- 5 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.reduced,centers = centers, nstart = 1)
end_time <- Sys.time()
clusters.kmeans <- model.kmeans$cluster
total_time <- end_time - start_time

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))

h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)

##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))








##data.ica (HIGH KURTOSIS COMPONENTS)
#############################
data.reduced <- data.ica.high
str(data.reduced)

##SET SEED
set.seed(9877)

##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.reduced)^2
for (i in 2:15)
{
  model.kmeans <- kmeans(data.reduced,iter.max = 25, centers=i,nstart = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]

plot(2:15, wss, type="b", xlab="Number of Clusters",ylab="Within Group SSE", main = "Elbow Graph") 
plot(2:15, sil, type="b", xlab = "Number of Clusters", ylab = "Average Silhouette Score", main = "Silhouette Score")


##CREATE MODEL
centers <- 4 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.reduced,centers = centers, nstart = 1)
end_time <- Sys.time()
total_time <- end_time - start_time
clusters.kmeans <- model.kmeans$cluster

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))

h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)

##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))















##data.ica (LOW KURTOSIS COMPONENTS)
#############################
data.reduced <- data.ica.low
str(data.reduced)

##SET SEED
set.seed(9877)

##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.reduced)^2
for (i in 2:15)
{
  model.kmeans <- kmeans(data.reduced,centers=i,nstart = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]

plot(2:15, wss, type="b", xlab="Number of Clusters",ylab="Within Group SSE", main = "Elbow Graph") 
plot(2:15, sil, type="b", xlab = "Number of Clusters", ylab = "Average Silhouette Score", main = "Silhouette Score")


##CREATE MODEL
centers <- 9 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.reduced,centers = centers, nstart = 1)
end_time <- Sys.time()
total_time <- end_time - start_time
clusters.kmeans <- model.kmeans$cluster

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))

h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)

##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))












##data.random
#############################
data.reduced <- data.random
str(data.reduced)

##SET SEED
set.seed(56837)

##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.reduced)^2
for (i in 2:15)
{
  model.kmeans <- kmeans(data.reduced,centers=i,nstart = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]

plot(2:15, wss, type="b", xlab="Number of Clusterss",ylab="Within Group SSE", main = "Elbow Graph") 
plot(2:15, sil, type="b", xlab = "Number of Clusterss", ylab = "Average Silhouette Score", main = "Silhouette Score")


##CREATE MODEL
centers <- 6 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.reduced,centers = centers, nstart = 1)
end_time <- Sys.time()
total_time <- end_time - start_time
clusters.kmeans <- model.kmeans$cluster

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))

h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)


##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))






##data.rfe
#############################
data.reduced <- data.rfe
data.reduced <- as.data.frame(scale(data.reduced))

##SET SEED
set.seed(56347)


##DEFINE NUMBER OF CLUSTERS
wss <- NULL
sil <- NULL
distance <- dist(data.reduced)^2
for (i in 2:15)
{
  model.kmeans <- kmeans(data.reduced,centers=i,nstart = 20)
  wss[i] <- sum(model.kmeans$withinss)
  sil[i] <- summary(silhouette(model.kmeans$cluster,distance))$avg.width
}

wss <- wss[!(is.na(wss))]
sil <- sil[!(is.na(sil))]

plot(2:15, wss, type="b", xlab="Number of Clusterss",ylab="Within Group SSE", main = "Elbow Graph") 
plot(2:15, sil, type="b", xlab = "Number of Clusterss", ylab = "Average Silhouette Score", main = "Silhouette Score")


##CREATE MODEL
centers <- 6 ##IDEAL NUMBER OF CENTERS
start_time <- Sys.time()
model.kmeans <- kmeans(data.reduced,centers = centers, nstart = 1)
end_time <- Sys.time()
total_time <- end_time - start_time
clusters.kmeans <- model.kmeans$cluster

##VALIDATE RESULTS
(confusion.matrix <- table(clusters.kmeans, data.y))

h <- homogeneity_score(data.y, clusters.kmeans)
c <- completeness_score(data.y, clusters.kmeans)
v <- v_measure(data.y, clusters.kmeans)


##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))




## EM
##########################################################


##data.pca
#############################

##SET SEED
set.seed(367747)

##GET MODEL
data.reduced <- data.pca
model.em <- Mclust(data.reduced, G = c(2:7))
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT

##GET TIME
start_time <- Sys.time()
model.time <- Mclust(data.reduced, G = model.em$G, modelNames = model.em$modelName)
end_time <- Sys.time()
total_time <- end_time - start_time


##VALIDATE RESULTS
clusters.em <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))

h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)

##PRINT RRESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))



##data.ica.high
#############################

##SET SEED
set.seed(367747)

##GET MODEL
data.reduced <- data.ica.high
str(data.reduced)
model.em <- Mclust(data.reduced)
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT


##GET TIME
start_time <- Sys.time()
model.time <- Mclust(data.reduced, G = model.em$G, modelNames = model.em$modelName)
end_time <- Sys.time()
total_time <- end_time - start_time


##VALIDATE RESULTS
clusters.em <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))

h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)

#PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))



##data.ica.lowK
#############################

##SET SEED
set.seed(367747)

##GET MODEL
data.reduced <- as.data.frame(data.ica.low)
str(data.reduced)
model.em <- Mclust(data.reduced)
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT

##GET TIME
start_time <- Sys.time()
model.time <- Mclust(data.reduced, G = model.em$G, modelNames = model.em$modelName)
end_time <- Sys.time()
total_time <- end_time - start_time


##VALIDATE RESULTS
clusters.em <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))

h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)


#PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))







##data.random
#############################

##SET SEED
set.seed(367747)

##GET MODEL
data.reduced <- data.random
model.em <- Mclust(data.reduced, modelNames = c("VVV"), G = 3:7)
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT

##GET TIME
start_time <- Sys.time()
model.time <- Mclust(data.reduced, G = model.em$G, modelNames = model.em$modelName)
end_time <- Sys.time()
total_time <- end_time - start_time


##VALIDATE RESULTS
clusters.em <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))

h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)

##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))


##data.rfe
#############################

##SET SEED
set.seed(367747)

##GET MODEL
data.reduced <- data.rfe
str(data.reduced)
model.em <- Mclust(data.reduced)
summary(model.em)
##plot(model.em) ##UNCOMMENT TO PLOT

##GET TIME
start_time <- Sys.time()
model.time <- Mclust(data.reduced, G = model.em$G, modelNames = model.em$modelName)
end_time <- Sys.time()
total_time <- end_time - start_time


##VALIDATE RESULTS
clusters.em <- model.em$classification
(confusion.matrix <- table(clusters.em, data.y))

h <- homogeneity_score(data.y, clusters.em)
c <- completeness_score(data.y, clusters.em)
v <- v_measure(data.y, clusters.em)


##PRINT RESULTS
print(paste("Homogeneity Score: ", h, sep = ""))
print(paste("Completeness Score: ", c, sep = ""))
print(paste("V Measure: ", v, sep = ""))
print(paste("Time: ", total_time))




##########################################################
#######  PART 4: DIMENSIONALITY REDUCTION AND ANN  #######
##########################################################


##DATA
##########################################################
#data.pca
#data.ica
#data.random
#data.rfe





## BASELINE
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(data.x, "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)




## PRINCIPLE COMPONENT ANALYSIS (PCA)
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(as.data.frame(data.pca), "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)









## INDEPENDENT COMPONENT ANALYSIS (ICA - High Kurtosis)
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(as.data.frame(data.ica.high), "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)









## INDEPENDENT COMPONENT ANALYSIS (ICA - Low Kurtosis)
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(as.data.frame(data.ica.low), "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)








## RANDOM PROJECCTIONS (RP)
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(as.data.frame(data.random), "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)









## RECURSIVE FEATURE ELIMINATION (RFE)
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
#############################
data.reduced <- cbind(as.data.frame(data.rfe), "diagnosis" = data.y)
data.reduced <- cbind(as.data.frame(scale(data.rfe)), "diagnosis" = data.y)
str(data.reduced)


##CREATE TRAIN AND TEST DATA 
train.index <- createDataPartition(data.reduced$diagnosis, p = .7, list = FALSE)
data.train <- data.reduced[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.reduced[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

##CREATE AND TEST BASELINE MODEL
ctrl <- trainControl(method="cv",number = 10) 
size <- 30
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
(time.cv <- model.ann$time$everything["user.self"])


predict.ann <- predict(model.ann, newdata = data.testing)
confusionMatrix(predict.ann, data.testing$diagnosis)








##########################################################
#############  PART 5: CLUSTERING AND ANN  ###############
##########################################################


## KMEANS
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
classification <- as.factor(clusters.kmeans)
data.classification <- cbind(data.x, "classification" = classification)
str(data.classification)


##ENCODE FACTORS
features <- dummyVars(" ~ classification", data = data.classification)
features.encoded <- data.frame(predict(features, newdata = data.classification))
data.encoded <- cbind(data.classification[,!(colnames(data.classification) %in% c("classification"))],features.encoded)
str(data.encoded)


##RELABEL DATA
data.kmeans <- cbind(as.data.frame(data.encoded), "diagnosis" = data.y)
str(data.kmeans)

##CREATE DATA
data.kmeans.classifications <- data.kmeans[,colnames(data.kmeans) %in% c("classification.1", "classification.2", "classification.3", "diagnosis")]  
str(data.kmeans.classifications)

data.kmeans.all <- data.kmeans  
str(data.kmeans.all)

##CREATE PARTITIONS 
train.index <- createDataPartition(data.kmeans$diagnosis, p = .7, list = FALSE)

##CREATE MODEL - ALL 
####################
data.kmeans <- data.kmeans
str(data.kmeans)

data.train <- data.kmeans[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.kmeans[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

ctrl <- trainControl(method="cv",number = 10) 
size <- c(5, 10, 15, 20, 25, 30, 35, 40)
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))

##TEST MODEL
predict.ann <- predict(model.ann, newdata = data.testing)
cf <- confusionMatrix(predict.ann, data.testing$diagnosis)

##TIME MODEL
ctrl <- trainControl(method="cv",number = 10) 
(size <- summary(model.ann)$tuneValue$size)
(decay <- summary(model.ann)$tuneValue$decay)
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
start_time <- Sys.time()
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
end_time <- Sys.time()
total_time <- end_time - start_time

##RESULTS
print(cf)
print(paste("Time: ", total_time, sep = ""))



##CREATE MODEL - CLASSIFICATIONS
################################
data.kmeans <- data.kmeans.classifications
str(data.kmeans)

data.train <- data.kmeans[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.kmeans[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

ctrl <- trainControl(method="cv",number = 10) 
size <- c(5, 10, 15, 20, 25, 30, 35, 40)
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
start_time <- Sys.time()
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
end_time <- Sys.time()
total_time <- end_time - start_time

##TEST MODEL
predict.ann <- predict(model.ann, newdata = data.testing)
cf <- confusionMatrix(predict.ann, data.testing$diagnosis)

##TIME MODEL
ctrl <- trainControl(method="cv",number = 10) 
(size <- summary(model.ann)$tuneValue$size)
(decay <- summary(model.ann)$tuneValue$decay)
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
start_time <- Sys.time()
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
end_time <- Sys.time()
total_time <- end_time - start_time

##RESULTS
print(cf)
print(paste("Time: ", total_time, sep = ""))







## EXPECTATION MAXIMIZATIONS
##########################################################

##SET SEED
set.seed(56347)

##GET DATA
classification <- as.factor(clusters.em)
data.classification <- cbind(data.x, "classification" = classification)
str(data.classification)


##ENCODE FACTORS
data.classification$classification <- as.integer(data.classification$classification)
data.encoded <- data.classification
str(data.encoded)

##RELABEL DATA
data.em <- cbind(as.data.frame(data.encoded), "diagnosis" = data.y)
str(data.em)

##CREATE DATA
data.em.classifications <- data.em[,colnames(data.em) %in% c("classification","diagnosis")]  
str(data.em.classifications)

data.em.all <- data.em  
str(data.em.all)

##CREATE PARTITIONS 
train.index <- createDataPartition(data.em$diagnosis, p = .7, list = FALSE)

##CREATE MODEL - ALL 
####################
data.em <- data.em.all
str(data.em)

data.train <- data.em[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.em[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

ctrl <- trainControl(method="cv",number = 10) 
size <- c(5, 10, 15, 20, 25, 30, 35, 40)
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))


##TEST MODEL
predict.ann <- predict(model.ann, newdata = data.testing)
cf <- confusionMatrix(predict.ann, data.testing$diagnosis)

##TIME MODEL
ctrl <- trainControl(method="cv",number = 10) 
(size <- summary(model.ann)$tuneValue$size)
(decay <- summary(model.ann)$tuneValue$decay)
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
start_time <- Sys.time()
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
end_time <- Sys.time()
total_time <- end_time - start_time

##RESULTS
print(cf)
print(paste("Time: ", total_time, sep = ""))



##CREATE MODEL - CLASSIFICATIONS
################################
data.em <- data.em.classifications
str(data.em)

data.train <- data.em[train.index,]
table(data.train$diagnosis)
prop.table(table(data.train$diagnosis))
data.training <- data.train
str(data.training)

data.test <- data.em[-train.index,]
table(data.test$diagnosis)
prop.table(table(data.test$diagnosis))
data.testing <- data.test
str(data.testing)

ctrl <- trainControl(method="cv",number = 10) 
size <- c(5, 10, 15, 20, 25, 30, 35, 40)
decay <- 0
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))


##TEST MODEL
predict.ann <- predict(model.ann, newdata = data.testing)
cf <- confusionMatrix(predict.ann, data.testing$diagnosis)

##TIME MODEL
ctrl <- trainControl(method="cv",number = 10) 
(size <- summary(model.ann)$tuneValue$size)
(decay <- summary(model.ann)$tuneValue$decay)
tuneGrid = expand.grid(size = size, decay = decay) ##modelLookup("gbm") ##names(getModelInfo())
start_time <- Sys.time()
(model.ann <- train(diagnosis ~ ., data.training, method = "nnet", verbose = FALSE, trControl = ctrl, tuneGrid = tuneGrid, metric = "Accuracy", maxit = 300))
end_time <- Sys.time()
total_time <- end_time - start_time

##RESULTS
print(cf)
print(paste("Time: ", total_time, sep = ""))


