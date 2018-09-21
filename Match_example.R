setwd("~/Dropbox/Research/Code/Matching_Algorithms")
library(MASS)
library(sn)
library(mnormt)
library(nnet)
library(car)
library(Matching)
library(dplyr)
library(cluster)
library(e1071)

source("Match_example_functions.R") #for FM function and MaxMax2SB

# Function to calculate MaxMax2SB for the pre-matched cohort
match5.pre <- function(SD, x1, x2, x3, x4, x5){
  
  max(abs((mean(x1)-mean(x2))/SD), abs((mean(x1)-mean(x3))/SD), abs((mean(x1)-mean(x4))/SD), 
      abs((mean(x1)-mean(x5))/SD), abs((mean(x2)-mean(x3))/SD), abs((mean(x2)-mean(x4))/SD), 
      abs((mean(x2)-mean(x5))/SD), abs((mean(x3)-mean(x4))/SD), abs((mean(x3)-mean(x5))/SD), 
      abs((mean(x4)-mean(x5))/SD))
}

# A Simulated data example for Z=5 treatment groups
# using VM, FM, FMnc, GPSnc, and COVnc

# Some configurations 
n1s <- 1200
g <- 2
b <- 0.50
lambda <- 0
s2 <- 1
s3 <- 1
P <- 10
df <- 7
eta <- 0

mu <- b*sqrt((1+2*s2+2*s3)/5)

mu1 <- as.matrix(c(mu,0,0,0,0,mu,0,0,0,0))
mu2 <- as.matrix(c(0,mu,0,0,0,0,mu,0,0,0))
mu3 <- as.matrix(c(0,0,mu,0,0,0,0,mu,0,0))
mu4 <- as.matrix(c(0,0,0,mu,0,0,0,0,mu,0))
mu5 <- as.matrix(c(0,0,0,0,mu,0,0,0,0,mu))
cov1 <- matrix(lambda, nrow = 10, ncol = 10); diag(cov1) <- 1
cov2 <- matrix(lambda, nrow = 10, ncol = 10); diag(cov2) <- s2
cov3 <- matrix(lambda, nrow = 10, ncol = 10); diag(cov3) <- s3
cov4 <- matrix(lambda, nrow = 10, ncol = 10); diag(cov4) <- 1
cov5 <- matrix(lambda, nrow = 10, ncol = 10); diag(cov5) <- 1
n1 <- n1s
n2 <- g*n1
n3 <- (g^2)*n1
n4 <- g*n1
n5 <- (g^2)*n1
N <- n1 + n2 + n3 + n4 + n5

set.seed(12)

C1 <- rmst(n1, as.vector(mu1), cov1, rep(eta, 10), nu=df)
C2 <- rmst(n2, as.vector(mu2), cov2, rep(eta, 10), nu=df)
C3 <- rmst(n3, as.vector(mu3), cov3, rep(eta, 10), nu=df)
C4 <- rmst(n4, as.vector(mu4), cov4, rep(eta, 10), nu=df)
C5 <- rmst(n5, as.vector(mu5), cov5, rep(eta, 10), nu=df)
X <- rbind(C1, C2, C3, C4, C5)
data <- data.frame(X)
colnames(data) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")

treat <- as.matrix(rep(c(1,2,3,4,5), c(n1,n2,n3,n4,n5)))
treat <- factor(treat, levels = c(1, 2, 3, 4, 5), 
                labels = c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"))

data <- cbind(data, treat)

# Generalized propensity scores for three treatments
fit <- multinom(treat ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, 
                data = data, trace = F)
Rx <- fitted(fit) 
colnames(Rx) <- c("p1", "p2", "p3", "p4", "p5")
data <- cbind(data, Rx)

# Determine eligibility
min.max.Ps <- data %>%
  group_by(treat) %>%
  summarise(min1 = min(p1), max1 = max(p1), 
            min2 = min(p2), max2 = max(p2), 
            min3 = min(p3), max3 = max(p3),
            min4 = min(p4), max4 = max(p4), 
            min5 = min(p5), max5 = max(p5))

data$Eligible <- 
  data$p1 >= max(min.max.Ps$min1) & data$p1 <= min(min.max.Ps$max1) &
  data$p2 >= max(min.max.Ps$min2) & data$p2 <= min(min.max.Ps$max2) &
  data$p3 >= max(min.max.Ps$min3) & data$p3 <= min(min.max.Ps$max3) & 
  data$p4 >= max(min.max.Ps$min4) & data$p4 <= min(min.max.Ps$max4) & 
  data$p5 >= max(min.max.Ps$min5) & data$p5 <= min(min.max.Ps$max5)

data <- dplyr::filter(data, Eligible)

# Calculate new propensity scores for eligible subjects
fit.E <- multinom(treat ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, 
                  data = data, trace = F)
Rx.E <- fitted(fit.E) 
colnames(Rx.E) <- c("p1", "p2", "p3", "p4", "p5")
data <- data %>% 
  dplyr :: select(-p1, -p2, -p3, -p4, -p5)
data <- cbind(data, Rx.E)

data$T1 <- data$treat == "Treatment 1"
data$T2 <- data$treat == "Treatment 2"
data$T3 <- data$treat == "Treatment 3"
data$T4 <- data$treat == "Treatment 4"
data$T5 <- data$treat == "Treatment 5"

n1 <- sum(data$treat == "Treatment 1")
n2 <- sum(data$treat == "Treatment 2")
n3 <- sum(data$treat == "Treatment 3")
n4 <- sum(data$treat == "Treatment 4")
n5 <- sum(data$treat == "Treatment 5")

log <- data.frame(cbind(logit(data$p1), logit(data$p2), logit(data$p3), logit(data$p4), logit(data$p5)))

# Matching performance
# Pre-matched cohort
covariates <- data %>%
  dplyr::select(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)

sds <- c()
for (i in 1:ncol(covariates)){
  sds[i] <- sd(covariates[data$T1==T, i], na.rm = T)
}
msbs.pre <- c()
for (i in 1:ncol(covariates)){
  msbs.pre[i] <- match5.pre(sds[i], 
                            covariates[,i][which(data$T1==T)], 
                            covariates[,i][which(data$T2==T)], 
                            covariates[,i][which(data$T3==T)], 
                            covariates[,i][which(data$T4==T)], 
                            covariates[,i][which(data$T5==T)])
}
max2sb.max.pre <- max(msbs.pre, na.rm = T)

# Vector matching (Lopez and Gutman, 2016)
clustnum <- 5

temp345 <- kmeans(cbind(logit(data$p3), logit(data$p4), logit(data$p5)), clustnum)
data$Quint345 <- temp345$cluster

temp245 <- kmeans(cbind(logit(data$p2), logit(data$p4), logit(data$p5)), clustnum)
data$Quint245 <- temp245$cluster

temp235 <- kmeans(cbind(logit(data$p2), logit(data$p3), logit(data$p5)), clustnum)
data$Quint235 <- temp235$cluster

temp234 <- kmeans(cbind(logit(data$p2), logit(data$p3), logit(data$p4)), clustnum)
data$Quint234 <- temp234$cluster

temp12 <- dplyr::filter(data, treat != "Treatment 3" & treat != "Treatment 4" & treat != "Treatment 5")
temp13 <- dplyr::filter(data, treat != "Treatment 2" & treat != "Treatment 4" & treat != "Treatment 5")
temp14 <- dplyr::filter(data, treat != "Treatment 2" & treat != "Treatment 3" & treat != "Treatment 5")
temp15 <- dplyr::filter(data, treat != "Treatment 2" & treat != "Treatment 3" & treat != "Treatment 4")

# t1 = reference treatment
match12 <- Matchby(Y = NULL, Tr = temp12$treat == "Treatment 1", 
                   X = logit(temp12$p1), by = temp12$Quint345, 
                   caliper = 0.5, replace = T, estimand = "ATT")
match13 <- Matchby(Y = NULL, Tr = temp13$treat == "Treatment 1", 
                   X = logit(temp13$p1), by = temp13$Quint245, 
                   caliper = 0.5, replace = T, estimand = "ATT")
match14 <- Matchby(Y = NULL, Tr = temp14$treat == "Treatment 1", 
                   X = logit(temp14$p1), by = temp14$Quint235, 
                   caliper = 0.5, replace = T, estimand = "ATT")
match15 <- Matchby(Y = NULL, Tr = temp15$treat == "Treatment 1", 
                   X = logit(temp15$p1), by = temp15$Quint234, 
                   caliper = 0.5, replace = T, estimand = "ATT")

vm.results <- match5.results()
pm.vm <- vm.results[1]
vm.maxmax <- vm.results[2]

# FM Matching
fuzz345 <- cmeans(data.frame(logit(data$p3), logit(data$p4), logit(data$p5)), 4)
membership345 <- data.frame(fuzz345$membership)

cutoff <- 0.25
clust345 <- list()
for (i in 1:nrow(membership345)){
  clust345[[i]] <- which(membership345[i,] >= cutoff)
}

fuzz245 <- cmeans(data.frame(logit(data$p2), logit(data$p4), logit(data$p5)), 4)
membership245 <- data.frame(fuzz245$membership)

clust245 <- list()
for (i in 1:nrow(membership245)){
  clust245[[i]] <- which(membership245[i,] >= cutoff)
}

fuzz235 <- cmeans(data.frame(logit(data$p2), logit(data$p3), logit(data$p5)), 4)
membership235 <- data.frame(fuzz235$membership)

clust235 <- list()
for (i in 1:nrow(membership235)){
  clust235[[i]] <- which(membership235[i,] >= cutoff)
}

fuzz234 <- cmeans(data.frame(logit(data$p2), logit(data$p3), logit(data$p4)), 4)
membership234 <- data.frame(fuzz234$membership)

clust234 <- list()
for (i in 1:nrow(membership234)){
  clust234[[i]] <- which(membership234[i,] >= cutoff)
}

match12 <- Matchby.fuzzy(Y = NULL, Tr = temp12$treat == "Treatment 1", 
                         X = cbind(logit(temp12$p1), logit(temp12$p2)), 
                         by = clust345[c(1:(n1+n2))], 
                         caliper = 0.5, replace = T, estimand = "ATT")
match13 <- Matchby.fuzzy(Y = NULL, Tr = temp13$treat == "Treatment 1", 
                         X = cbind(logit(temp13$p1), logit(temp13$p3)), 
                         by = clust245[c(1:n1, (n1+n2+1):(n1+n2+n3))], 
                         caliper = 0.5, replace = T, estimand = "ATT")
match14 <- Matchby.fuzzy(Y = NULL, Tr = temp14$treat == "Treatment 1", 
                         X = cbind(logit(temp14$p1), logit(temp14$p4)), 
                         by = clust235[c(1:n1, (n1+n2+n3+1):(n1+n2+n3+n4))], 
                         caliper = 0.5, replace = T, estimand = "ATT")
match15 <- Matchby.fuzzy(Y = NULL, Tr = temp15$treat == "Treatment 1", 
                         X = cbind(logit(temp15$p1), logit(temp15$p5)), 
                         by = clust234[c(1:n1, (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5))], 
                         caliper = 0.5, replace = T, estimand = "ATT")

fm.results <- match5.results()
pm.fm <- fm.results[1]
fm.maxmax <- fm.results[2]

# FMnc Matching
match12 <- Matchby.fuzzy(Y = NULL, Tr = temp12$treat == "Treatment 1", 
                         X = cbind(logit(temp12$p1), logit(temp12$p2)), 
                         by = clust345[c(1:(n1+n2))], 
                         caliper = NULL, replace = T, estimand = "ATT")
match13 <- Matchby.fuzzy(Y = NULL, Tr = temp13$treat == "Treatment 1", 
                         X = cbind(logit(temp13$p1), logit(temp13$p3)), 
                         by = clust245[c(1:n1, (n1+n2+1):(n1+n2+n3))], 
                         caliper = NULL, replace = T, estimand = "ATT")
match14 <- Matchby.fuzzy(Y = NULL, Tr = temp14$treat == "Treatment 1", 
                         X = cbind(logit(temp14$p1), logit(temp14$p4)), 
                         by = clust235[c(1:n1, (n1+n2+n3+1):(n1+n2+n3+n4))], 
                         caliper = NULL, replace = T, estimand = "ATT")
match15 <- Matchby.fuzzy(Y = NULL, Tr = temp15$treat == "Treatment 1", 
                         X = cbind(logit(temp15$p1), logit(temp15$p5)), 
                         by = clust234[c(1:n1, (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5))], 
                         caliper = NULL, replace = T, estimand = "ATT")

fmnc.results <- match5.results()
pm.fmnc <- fmnc.results[1]
fmnc.maxmax <- fmnc.results[2]

# GPSnc Matching
match12 <- Match(Y = NULL, Tr = temp12$treat == "Treatment 1", 
                 X = log[1:(n1+n2), ],  
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match13 <- Match(Y = NULL, Tr = temp13$treat == "Treatment 1", 
                 X = log[c(1:n1, (n1+n2+1):(n1+n2+n3)), ], 
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match14 <- Match(Y = NULL, Tr = temp14$treat == "Treatment 1", 
                 X = log[c(1:n1, (n1+n2+n3+1):(n1+n2+n3+n4)), ],  
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match15 <- Match(Y = NULL, Tr = temp15$treat == "Treatment 1", 
                 X = log[c(1:n1, (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)), ], 
                 replace = T, ties = F, version = "fast", estimand = "ATT")

gpsnc.results <- match5.results()
pm.gpsnc <- gpsnc.results[1]
gpsnc.maxmax <- gpsnc.results[2]

# COVnc Matching
match12 <- Match(Y = NULL, Tr = temp12$treat == "Treatment 1", 
                 X = covariates[1:(n1+n2), ],  
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match13 <- Match(Y = NULL, Tr = temp13$treat == "Treatment 1", 
                 X = covariates[c(1:n1, (n1+n2+1):(n1+n2+n3)), ], 
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match14 <- Match(Y = NULL, Tr = temp14$treat == "Treatment 1", 
                 X = covariates[c(1:n1, (n1+n2+n3+1):(n1+n2+n3+n4)), ], 
                 replace = T, ties = F, version = "fast", estimand = "ATT")
match15 <- Match(Y = NULL, Tr = temp15$treat == "Treatment 1", 
                 X = covariates[c(1:n1, (n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)), ], 
                 replace = T, ties = F, version = "fast", estimand = "ATT")

covnc.results <- match5.results()
pm.covnc <- covnc.results[1]
covnc.maxmax <- covnc.results[2]



