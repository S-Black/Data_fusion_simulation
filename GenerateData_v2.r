################################################################################
# Generate multinomial data sets that will be used to test and
# compare statistical matching algorithms
#
# Scenarioes:
# - max. category: 5
# - variables: matching: 2, 8 all 5 categories
# - correlation: low/high between common and matching variables
##############################################################################

#Libraries
library(matrixcalc) #matrix calc
library(simstudy) #data generation library

cv.test <- function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cramér V:")
  return(as.numeric(CV))
}

#function to generate initialisation for categories
RandomInit <- function(j){
  # j: number of categories for the variable
  x <- runif(j)
  z <- round(x/sum(x),3)
  while(sum(z)!=1){ #make sure rounded proportions sum exactly to 1, else data can't be generated
    x <- runif(j)
    z <- round(x/sum(x),3)
  }
  return(z)  #return vector with random numbers that sum to 1
}

set.seed(111)
########################

# COMMON VARIABLES (X): 5
#proportion of observations per category (random)
#need to sum to one
#x1 <- RandomInit(5) #X1 (age)
x1 <- c(0.179, 0.118, 0.142, 0.16, 0.401)
#x2<-RandomInit(2) #X2 (gender) BIASED
x2<-c(0.6,0.4) #LCF
#x3<-RandomInit(5) #X3 (income) BIASED
x3 <- c(0.122, 0.247, 0.107, 0.298, 0.226)
#x4<-RandomInit(3) #X4 (region)
x4 <- c(0.45,0.3,0.25)
#x5<-RandomInit(4) #X5
x5 <- c(0.266, 0.333, 0.265, 0.136)

# UNIQUE VARIABLES (Y): 2
#y1 <- RandomInit(3)
y1 <- c(0.284, 0.419, 0.297)
#y2 <- RandomInit(4)
y2 <- c(0.284, 0.358, 0.211, 0.147)

# MATCHING VARIABLES (Z): 2
#z1 <- RandomInit(5)
z1 <- c(0.12, 0.296, 0.165, 0.133, 0.286)
#RandomInit(5)
z2 <- c(0.179, 0.118, 0.142, 0.16, 0.401)
#z3 <- RandomInit(5)
z3 <- c(0.122, 0.247, 0.107, 0.298, 0.226)
#z4 <- RandomInit(5)
z4 <- c(0.031, 0.274, 0.25, 0.252, 0.193)
#z5 <- RandomInit(5)
z5 <- c(0.368, 0.205, 0.145, 0.105, 0.177)
#z6 <- RandomInit(5)
z6 <- c(0.206, 0.106, 0.204, 0.299, 0.185)
#z7 <- RandomInit(5)
z7 <- c(0.119, 0.218, 0.234, 0.245, 0.184)
#z8 <- RandomInit(5)
z8 <- c(0.198, 0.254, 0.19, 0.223, 0.135)

#CORRELATION
GenCorData <- function(N,r,s){
  #x1,y1,z1
  baseprobs <- matrix(c(y1[1],   y1[2],  y1[3],     0,     0, #Y1
                        x1[1],   x1[2],  x1[3], x1[4], x1[5], #X1
                        z1[1],   z1[2],  z1[3], z1[4], z1[5]) #Z1
                        ,nrow = 3, byrow = TRUE)

                     #Y1 X1 Z1
  corMat <- matrix(c(1, r, s, #Y1
                     r, 1, r, #X1
                     s, r, 1),#Z1
                     nrow = 3)

  is.positive.definite(corMat)
  is.positive.semi.definite(corMat)

  dT <- genData(N)
  dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs, prefix = "q", corMatrix = corMat, corstr = "cs")
  gendat<-as.data.frame(dX)
  drops <- c("id")
  gendat<-gendat[ , !(names(gendat) %in% drops)]
  colnames(gendat) <- c('Y1','X1','Z1')

  #x3,y2,z2
  baseprobs <- matrix(c(y2[1], y2[2], y2[3], y2[4],     0, #Y2
                        x3[1], x3[2], x3[3], x3[4], x3[5], #X3
                        z2[1], z2[2], z2[3], z2[4], z2[5]) #Z2
                        ,nrow = 3, byrow = TRUE)

  dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs, prefix = "q", corMatrix = corMat, corstr = "cs")
  gendat <- cbind(gendat, as.data.frame(dX))
  drops <- c("id")
  gendat<-gendat[ , !(names(gendat) %in% drops)]
  colnames(gendat) <- c('Y1','X1','Z1','Y2','X3','Z2')


  #x2,z3,z4
  baseprobs <- matrix(c(x2[1], x2[2],     0,     0,     0, #X2
                        z3[1], z3[2], z3[3], z3[4], z3[5], #Z3
                        z4[1], z4[2], z4[3], z4[4], z4[5]) #Z4
                        ,nrow = 3, byrow = TRUE)

                    #X2 Z3 Z4
  corMat <- matrix(c(1, r, r, #X3
                     r, 1, 0, #Z3
                     r, 0, 1),#Z4
                     nrow = 3)

  is.positive.definite(corMat)
  is.positive.semi.definite(corMat)

  dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs, prefix = "q", corMatrix = corMat, corstr = "cs")
  gendat <- cbind(gendat, as.data.frame(dX))
  drops <- c("id")
  gendat<-gendat[ , !(names(gendat) %in% drops)]
  colnames(gendat) <- c('Y1','X1','Z1','Y2','X3','Z2','X2','Z3','Z4')

  #x4,z5,z6
  baseprobs <- matrix(c(x4[1], x4[2], x4[3],     0,     0, #X4
                        z5[1], z5[2], z5[3], z5[4], z5[5], #Z5
                        z6[1], z6[2], z6[3], z6[4], z6[5]) #Z6
                        ,nrow = 3, byrow = TRUE)

  dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs, prefix = "q", corMatrix = corMat, corstr = "cs")
  gendat <- cbind(gendat, as.data.frame(dX))
  drops <- c("id")
  gendat<-gendat[ , !(names(gendat) %in% drops)]
  colnames(gendat) <- c('Y1','X1','Z1','Y2','X3','Z2','X2','Z3','Z4','X4','Z5','Z6')

  #x5,z7,z8
  baseprobs <- matrix(c(x5[1], x5[2], x5[3], x5[4],     0, #X5
                        z7[1], z7[2], z7[3], z7[4], z7[5], #Z7
                        z8[1], z8[2], z8[3], z8[4], z8[5]) #Z8
                        ,nrow = 3, byrow = TRUE)

  dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs, prefix = "q", corMatrix = corMat, corstr = "cs")
  gendat <- cbind(gendat, as.data.frame(dX))
  drops <- c("id")
  gendat<-gendat[ , !(names(gendat) %in% drops)]
  colnames(gendat) <- c('Y1','X1','Z1','Y2','X3','Z2','X2','Z3','Z4','X4','Z5','Z6','X5','Z7','Z8')

  return(gendat)
}

#########
# Save data to file

N <- 100000 #number of data samples to generate
#r <- .7 #correlation between common variables
#s <- .7 #correlation between matching variables (Y,Z)
#GenCorData(N,r,s)
high <- 0.7
low <- 0.1
t1 <- GenCorData(N,high,high)
saveRDS(t1, file = 'HD_low_corr_highhigh_v2.rds') #high correlation between X,Z and X,Y - high correlations between Y and Z
t2 <- GenCorData(N,high,low)
saveRDS(t2, file = 'HD_low_corr_highlow_v2.rds') #high correlation between X,Z and X,Y - low correlations between Y and Z
t3 <- GenCorData(N,low,high)
saveRDS(t3, file = 'HD_low_corr_lowhigh_v2.rds') #low correlation between X,Z and X,Y - high correlations between Y and Z
t4 <- GenCorData(N,low,low)
saveRDS(t4, file = 'HD_low_corr_lowlow_v2.rds') #low correlation between X,Z and X,Y - low correlations between Y and Z
t5 <- GenCorData(N,high,0)
saveRDS(t5, file = 'HD_low_corr_high_v2.rds') #high correlation between X,Z and X,Y - no correlations between Y and Z
t6 <- GenCorData(N,low,0)
saveRDS(t6, file = 'HD_low_corr_low_v2.rds') #low correlation between X,Z and X,Y - no correlations between Y and Z
