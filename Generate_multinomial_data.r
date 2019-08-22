################################################################################
# Generate multinomial data sets that will be used to test and
# compare data fusion algorithms. The generated data set will be
# split in two to test the fusion algorithms.
#
# Data set (C & X and S & X are correlated):
# + four common variables X (e.g. age, sex, ethnicity and income)
# + two C variables (not in the S dataset)
# + one S variable  (not in the C dataset)
##############################################################################

#Libraries
library(feather) #save data files that can be read by Python
library(matrixcalc) #matrix calc
library(simstudy) #data generation library

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


N=8000 #number of generated data points
P=7 #number of variables

#proportion of observations per category (random)
#need to sum to one
set.seed(111)
c1<-RandomInit(2) #C1
#sum(c1)
c2<-RandomInit(4) #C2
#sum(c2)
x1<-RandomInit(4) #X1
#sum(x1)
x2<-RandomInit(2) #X2
#sum(x2)
x3<-RandomInit(6) #X3
#sum(x3)
x4<-RandomInit(4) #X4
#sum(x4)
s1<-RandomInit(3) #S1
#sum(s1)

#Categorical variables proportions
baseprobs <- matrix(c(c1[1], c1[2],    0,      0,     0,     0, #C1
                      c2[1], c2[2], c2[3], c2[4],     0,     0, #C2
                      x1[1], x1[2], x1[3], x1[4],     0,     0, #X1
                      x2[1], x2[2],     0,     0,     0,     0, #X2
                      x3[1], x3[2], x3[3], x3[4], x3[5], x3[6], #X3
                      x4[1], x4[2], x4[3], x4[4],     0,     0, #X4
                      s1[1], s1[2], s1[3],     0,     0,     0) #S1
                      ,nrow = P, byrow = TRUE)

# Correlation matrix (C and S are correlated with X but uncorrelated with each
# other (conditionally independent); X variables are .1 correlated with each other)
r<-.1 #correlation between unrelated variables
q<-.8 #correlation between related variables

                  #C1 C2 X1 X2 X3 X4 S1
corMat <- matrix(c(1, r, r, q, r, r, 0, #C1
                   r, 1, r, r, q, r, 0, #C2
                   r, r, 1, r, r, r, q, #X1
                   q, r, r, 1, r, r, r, #X2
                   r, q, r, r, 1, r, r, #X3
                   r, r, r, r, r, 1, r, #X4
                   0, 0, q, r, r, r, 1) #S1
                  , nrow = P)

#Correlation matrix needs to be positive definite
is.positive.definite(corMat)
is.positive.semi.definite(corMat)

set.seed(333)
dT <- genData(N)
dX <- genCorOrdCat(dT, adjVar = NULL, baseprobs = baseprobs,
prefix = "q", corMatrix = corMat, corstr = "cs")

gendat<-as.data.frame(dX)
drops <- c("id")
gendat<-gendat[ , !(names(gendat) %in% drops)]

colnames(gendat) <- c('C1','C2','X1','X2','X3','X4','S1')

# Paths
current_path<-getwd()
Data_path<-paste0(current_path,'/data/')
write_feather(gendat, paste0(Data_path,'Multdata.feather'))
