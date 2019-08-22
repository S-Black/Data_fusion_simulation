####################################################################################
# This program will simulate statistical matching with generated multinomial data
#
# The following procedures will be employed:
# 1) StatMatch (Hot Deck) - StatMatch in R
# 2) Regression - Regression in Stan
# 3) DPMPM - NPBayesImpute R package
# 4) Plain imputation with MICE R package
#
###################################################################################

#Libraries
library(feather)
library(StatMatch)
#library(rstan)
library(rstanarm)
#library(shinystan)
library(NPBayesImpute)
library(mice)

#Functions
# Function creates a concatenated data set with nans if there are no data
concat_nan <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)

  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)

  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)

  # add blank columns to d2 with missing values
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA #as.double(NA) if mode is not categorical
    }
  }

  # add blank columns to d1 with missing values
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }
  return(rbind(d1, d2))
}

# Paths
current_path <- getwd()
Data_path <- paste0(current_path,'/data/')
Stan_path <- paste0(current_path,'/Stan Code/')
#Read generated dataset
OrigData <- read_feather(paste(Data_path,'Multdata.feather',sep = ""))

#Make S and C data sets from generated data
CData <- SData <- OrigData

drops_s <- c("S1")
CData<-CData[ , !(names(CData) %in% drops_s)]
drops_c <- c("C1","C2")
SData<-SData[ , !(names(SData) %in% drops_c)]

#--------------------------------------------------------------------------------
# Regression - multinomial (in STAN)
# Calculate regression coefficients for each variable
# C1 ~ bernoulli(X)
# C2 ~ ordered logistic(X)
# S1 ~ ordered logistic(X)
# Apply regression coefficients to X data to derive synthetich dataset
#------------

SEED=222
CHAINS=4
ITER=2000

#Model C1 (binary variable)
CData <- CData
CData$C1 <- CData$C1-1
C1Fit <- stan_glm(C1 ~ X1 + X2 + X3 + X4, data = CData,
                 family = binomial(link = "logit"),
                 chains = CHAINS, seed = SEED, iter = ITER)
PPD_c1<- posterior_predict(C1Fit, newdata = as.data.frame([c('X1','X2','X3','X4')]), draws = NULL,seed = SEED)
C1_hat <- lapply(PPD_c1[1:ncol(PPD_c1)], median)

#Model C2 (ordered categorical variable)
CData <- CData
CData$C2 <- factor(CData$C2)
C2Fit <- stan_polr(C2 ~ X1 + X2 + X3 + X4, data = CData,prior = R2(0.25),chains = CHAINS, seed = SEED, iter = ITER)
PPD_c2<- posterior_predict(C2Fit, newdata = as.data.frame(SData[c('X1','X2','X3','X4')]), draws = NULL,seed = SEED)
C2_hat <- lapply(PPD_c2[1:ncol(PPD_c2)], median)

#Model S1 (ordered categorical variable)
SData_t <- SData
SData_t$S1 <- factor(SData_t$S1)
S1Fit <- stan_polr(S1 ~ X1 + X2 + X3 + X4, data = SData_t,prior = R2(0.25),chains = CHAINS, seed = SEED, iter = ITER)
PPD_s1<- posterior_predict(S1Fit, newdata = as.data.frame(CData[c('X1','X2','X3','X4')]), draws = NULL,seed = SEED)
S1_hat <- lapply(PPD_s1[1:ncol(PPD_s1)], median)

# Synthetic dataset
#Match C variables to S data
SD_RegC = data.frame(t(as.data.frame(C1_hat)),t(as.data.frame(C2_hat)),SData)
colnames(SD_RegC) <- c('C1','C2','X1','X2','X3','X4','S1')
SD_RegC$C1 <- SD_RegC$C1 + 1
write_feather(SD_RegC, paste0(Data_path,'SD_RegC.feather'))

#Match S variables to C data
SD_RegS = data.frame(CData,t(as.data.frame(S1_hat)))
colnames(SD_RegS) <- c('C1','C2','X1','X2','X3','X4','S1')
write_feather(SD_RegS, paste0(Data_path,'SD_RegS.feather'))


#--------------------------------------------------------------------------------
# StatMatch
# Hot Deck recipient/donor: NND.hotdeck
#------------
#make variables factor variables
CData_t <- lapply(CData[1:ncol(CData)], factor)
CData_t <- as.data.frame(CData_t)

SData_t <- lapply(SData[1:ncol(SData)], factor)
SData_t <- as.data.frame(SData_t)

#Match S variable to C data
out.NND <- NND.hotdeck(data.rec=CData_t, data.don=SData_t,
match.vars=c("X1","X3","X4"), don.class="X2",# don.class has no missing values and matching variables are separated for each class
dist.fun='Gower')

#Match C variables to S data
out.NNDC <- NND.hotdeck(data.rec=SData_t, data.don=CData_t,
match.vars=c("X1","X3","X4"), don.class="X2",# don.class has no missing values and matching variables are separated for each class
dist.fun='Gower')

# create synthetic data.set
#Match S variables to C data
SD_NNDS <- create.fused(data.rec=CData_t, data.don=SData_t,
mtc.ids=out.NND$mtc.ids, z.vars=c('S1'))
write_feather(SD_NNDS, paste(Data_path,'SD_NNDS.feather',sep = ""))

#Match C variables to S data
SD_NNDC <- create.fused(data.rec=SData_t, data.don=CData_t,
mtc.ids=out.NNDC$mtc.ids, z.vars=c('C1','C2'))
write_feather(SD_NNDC, paste(Data_path,'SD_NNDC.feather',sep = ""))

#--------------------------------------------------------------------------------
# DPMPM - NPBayesImpute
# Concatenate files and impute missing values
#------------

X=concat_nan(CData,SData)
X <- lapply(X[1:ncol(X)], factor)
X <- as.data.frame(X)
DPMPM <- CreateModel(X,MCZ=NULL,50,200000,0.25,0.25)
#Run model for
DPMPM$Run(ITER/2,ITER,2)
#retrieve parameters from the final iteration
result <- DPMPM$snapshot
# Synthetic dataset
SD_DPMPM <- GetDataFrame(result$ImputedX,X)

#Match S variables to C data
SD_DPMPMS <- SD_DPMPM[1:nrow(CData),]
write_feather(SD_DPMPMS, paste(Data_path,'SD_DPMPMS.feather',sep = ""))

#Match C variables to S data
start <-nrow(CData)+1
end <- nrow(CData)+nrow(SData)
SD_DPMPMC <- SD_DPMPM[start:end,]
write_feather(SD_DPMPMC, paste(Data_path,'SD_DPMPMC.feather',sep = ""))


#--------------------------------------------------------------------------------
# MICE - impute missing values
# Concatenate files and impute missing values
# Methods:
# logreg - binary data
# polyreg - unordered categorical variables
# polyr - ordered categorical variables
#------------

                                              #C1      #C2  X1 X2 X3 X4   S1
tempData <- mice(data=X,m=1,maxit=100,meth=c('logreg','polr','','','','','polr'),seed=500)
SD_mice <- complete(tempData,1)
#write_feather(SD_mice, paste(Data_path,'SD_mice.feather',sep = ""))

#Match S variables to C data
SD_miceS <- SD_mice[1:nrow(CData),]
write_feather(SD_miceS, paste(Data_path,'SD_miceS.feather',sep = ""))

#Match C variables to S data
start <-nrow(CData)+1
end <- nrow(CData)+nrow(SData)
SD_miceC <- SD_mice[start:end,]
write_feather(SD_miceC, paste(Data_path,'SD_miceC.feather',sep = ""))
