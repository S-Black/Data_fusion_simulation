##############################################################################
# This program runs two statistical matching algorithms: DPMPM and hot deck
# Inputs: Generated data sets for various scenarioes (6 corr x 3 HD)
# Steps:
#       1) Read in data
#       2) Create data sets with various sample sizes, variable numbers
#       3) Run statistical matching algorithms
#       4) Assess statistical matching with Raessler checks
#       5) Print charts with results
##############################################################################
library(NPBayesImpute)
library(StatMatch)

#Functions
#-------------------------------------------------

# Functions creates a concatenated data set with nans if there are no data
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

#Calculates the percentage of wrongly matched values
AvgIncorrect <- function(Orig,Est){
  t <- table(Orig, Est)
  error <- sum(t*toeplitz(c(0,rep(1,nrow(t)-1))))/sum(t)
  return(error)
}

#Hellinger distance between two discrete distributions.
hellinger <- function(p, q){
  if(nrow(p)>50){
    p = as.vector(as.matrix(ftable(prop.table(table(p)))))
    q = as.vector(as.matrix(ftable(prop.table(table(q)))))
  }
  z = sqrt(p) - sqrt(q)
  return(sqrt(z%*%z/2))
}

#Cramers V
cv.test <- function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  #print.noquote("Cramér V:")
  return(as.numeric(CV))
}

#Correlation matrix
Corr_Cramer <- function(df){
  df = as.data.frame(df)
  corr = as.data.frame(matrix(0, ncol = ncol(df), nrow = ncol(df)))
  names(corr) <- names(df)
  for(i in 1:ncol(corr)){
    for(j in 1:ncol(corr)){
      if(i == j)
        corr[i,j] <- 1.0
      else{
        cv <- cv.test(df[[i]],df[[j]])
        corr[i,j] <- cv
        corr[j,i] <- cv
      }
    }
  }
  return(corr)
}


##################################

Run_Statmatch <- function(z_vars){
  # Data settings
  # 6 correlation scenarioes
  corr_scen <- list('highhigh','highlow','lowhigh','lowlow','high','low')
  n_base <- 2000
  n_donor <- 2500
  n_aux <- 1000

  #DPMPM Model settings
  iterations <- 30000
  burn <- 20000 #burnins
  thin <- 10 #thin every # iteration
  Kmax <- 50 #maximum number of mixture components
  aalpha <- 0.25 #the hyper parameter ’a’ for alpha in stick-breaking prior distribution.
  balpha <- aalpha #the hyper parameter ’b’ for alpha in stick-breaking prior distribution.

  to_save <- NULL

for (c in 1:length(corr_scen)){ #6 correlation scanrioes

  # Data prep
  ########################################################################
  #Read data
  x <-  readRDS(file = paste0('HD_low_corr_',corr_scen[[c]],'_v2.rds'))

  #draw sample for all data sets because we don't want accidentally double draw an observation
  D <- x[sample(nrow(x), n_base+n_donor+n_aux, replace = FALSE, prob = NULL),]
  D1 <- D[1:n_base,] #base/recipient data
  D2 <- D[(n_base+1):(n_base+n_donor),] #donor/supplement data
  aux_data <- D[(n_base+n_donor+1):(n_base+n_donor+n_aux),(names(D) %in% c('Y1','Y2','Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8'))] #glue/auxiliary information

  #keep Z variables from base file
  drops <- c('X1','X2','X3','X4','X5','Y1','Y2')
  Z <- D1[, !(names(D1) %in% drops)]

  # delete Z variables from base data set
  D1 <- D1[,!(names(D1) %in% names(Z))]

  # delete Y variables from donor file
  D2 <- D2[,!(names(D2) %in% c('Y1','Y2'))]
  # number of matching variables
  if(z_vars==2){
    D2 <- D2[,!(names(D2) %in% c('Z3','Z4','Z5','Z6','Z7','Z8'))]
    aux_data <- aux_data[,!(names(aux_data) %in% c('Z3','Z4','Z5','Z6','Z7','Z8'))]
  }

  # Statistical Matching
  ##################################################################
  #DPMPM
  #######

  # Concatenate three datasets
  temp <- concat_nan(D2,aux_data)
  cX <- concat_nan(D1,temp)
  cX <- lapply(cX[1:ncol(cX)], factor)
  cX <- as.data.frame(cX)

  # run model and save output
  model <- CreateModel(cX,MCZ=NULL,Kmax,200000,aalpha,balpha)
  start_time <- Sys.time()
  model$Run(burn,iterations,thin)
  end_time <- Sys.time()
  end_time - start_time

  #retrieve parameters from the final iteration
  result <- model$snapshot

  #convert ImputedX matrix to dataframe, using proper factors/names etc.
  IX <- GetDataFrame(result$ImputedX,cX)

  #only keep base file data
  start <-1
  end <- nrow(D1)
  SD_DPMPM <- IX[start:end,]
  gc()

  #StatMatch
  ##########
  # Step 1) Match Z from aux to D1 using Y distance
  out.NND.1 <- NND.hotdeck(data.rec=D1, data.don=aux_data,
               match.vars=c('Y1','Y2'),dist.fun='Gower')

  # create the synthetic data.set:
  fused.1 <- create.fused(data.rec=D1, data.don=aux_data,mtc.ids=out.NND.1$mtc.ids,
             z.vars=names(aux_data[ , grepl("Z",names(aux_data))]))

  # Step 2) Replace Z by matching Z from D2 to D1 using (x,Z) distance
  out.NND.2 <- NND.hotdeck(data.rec=fused.1, data.don=D2,
               match.vars=c(names(D2[ , grepl("X",names(D2))]),names(D2[ , grepl("Z",names(D2))])),dist.fun="Gower")

  SD_HD <- create.fused(data.rec=D1, data.don=D2,
           mtc.ids=out.NND.2$mtc.ids, z.vars=names(D2[ , grepl("Z",names(D2))]))

  gc()
  #############################################################################
  # Assess statistical matching performance
  # Calculate Raessler (2014) checks for original data and all synthetic datasets
  # Compare results
  #################
  n_matchvars <- names(SD_DPMPM[ , grepl("Z",names(SD_DPMPM))])
  orig <- cbind(D1,Z)

  #Joint distribution Y and Z
  v_unique <- c(n_matchvars, names(orig[ , grepl("Y",names(orig))]))
  J_HD <- c(hellinger(orig[,v_unique],SD_DPMPM[,v_unique]),hellinger(orig[,v_unique],SD_HD[,v_unique]))

  for(j in 1:length(n_matchvars)){

    var=n_matchvars[j] #matched variable we want to assess

    # i.   Preserve individual values of Z
    #----------------------------------
    #Average number of individuals in the incorrect cell of the contingency table
    e <- c(AvgIncorrect(orig[[var]],SD_DPMPM[[var]]),AvgIncorrect(orig[[var]],SD_HD[[var]]))

    # ii.  Preserve joint distribution
    #------------------------------------

    #Hellinger distance (0 means distributions are the same)
    common_vars <- c(names(D1[ , grepl("X",names(D1))]),var)
    HD <- c(hellinger(orig[,common_vars],SD_DPMPM[,common_vars]),hellinger(orig[,common_vars],SD_HD[,common_vars]))

    # iii. Preserve correlation structure
    #------------------------------------

    # Cramers V
    Cor_D1 <- Corr_Cramer(orig[,common_vars])
    Cor1 <- Corr_Cramer(SD_DPMPM[,common_vars])
    Cor2 <- Corr_Cramer(SD_HD[,common_vars])
    #print(Cor)
    Cor = c(sum(sum((Cor1-Cor_D1)**2)),sum(sum((Cor2-Cor_D1)**2)))

    # iv.  Preserve marginal distributions
    #-------------------------------------
    MargOrig <- table(orig[var])/nrow(D1)
    MargDPMPM <- table(SD_DPMPM[var])/nrow(SD_DPMPM)
    MargHD <- table(SD_HD[var])/nrow(SD_HD)

    Marg_HD <- c(hellinger(MargOrig,MargDPMPM),hellinger(MargOrig,MargHD))

    #Save calcualted checks
    SummaryMeasures <- data.frame(e,HD,Cor,Marg_HD)
    colnames(SummaryMeasures) <- c('incorrect','joint HD','Correlation','marginal HD')
    rownames(SummaryMeasures) <- c('DPMPM','Hot Deck')

    temp<-data.frame(SummaryMeasures,var,z_vars,corr_scen[[c]],J_HD)
    colnames(temp) <- c('incorrect','joint','correlation','marginal','z vars','corr scen','yz_joint')
    to_save <- rbind(to_save,temp)
    gc()
    }
  }
  return(to_save)
}

#vars5 <- TRUE #10 or 5 variables
#z_vars <- 2 #2, 4 or 8
to_save <- Run_Statmatch(2)
saveRDS(to_save, file = paste0('StatMatch_v2_Checks_z2.rds'))
rm(to_save)
gc()
to_save <- Run_Statmatch(8)
saveRDS(to_save, file = paste0('StatMatch__v2_Checks_z8.rds'))
rm(to_save)
gc()
