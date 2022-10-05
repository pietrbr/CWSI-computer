rm(list=ls())
source('CWSI_fun_noTd.R')
library(sensitivity)
library(parallel)
library(pbapply)

inputA <- read.table('data_calosso/A.txt', sep = ';', header=T)
inputB <- read.table('data_calosso/B.txt', sep = ';', header=T)
inputC <- read.table('data_calosso/C.txt', sep = ';', header=T)
inputD <- read.table('data_calosso/D.txt', sep = ';', header=T)
nA <- nrow(inputA)
nB <- nrow(inputB)
nC <- nrow(inputC)
nD <- nrow(inputD)

#####################################
# outlier correction A
inputA$CANOPY_TEMP[2] <- 22
inputA$CANOPY_TEMP[31] <- 26
inputA$CANOPY_TEMP[37] <- 27.5
inputA$CANOPY_TEMP[10] <- 23
#inputA$WIND_SPEED <- rep(mean(inputA$WIND_SPEED), nA)


#####################################
# outlier correction B
inputB$CANOPY_TEMP[42] <- mean(inputB$CANOPY_TEMP, na.rm=T)
inputB$WIND_SPEED <- rep(mean(inputB$WIND_SPEED, na.rm=T), nB)
inputB$UV_RAD[35] <- mean(inputB$UV_RAD, na.rm=T)


######################################
# outlier correction C
inputC$CANOPY_TEMP[1] <- 22
inputC$CANOPY_TEMP[5] <- 22.5
inputC$CANOPY_TEMP[13] <- 23
inputC$CANOPY_TEMP[23] <- 24
#inputC$WIND_SPEED <- rep(mean(inputC$WIND_SPEED, na.rm=T), nC)


# unite the tables
input <- rbind(inputA, inputB, inputC, inputD)
n <- nrow(input)



#################################################
# Ta temperature setting with Calosso data

Ta_minA <- 22 #A
Ta_maxA <- 25 #A
Ta_correttoA <- seq(Ta_minA, Ta_maxA, length=nA)

Ta_minB <- 27.3  #B
Ta_maxB <- 28.5  #B
Ta_correttoB <- seq(Ta_minB, Ta_maxB, length=nB)

Ta_minC <- 25.0  #C
Ta_maxC <- 25.8  #C
Ta_correttoC <- seq(Ta_minC, Ta_maxC, length=nC)

Ta_minD <- 28  #D
Ta_maxD <- 28  #D
Ta_correttoD <- seq(Ta_minD, Ta_maxD, length=nD)

Ta_corretto <- c(Ta_correttoA, Ta_correttoB, Ta_correttoC, Ta_correttoD)


# build the dataframe
df <- data.frame(
  Ta = Ta_corretto,
  Tc = input$CANOPY_TEMP,
  UV = input$UV_RAD,
  p = input$PRESSURE,
  Ur = input$HUM,
  u2 = input$WIND_SPEED,
  h = rep(1.5, n)
)

# remove NA and outlier
df <- df[-c(67,126),]
n <- nrow(df)

# function to compute the CWSI
CWSI_final_fun <- function(X){
  y <- CWSI_fun_noTd(X[,1], X[,2], X[,3], X[,4], X[,5], X[,6], X[,7])
}

# CWSI 1.5m for all the measurments
ris <- CWSI_final_fun(df)



##########################################################################################
# Compute Sobol Indices (Single Split)

# split the data in two equal size samples (needed to compute Sobol Indices)
id1 <- sample(1:n, n/2, replace = F)
id2 <- setdiff(1:n, id1)

# Different Algorithms to compute Sobol Indices, Martinez provides the best results
# x <- sobolEff(model=CWSI_final_fun, df[id1,], df[id2,], nboot=10000, order=1)
# x <- soboljansen(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, conf = 0.95)
x <- sobolmartinez(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, conf=0.95)

plot(x)
title(main='Sensitivity Analysis, Sobol Indices (Single Split)')
grid()

print(x)



###########################################################################################
# Compute Sobol Indices (Mean of Multiple Splits)


# function that returns matrices of First Order Indices and Total Effect Indices
# S: first order; T: total
sobol_compute <- function(b){
  id1 <- sample(1:n, n/2, replace = F)
  id2 <- setdiff(1:n, id1)
  
  # Sobol Indices with Martinez Algo
  x_b <- sobolmartinez(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, conf=0.95)
  final_S <- x_b$S
  final_T <- x_b$T
  
  # Sobol Indices with Asymptotic Efficient Formulas
  # x_b_S <- sobolEff(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, order=1)
  # x_b_T <- sobolEff(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, order=0)
  # final_S <- x_b_S$S
  # final_T <- x_b_T$S
  
  return(list('S'=final_S, 'T'=final_T))
}


# functions to parallelize the process
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list( "n", "CWSI_final_fun", 'CWSI_fun_noTd', 'df', 'sobolmartinez', 'sobol_compute', 'soboljansen', 'sobolEff'))
wrapper <- function(b){sobol_compute(b)} 


# compute Sobol Indices for B different splits, collect them in S_T_list
set.seed(1234)
B <- 1000
b_list <- 1:B
S_T_list <- pbsapply(b_list, wrapper, cl=cl)

S_list <- S_T_list[1,]
T_list <- S_T_list[2,]


# Do the mean of the computed Sobol Indices
S_final <- Reduce('+', S_list)/B
T_final <- Reduce('+', T_list)/B
rownames(S_final) <- rownames(x$S)
rownames(T_final) <- rownames(x$T)


# Visualize the result
x$S <- S_final
x$T <- T_final

plot(x)
title(main='Sensitivity Analysis, Sobol Indices')
grid()


print('First Order Effect')
print(S_final)

print('Total Effect')
print(T_final)




#########################################################################################
# Second Order Effect


# Find Correct names of the pairs of variables
list_names <- rownames(x$S)
nm <- list()
for(i in 1:6){
  for(j in (i+1):7){
    nm <- append(nm, paste0(list_names[i],'-',list_names[j]))
  }
}

# Second Order Effects
x2 <- sobolEff(model=CWSI_final_fun, df[id1,], df[id2,], nboot=10000, order=2)
rownames(x2$S) <- nm

par(las=2)
plot(x2)
title(main='Second Order Effects (Single Split)')
#grid()
print(x2$S)



###########################################################################################
# Compute Second Order Effect (Mean of Multiple Splits)


# function that returns matrices of Second Order Effects
sobol_2_compute <- function(b){
  id1 <- sample(1:n, n/2, replace = F)
  id2 <- setdiff(1:n, id1)
  
  
  # Sobol Indices with Asymptotic Efficient Formulas
  x_b_2 <- sobolEff(model=CWSI_final_fun, df[id1,], df[id2,], nboot=1000, order=2)
  final_S <- x_b_2$S
  
  return(list(final_S))
}


# functions to parallelize the process
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = list( "n", "CWSI_final_fun", 'CWSI_fun_noTd', 'df', 'sobolmartinez', 'sobol_2_compute', 'soboljansen', 'sobolEff'))
wrapper2 <- function(b){sobol_2_compute(b)} 


# compute Sobol Indices for B different splits, collect them in S_T_list
set.seed(1234)
B <- 1000
b_list <- 1:B
S2_list <- pbsapply(b_list, wrapper2, cl=cl)


# Do the mean of the computed Sobol Indices
S2_final <- Reduce('+', S2_list)/B


# Visualize the result
rownames(S2_final) <- nm
x2$S <- S2_final

par(las=2)
#grid(lty='solid')
plot(x2)
title(main='Sensitivity Analysis, Second Order Effects')



print('Second Order Effect')
print(S2_final)
