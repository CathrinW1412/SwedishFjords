<Applying PbModel.R to data>
    Copyright (C) <2026>  <Cathrin Wittig>
#=============================================================================
# load data
#=============================================================================

setwd("D:/Sweden_2026/R_Scripts/Biomixing_results_S3")
load("SwedishFjords_Data.Rdata")

source("PbModel.R")

#load packages

library(FME)
library(viridis)
library(marelac)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)

#=============================================================================
# aux functions
#=============================================================================
<aux functions>
    Copyright (C) <2026>  <Sebastiaan van de Velde>
#-----------------------------------------------------------------------------
# set params
#-----------------------------------------------------------------------------

set.params <- function(df,PL,fjord,site){
  
  # Set-up for porosity
  
  selection <- (df$fjord==fjord) & (df$site == site) & (df$var=="porosity") & (df$replicate == 1)
  
  # setup porosity
  PL$por.grid <- setup.prop.1D(xy=cbind(df$depth[selection],df$value[selection]),grid=PL$grid,interpolate="linear")
  PL$svf.grid <- setup.prop.1D(xy=cbind(df$depth[selection],1-as.numeric(df$value[selection])),grid=PL$grid,interpolate="linear")
  
  PL$u.inf   <- PL$SedFlux/PL$rho.sed/(1-PL$por.grid$int[PL$N+1])  # sedimentation velocity pore water [cm yr-1]
  PL$v.inf   <- PL$SedFlux/PL$rho.sed/(1-PL$por.grid$int[PL$N+1])  # sedimentation velocity solids [cm yr-1]
  PL$v.grid  <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.grid$int[1],por.inf=PL$por.grid$int[PL$N+1],por.grid=PL$por.grid)$v # advection velocity solid
  PL$u.grid  <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.grid$int[1],por.inf=PL$por.grid$int[PL$N+1],por.grid=PL$por.grid)$u # advection velocity pore water
  PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  
  return(PL)
}

#-----------------------------------------------------------------------------
# plot function
#-----------------------------------------------------------------------------

plot.function <- function(df,fjord,site,PL=PL,Pb210.res,Chla.res){
  
  
  conv    <- PL$por.grid$mid + (1.-PL$por.grid$mid)*PL$rho.sed
  
  par(mfrow=c(1,2))
  
  # Pb210 
  
  sel <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
  
  plot(y=df$depth[sel],x=df$value[sel],
       ylim=c(20.,0.),xlim=c(0,max(df$value[sel],(Pb210.res+PL$Pb210.supp),na.rm=T)),pch=16,col="black",
       xlab="Pb210 (Bq kg-1 dwt)",ylab="depth (cm)",axes=T,cex=1.1,
       main=paste(fjord,site))
  sel.sd <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210_sd") & (df$replicate == 1)
  for (i in 1:length(df$depth[sel.sd])){
    lines(y=c(df$depth[sel.sd][i],df$depth[sel.sd][i]),
          x=c((df$value[sel][i]+df$value[sel.sd][i]),(df$value[sel][i]-df$value[sel.sd][i])),col="black",lty=1,lwd=2)
  }
  lines(x=Pb210.res+PL$Pb210.supp,y=PL$grid$x.mid,lwd=2,col="red",lty=2)
  
  # ChlA
  
  sel <- (df$fjord==fjord) & (df$site == site) & (df$var == "ChlA")
  
  plot(y=df$depth[sel],x=df$value[sel],
       ylim=c(20.,0.),xlim=c(0,max(df$value[sel],Chla.res/conv,na.rm=T)),pch=16,col="black",
       xlab="Chla (µg g-1 wwt)",ylab="depth (cm)",axes=T,cex=1.2,
       main=paste(fjord,site))
  lines(x=Chla.res/conv + PL$Chla.bg,y=PL$grid$x.mid,lwd=2,col="red",lty=2)
  
}


#---------------------------------------------------
# Best fit with CRAN:modcost (Soetaert & Petzoldt, 2010) 
#---------------------------------------------------

#---------------------------------------------------
# Modcost_Hakefjord_Head
#---------------------------------------------------

fjord <- "Hake"
site  <- "Head"  

# setup boundary conditions

PL$Pb210.supp <- 35     # Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.115  # MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 200./PL$conv.Pb210F # flux excess  Pb [Bq m-2 yr-1]

PL$F.Chla  <- 2 
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

#Prepare Lists and Matrix for results 

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # vary Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_HH <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_HH

ModCost_Result_HH <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_HH




#---------------------------------------------------
# Modcost_Hakefjord_Center
#---------------------------------------------------

fjord <- "Hake"
site  <- "Center"  

# setup boundary conditions

PL$Pb210.supp <- 48     # Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.164  # MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 350./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$F.Chla  <- 2.3 
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <-  0+0.1*j # vary x.L
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_HC <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_HC

ModCost_Result_HC <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_HC





#---------------------------------------------------
# Modcost_Hakefjord_Mouth 
#---------------------------------------------------

fjord <- "Hake"
site  <- "Mouth"  

# setup boundary conditions

PL$Pb210.supp <- 58     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.470  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 600./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$F.Chla  <- 3
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_HM <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_HM

ModCost_Result_HM <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_HM





#---------------------------------------------------
# Modcost_Gullmarsfjord_Head
#---------------------------------------------------

fjord <- "Gullmar"
site  <- "Head"  

# setup boundary conditions

PL$Pb210.supp <- 35     #Supported Pb [Bq kg-1] 
PL$SedFlux    <- 0.067 #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 430./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$F.Chla  <- 7
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

write.csv(R2_matrix, "R2.csv")
write.csv(result_matrix, "result.csv")


min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_GH <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_GH

ModCost_Result_GH <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_GH


#---------------------------------------------------
# Modcost_Gullmarsfjord_Center
#---------------------------------------------------

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    fjord <- "Gullmar"
    site  <- "Center"  
    
    # setup boundary conditions
    
    PL$Pb210.supp <- 75     #Supported Pb [Bq kg-1] 
    PL$SedFlux    <- 0.150  #MAR [g cm-2 yr-1] from Watts et al., 2024
    
    PL$F.Pb210    <- 300./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]
    
    PL$F.Chla  <- 0.2
    PL$Chla.bg <- 0.3
    
    PL$k.1 <- 0.005*365.25
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    #Due to a subsurface increase only the profile until 8 cm is considered
    sel6 <- (Data$depth < 8) 
    Data <- data.frame (name = c(rep("Pb210", length(Data$depth[sel6]))), depth = Data$depth[sel6], Pb210 = Data$Pb210[sel6])
    
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

write.csv(R2_matrix, "R2.csv")
write.csv(result_matrix, "result.csv")


min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_GC <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_GC

ModCost_Result_GC <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_GC

#---------------------------------------------------
# Modcost_Gullmarsfjord_Mouth
#---------------------------------------------------

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    fjord <- "Gullmar"
    site  <- "Mouth"  
    
    # setup boundary conditions
    
    PL$Pb210.supp <- 65     #Supported Pb [Bq kg-1]
    PL$SedFlux    <- 0.222  #MAR [g cm-2 yr-1] from Watts et al., 2024
    
    PL$F.Pb210    <- 450./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]
    
    PL$F.Chla  <- 2.5
    PL$Chla.bg <- 0.3
    
    PL$k.1 <- 0.005*365.25
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_GM <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_GM

ModCost_Result_GM <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_GM

#---------------------------------------------------
# Modcost_Byfjord_Head
#---------------------------------------------------

#Prepare Lists and Matrix for results

rows = 200
cols = 200
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_Db.0_list <- list()
Par_x.L_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL$Db.0 <- 0+0.1*i   # Db.0
    PL$x.L <- 0+0.1*j # vary x.L
    
    fjord <- "By"
    site  <- "Head"  
    
    # setup boundary conditions
    
    PL$Pb210.supp <- 80     #Supported Pb [Bq kg-1]
    PL$SedFlux    <- 0.047  #MAR [g cm-2 yr-1] from Watts et al., 2024
    
    PL$F.Pb210    <- 100./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]
    
    PL$F.Chla  <- 3.5
    PL$Chla.bg <- 0.3
    
    PL$k.1 <- 0.005*365.25
    
    #Run Model with current Db.0 and x.L
    
    PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
    output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
    Pb210.res <- output$Pb.decay
    
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site == site)  & (df$var == "Pb210") & (df$replicate == 1)
    Data <- data.frame (name = c(rep("Pb210", length(df$depth[sel.dat]))), depth = df$depth[sel.dat], Pb210 = df$value[sel.dat])
    
    out <- data.frame(depth = PL$grid$x.mid,Pb210 = Pb210.res+PL$Pb210.supp)
    
    colnames(out)=c("depth","Pb210")
    colnames(Data)=c("name","depth","Pb210")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "Pb210")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_x.L_list[[j]] <- paste(PL$x.L)
  }
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

#Print Results Matrix
rownames(result_matrix) <- Par_Db.0_list
colnames(result_matrix) <- Par_x.L_list
result_matrix

rownames(R2_matrix) <- Par_Db.0_list
colnames(R2_matrix) <- Par_x.L_list
R2_matrix

min_value <- min(result_matrix)
max_value <- max(R2_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

max_indices <- which(R2_matrix == max_value, arr.ind = TRUE)
row <- max_indices[1,1]
col <- max_indices[1,2]
Db.0 <- Par_Db.0_list[[row]]
x.L <- Par_x.L_list[[col]]

# Output

ModCost_Result_BH <- as.numeric(list(Db.0,x.L,min_value))
ModCost_Result_BH

ModCost_Result_BH <- as.numeric(list(Db.0,x.L,max_value))
ModCost_Result_BH

#---------------------------------------------------
# Modcost_Byfjord_Center
#---------------------------------------------------
#No Bioturbation observed.

#---------------------------------------------------
# Modcost_Byfjord_Mouth
#---------------------------------------------------
#No Bioturbation observed.


#---------------------------------------------------
# Extract All ModCost Results
#---------------------------------------------------

AllModCost_results_Db.0 <- data.frame (ModCost_Result_HH,ModCost_Result_HC,ModCost_Result_HM,ModCost_Result_GH,ModCost_Result_GC,ModCost_Result_GM,ModCost_Result_BH)
rownames(AllModCost_results_Db.0) <- c("Db.0","x.L","R2")
colnames(AllModCost_results_Db.0) <- c("HH","HC","HM","GH","GC","GM","BH")
write.csv(AllModCost_results_Db.0, file = "ModCost_results_Db.0_Sweden.csv", row.names = TRUE)


#-----------------------------------------------------------------------------
# Run Model with ModCost results and plot fits
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Hakefjord_Head
#-----------------------------------------------------------------------------

fjord <- "Hake"
site  <- "Head"  

# setup boundary conditions

PL$Pb210.supp <- 35     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.115  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 200./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_HH [1] # cm2 yr-1 
PL$x.L     <- ModCost_Result_HH [2]  # cm 
PL$F.Chla  <- 2 
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

# plot and save profiles     


plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("HF_H_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}
#-----------------------------------------------------------------------------
# Hakefjord_Center
#-----------------------------------------------------------------------------

fjord <- "Hake"
site  <- "Center"  

# setup boundary conditions

PL$Pb210.supp <- 48     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.164  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 350./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_HC [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_HC [2]  # cm
PL$F.Chla  <- 2.3 #2.3
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("HF_C_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#Extract parameters to plot Biodiffusion profile of HC

Depth_HC_Db <- PL$grid$x.int
Biodiffusion_HC_Db <- PL$Db.grid$int

#Prepare Plotting Area

x11(height=25,width=40)

#Plot 

plot_HC_Db = plot (y=Depth_HC_Db,x=Biodiffusion_HC_Db,
                   ylim=c(18,0),xlim=c(0,20), pch=16,col="darkgreen",
                   xlab="Bioturbation (yr-1)",ylab="Depth (cm)",axes=T,cex=2.0, cex.axis =3.0,cex.lab =1.5,cex.main=1.5,
                   main="Bioturbation_HC")

savePlot(filename=paste("Biodiffusion_profile_HC.pdf"),type=c("pdf"))
savePlot(filename=paste("Biodiffusion_profile_HC.png"),type=c("png"))

#-----------------------------------------------------------------------------
# Hakefjord_Mouth
#-----------------------------------------------------------------------------

fjord <- "Hake"
site  <- "Mouth"  

# setup boundary conditions

PL$Pb210.supp <- 58     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.470  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 600./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_HM [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_HM [2]  # cm
PL$F.Chla  <- 3
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("HF_M_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#-----------------------------------------------------------------------------
# Gullmarsfjord_Head
#-----------------------------------------------------------------------------


fjord <- "Gullmar"
site  <- "Head"  

# setup boundary conditions

PL$Pb210.supp <- 35     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.067#MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 430./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_GH [1] # cm2 yr-1 
PL$x.L     <- ModCost_Result_GH [2]  # cm 
PL$F.Chla  <- 7
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("GF_H_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#-----------------------------------------------------------------------------
# Gullmarsfjord_Center
#-----------------------------------------------------------------------------

fjord <- "Gullmar"
site  <- "Center"  

# setup boundary conditions

PL$Pb210.supp <- 75     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.150  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 300./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_GC [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_GC [2]  # cm
PL$F.Chla  <- 0.2
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay


# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("GF_C_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#Extract parameters to plot Biodiffusion profile of GC

Depth_GC_Db <- PL$grid$x.int
Biodiffusion_GC_Db <- PL$Db.grid$int

#Prepare Plotting Area

x11(height=25,width=40)

# Plot

plot_GC_Db = plot (y=Depth_GC_Db,x=Biodiffusion_GC_Db,
                   ylim=c(18,0),xlim=c(0,20), pch=16,col="darkgreen",
                   xlab="Bioturbation (yr-1)",ylab="Depth (cm)",axes=T,cex=2.0, cex.axis =3.0,cex.lab =1.5,cex.main=1.5,
                   main="Bioturbation_GC")

savePlot(filename=paste("Biodiffusion_profile_GC.pdf"),type=c("pdf"))
savePlot(filename=paste("Biodiffusion_profile_GC.png"),type=c("png"))

#-----------------------------------------------------------------------------
# Gullmarsfjord_Mouth
#-----------------------------------------------------------------------------

fjord <- "Gullmar"
site  <- "Mouth"  

# setup boundary conditions

PL$Pb210.supp <- 65     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.222  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 450./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_GM [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_GM [2]  # cm
PL$F.Chla  <- 2.5
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("GF_M_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#-----------------------------------------------------------------------------
# Byfjord_Head
#-----------------------------------------------------------------------------

fjord <- "By"
site  <- "Head"  

# setup boundary conditions

PL$Pb210.supp <- 80     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.047  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 100./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_BH [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_BH [2]  # cm
PL$F.Chla  <- 3.5
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay


# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("BF_H_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#-----------------------------------------------------------------------------
# Byfjord_Center
#-----------------------------------------------------------------------------

fjord <- "By"
site  <- "Center"  

# setup boundary conditions

PL$Pb210.supp <- 48    #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.038  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 350./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- 0 # cm2 yr-1
PL$x.L     <- 0  # cm
PL$F.Chla  <- 3.5
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay


# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("BF_C_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}

#Extract parameters to plot Biodiffusion profile of BC

Depth_BC_Db <- PL$grid$x.int
Biodiffusion_BC_Db <- PL$Db.grid$int

#Prepare Plotting Area

x11(height=25,width=40)

#Plot

plot_BC_Db = plot (y=Depth_BC_Db,x=Biodiffusion_BC_Db,
                   ylim=c(18,0),xlim=c(0,20), pch=16,col="darkgreen",
                   xlab="Bioturbation (yr-1)",ylab="Depth (cm)",axes=T,cex=2.0, cex.axis =3.0,cex.lab =1.5,cex.main=1.5,
                   main="Bioturbation_BC")

savePlot(filename=paste("Biodiffusion_profile_BC.pdf"),type=c("pdf"))
savePlot(filename=paste("Biodiffusion_profile_BC.png"),type=c("png"))

#-----------------------------------------------------------------------------
# Byfjord_Mouth
#-----------------------------------------------------------------------------


fjord <- "By"
site  <- "Mouth"  

# setup boundary conditions

PL$Pb210.supp <- 39    #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.084  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 300./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- 0 # cm2 yr-1
PL$x.L     <- 0  # cm
PL$F.Chla  <- 3.5
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay


# plot and save profiles     

plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
       lwd=2)

pdf = {
  pdf("BF_M_Chla_Pb210.pdf")
  plot.function(df=df,fjord=fjord,site=site,PL=PL,Pb210.res=Pb210.res,Chla.res=Chla.res)
  legend("bottomright",bty='n',legend=paste("Db.0=",PL$Db.0,"x.L=",PL$x.L),
         lwd=2)
  dev.off()
}


#-----------------------------------------------------------
#Sensitivity analysis
#-----------------------------------------------------------

#Prepare Plotting Area

x11(height=25,width=40)
par(mfrow=c(1,2))

# Example Hakefjord_Best Fit

fjord <- "Hake"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# setup boundary conditions

PL$Pb210.supp <- 48     #Supported Pb [Bq kg-1]
PL$SedFlux    <- 0.164  #MAR [g cm-2 yr-1] from Watts et al., 2024

PL$F.Pb210    <- 350./PL$conv.Pb210F #flux excess  Pb [Bq m-2 yr-1]

PL$Db.0    <- ModCost_Result_HC [1] # cm2 yr-1
PL$x.L     <- ModCost_Result_HC [2]  # cm
PL$F.Chla  <- 2.3 
PL$Chla.bg <- 0.3

PL$k.1 <- 0.005*365.25

PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)

# Model solution     

output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)

Chla.res  <- output$y[1:PL$N]
Phba.res  <- output$y[(PL$N+1):(2*PL$N)]
rest.res  <- output$y[(2*PL$N+1):(3*PL$N)]
Pb210.res <- output$Pb.decay

BestFit <- Pb210.res+PL$Pb210.supp

#Vary x.L from 1-20 cm

Sens_210Pb_list <- list()
Par_x.L_list <- list()

for (i in 1:20) {
  
  PL$x.L <- 0+1*i
  PL$Db.0 <- ModCost_Result_HC [1]   # cm2 yr-1
  
  PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
  output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
  Pb210.res <- output$Pb.decay
  
  Sens_210Pb_list[[i]] <- paste(Pb210.res+PL$Pb210.supp)  # Add a parameter to the list
  Par_x.L_list[[i]] <- paste(PL$x.L)
}

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == "Pb210") & (df$replicate == 1)

PL$x.L2 <- as.matrix(Par_x.L_list)

x.L_Plot<-plot(x=BestFit,y=PL$grid$x.mid,ylim=c(20,0.),type="l",lwd=4,lty=2,col="black",cex=1.5,xlab=expression("210Pb concentration (Bq kg-1 dwt)"),ylab="Depth (cm)",
                   xlim=c(40,200),cex.axis =1.5, cex.lab =1.5,cex.main=1.5,main = "Pb210-Model Fit_HC")
points(x=df$value[sel.dat],y=df$depth[sel.dat],pch=19,cex=1.5, col="darkblue")

for (i in 1:20) {
  output$y <- as.matrix(Sens_210Pb_list[[i]])
  lines(x=output$y,y=PL$grid$x.mid, col = i + 1,lty = 1, lwd = 2)
}
legend("topright",cex=0.8,legend=paste("x.L =",PL$x.L2), lty = 1, lwd = 2, col = c(2:11),bg="white")


#Vary Db.0 from 1 to 20 cm2 yr-1 

Sens_210Pb_list <- list()
Par_Db.0_list <- list()

for (i in 1:20) {
  
  PL$x.L <- ModCost_Result_HC [2] #cm
  PL$Db.0 <- 0+1*i  
  
  PL <- set.params(df=df,PL=PL,fjord=fjord,site=site)
  output <- steady.1D(y = state, func = Pb.model, parms = PL, nspec = PL$N.var, positive = TRUE)
  Pb210.res <- output$Pb.decay
  
  Sens_210Pb_list[[i]] <- paste(Pb210.res+PL$Pb210.supp)  # Add a parameter to the list
  Par_Db.0_list[[i]] <- paste(PL$Db.0)
}

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == "Pb210") & (df$replicate == 1)

PL$Db.02 <- as.matrix(Par_Db.0_list)

x.L_Plot<-plot(x=BestFit,y=PL$grid$x.mid,ylim=c(20,0.),type="l",lwd=4,lty=2,col="black",cex=1.5,xlab=expression("210Pb concentration (Bq kg-1 dwt)"),ylab="Depth (cm)",
               xlim=c(40,200),cex.axis =1.5, cex.lab =1.5,cex.main=1.5,main = "Pb210-Model Fit_HC")
points(x=df$value[sel.dat],y=df$depth[sel.dat],pch=19,cex=1.5, col="darkblue")

for (i in 1:20) {
  output$y <- as.matrix(Sens_210Pb_list[[i]])
  lines(x=output$y,y=PL$grid$x.mid, col = i + 1,lty = 1, lwd = 2)
}
legend("topright",cex=0.8,legend=paste("Db.0 =",PL$Db.02), lty = 1, lwd = 2, col = c(2:11),bg="white")

# Save Plots

savePlot(filename=paste("Sensitivity_Db_HC_210Pb.pdf"),type=c("pdf"))
savePlot(filename=paste("Sensitivity_Db_HC_210Pb.png"),type=c("png"))
