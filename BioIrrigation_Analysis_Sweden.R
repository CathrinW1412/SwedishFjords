<Applying FLIPPER (van de Velde, S. (2022). sevdevel/FLIPPER: v0.1.0 (Version v0.1.0.) Zenodo. https://doi.org/10.5281/zenodo.6624982) to pore water data>
    Copyright (C) <2026>  <Cathrin Wittig>
#=====================================================
# Load data and set working directory
#=====================================================
setwd("D:/Sweden_2026/R_Scripts/Bioirrigation_results_Figure 5")

#Load Packages
library(readr)

#Load required FLIPPER files 
#van de Velde, S. (2022). sevdevel/FLIPPER: v0.1.0 (Version v0.1.0.) Zenodo. https://doi.org/10.5281/zenodo.6624982
source("user interface function.R")
source("FLIPPER_plotfunction.R")
source("profile v01.R")

#Load Data to analyse
#Measured Data (Porewater results, sediment porosity)
load("SwedishFjords_Data.Rdata")
#Flux Data
load("SwedishFjords_Fluxes.Rdata")

#=====================================================
# Aux functions
#=====================================================

prepare.input.function <- function(input,var,rep,fjord,site,factor=1.){
  sel     <- (df$fjord==fjord) & (df$site==site) & (df$var == var)& (df$rep == rep)
  sel.por <- (df$fjord==fjord) & (df$site==site) & (df$var == "porosity")
  
  var.temp <- as.data.frame(cbind(df$value[sel]*factor,df$depth[sel]*1e-2))
  var.temp <- var.temp[(var.temp[,2]>0.),]
  por.temp <- cbind(df$value[sel.por],df$depth[sel.por]*1e-2)
  
  max.length  <- nrow(var.temp)#max(nrow(var.temp),nrow(por.temp))
  output.temp <- as.data.frame(cbind(rep(NA,max.length),rep(NA,max.length),rep(NA,max.length)))
  for (i in 1:max.length){
    depth.i         <- var.temp[i,2]
    #por.i           <- por.temp[por.temp[,2]==depth.i,1]
    por.i           <-  por.temp[order(abs(por.temp[,2]-depth.i)),][1,1]
    #por.temp        <- por.temp[-c(por.temp[,2]==depth.i),]
    val.i           <- var.temp[var.temp[,2]==depth.i,1]
    output.temp[i,] <- c(val.i,depth.i,por.i)
  } 
  #  por.temp <- por.temp[,]
  #  if (nrow(por.temp)!=0){
  #  output.temp <- rbind(output.temp,cbind(rep(NA,nrow(por.temp)),por.temp[,2],por.temp[,1]))
  #}
  
  colnames(output.temp) <- c("C","x","por")
  
  output.temp$tort <- 1-2*log(as.numeric(output.temp$por))
  
  return(output.temp)
}


Basic.model  <- function (t=0, C, parms, Prod) {
  with (as.list(parms),{
    
    # if (UBC == "flux.up" & LBC == "no.flux"){
     # tran.summ <- tran.1D(C = C, flux.up = flux.up, C.up = C.up, D = Ds,
     #                      v=v, VF = por.grid, dx = grid, full.output = T)
    # }
    # if (UBC == "conc.up" & LBC == "no.flux"){
      # tran.summ <- tran.1D(C = C, C.up = C.up, D = Ds,v = v,
      #                      VF = por.grid, dx = grid, full.output = T)
    # }
    # if (UBC == "flux.up" & LBC == "conc.down"){
    #   tran.summ <- tran.1D(C = C, flux.up = flux.up, C.up = C.up, C.down = C.down, 
    #                        D = Ds, v=v, VF = por.grid, dx = grid,
    #                        full.output = T)
    # }
    # if (UBC == "conc.up" & LBC == "conc.down"){
    tran.summ <- tran.1D(C = C, C.up = C.up, C.down = C.down,
                         D = Ds, v = v, VF = por.grid, dx = grid,
                         full.output = T)
    # }
    # if (UBC == "flux.up" & LBC == "flux.down"){
    #   tran.summ <- tran.1D(C = C, flux.up = flux.up, flux.down = flux.down,
    #                        D = Ds, v = v, VF = por.grid, 
    #                        dx = grid, full.output = T)
    # }
    # if (UBC == "conc.up" & LBC == "flux.down"){
      # tran.summ <- tran.1D(C = C, C.up = C.up, flux.down = flux.down,
      #                      D = Ds, v = v, VF = por.grid, dx = grid,
      #                      full.output = T)
    # }
    
    # Definition of transport parameters: diffusion + irrigation
    
    tran <- tran.summ$dC # M L-3 T-1 POREWATer
    
    irr <- irr.grid$mid*(C.up - C) # M L-3 T-1 POREWATer
    
    # Return differential equation and other
    
    return(list(dCdt       = tran + irr + Prod/por.grid$mid, # M L-3 T-1 POREWATER
                production = Prod,                           # M L-3 T-1 SEDIMENT
                dif.flux   = tran.summ$dif.flux,
                irr.flux   = sum(grid$dx*irr*por.grid$mid),
                flux.up    = tran.summ$flux.up,
                flux.down  = tran.summ$flux.down,
                adv.flux   = tran.summ$adv.flux
    ))
  })
}

set.parmlist.function <- function(PL,env.parms,input,fjord,site,spec="DIC",R.int){
  
  sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
  
  PL$N    <- 200
  PL$x.up <- min(input$x)
  PL$x.down <- max(input$x[!is.na(input$C)])
  
  PL$grid     <- setup.grid.1D(x.up = 0., x.down = PL$x.down, N = PL$N)
  PL$por.grid <- setup.prop.1D(xy=cbind(input$x,input$por),interpolate="linear", grid = PL$grid)
  PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr,y.inf=0,x.att=PL$irr.att) 
  PL$Ds       <- setup.prop.1D(xy=cbind(input$x,(env.parms$Dmol/input$tort)),interpolate="linear", grid = PL$grid)
  PL$v        <- 0.
  PL$flux.up  <- R.int
  PL$C.up     <- mean(df$value[sel.dat][which(df$depth[sel.dat]==min(df$depth[sel.dat]))])
  #mean(input$C[which(input$x == PL$x.up)])
  PL$C.down   <- mean(input$C[which(input$x == PL$x.down)])
  
  # R.int = - x.att*y.0* (exp(-x.down/x.att) - exp(0))
  # y.0 = R.int / (x.att * (exp(-x.down/x.att) -1.) )
  
  #Prod <- R.int/PL$grid$x.down # mmol m-3 d-1
  PL$prod.0   <-R.int/ (- PL$prod.att*(exp(-PL$grid$x.down/PL$prod.att) - exp(0)))
  PL$prod.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$prod.0,y.inf=0,x.att=PL$prod.att)
  #sum(PL$grid$dx*PL$prod.grid$mid)
  #plot(PL$prod.grid,PL$grid)
  
  return(PL)
}

#---------------------------------------------------
# Modcost_Hakefjord_Center
#---------------------------------------------------

#Hakefjord_Center

fjord <- "Hake"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# DIC

spec  <- "DIC"

sel <- (Flux.list$Flux$code == "HC") & (Flux.list$Flux$var == spec) & (Flux.list$Flux$p.F < 0.05)  
R.int <- mean(Flux.list$Flux$flux[sel])


env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species="HCO3")[["HCO3"]]*3600*24 # [m2 d-1]

input1 <- prepare.input.function(input=df,var=spec,rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var=spec,rep=2,fjord=fjord,site=site)
input <- rbind(input1,input2)

#Prepare Lists and Matrix for results

rows = 10
cols = 10
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_irr_list <- list()
Par_irr.att_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL <- list()
    
    PL$prod.att <- 13./100.
    PL$irr <- 0+0.1*i   # vary Irrigation
    PL$irr.att <- 0.01*j # vary irr.att
    
    #Run Model with current irr and irr.att
    
    PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)
    output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                        Prod = PL$prod.grid$mid, atol = 1e-8)
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
    
    Data <- data.frame (name = c(rep("DIC", length(df$depth[sel.dat]))), depth = df$depth[sel.dat]/100, DIC = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,DIC = output$y)
    
    colnames(out)=c("depth","DIC")
    colnames(Data)=c("name","depth","DIC")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "DIC")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_irr.att_list[[j]] <- paste(PL$irr.att)
  }
  Par_irr_list[[i]] <- paste(PL$irr)
}

#Print Results Matrix
rownames(result_matrix) <- Par_irr_list
colnames(result_matrix) <- Par_irr.att_list
result_matrix

rownames(R2_matrix) <- Par_irr_list
colnames(R2_matrix) <- Par_irr.att_list
R2_matrix

#Find Lowest ModCost
min_value <- min(result_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Irr <- Par_irr_list[[row]]
Irr.att <- Par_irr.att_list[[col]]

# Output

ModCost_Result_HC <- as.numeric(list(Irr,Irr.att,min_value))
ModCost_Result_HC

#---------------------------------------------------
# Modcost_Gullmarsfjord_Center
#---------------------------------------------------

#Gullmarsfjord_Center


fjord <- "Gullmar"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# DIC

spec  <- "DIC"

sel <- (Flux.list$Flux$code == "GC") & (Flux.list$Flux$var == spec) & (Flux.list$Flux$p.F < 0.05)  
R.int <- mean(Flux.list$Flux$flux[sel])
env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species="HCO3")[["HCO3"]]*3600*24 # [m2 d-1]

input1 <- prepare.input.function(input=df,var=spec,rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var=spec,rep=2,fjord=fjord,site=site)
input <- rbind(input1,input2)
input <- input[-30,]

#Prepare Lists and Matrix for results

rows = 10
cols = 10
result_matrix <- matrix(0, nrow = rows, ncol = cols)

Par_irr_list <- list()
Par_irr.att_list <- list()
R2_matrix <- matrix(0, nrow = rows, ncol = cols)

for (i in 1:rows) {
  
  for (j in 1:cols) {
    
    PL <- list()
    
    PL$prod.att <- 13./100.
    PL$irr <- 0+0.1*i   # vary Irrigation
    PL$irr.att <- 0.01*j # vary irr.att
    
    #Run Model with current irr and irr.att
    
    PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)
    output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                        Prod = PL$prod.grid$mid, atol = 1e-8)
    #Do a ModCost Analysis of current run
    #Prepare Data
    sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
    
    Data <- data.frame (name = c(rep("DIC", length(df$depth[sel.dat]))), depth = df$depth[sel.dat]/100, DIC = df$value[sel.dat])
    out <- data.frame(depth = PL$grid$x.mid,DIC = output$y)
    
    colnames(out)=c("depth","DIC")
    colnames(Data)=c("name","depth","DIC")
    
    #ModCost
    
    cost <- modCost(model = out, obs = Data, x = "depth", y = "DIC")
    #Calculate R-squared from fit
    sse <- sum((cost$res[,4] - cost$res[,3])^2)
    ssr <- sum(((cost$res[,4]) - mean(cost$res[,3]))^2)
    sst <- ssr + sse
    
    R2 <- ssr/sst
    
    #Save results in a matrix
    
    result_matrix[i, j] <- cost$model
    R2_matrix[i, j] <- R2
    Par_irr.att_list[[j]] <- paste(PL$irr.att)
  }
  Par_irr_list[[i]] <- paste(PL$irr)
}

#Print Results Matrix
rownames(result_matrix) <- Par_irr_list
colnames(result_matrix) <- Par_irr.att_list
result_matrix

rownames(R2_matrix) <- Par_irr_list
colnames(R2_matrix) <- Par_irr.att_list
R2_matrix

#Find Lowest ModCost
min_value <- min(result_matrix)

# Find the row and column indices of the minimum value
min_indices <- which(result_matrix == min_value, arr.ind = TRUE)
row <- min_indices[1,1]
col <- min_indices[1,2]
Irr <- Par_irr_list[[row]]
Irr.att <- Par_irr.att_list[[col]]

# Output

ModCost_Result_GC <- as.numeric(list(Irr,Irr.att,min_value))
ModCost_Result_GC


#---------------------------------------------------
# Extract All ModCost Results
#---------------------------------------------------

AllModCost_results_Irrigation <- data.frame (ModCost_Result_HC,ModCost_Result_GC)
rownames(AllModCost_results_Irrigation) <- c("Irr","Irr.att","R2")
colnames(AllModCost_results_Irrigation) <- c("HC","GC")

write.csv(AllModCost_results_Irrigation, file = "ModCost_results_Irr_Sweden.csv", row.names = TRUE)


#=====================================================
# Porewater data analysis - DIC
#=====================================================

data.analysis <- list()

#Prepare Plotting Area

x11(height=25,width=40)
par(mfrow=c(2,2))

#-----------------------------------------------------------
# Hake Center
#-----------------------------------------------------------

fjord <- "Hake"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# DIC

spec  <- "DIC"

sel <- (Flux.list$Flux$code == "HC") & (Flux.list$Flux$var == spec) & (Flux.list$Flux$p.F < 0.05)  
R.int <- mean(Flux.list$Flux$flux[sel])

env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species="HCO3")[["HCO3"]]*3600*24 # [m2 d-1]

input1 <- prepare.input.function(input=df,var=spec,rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var=spec,rep=2,fjord=fjord,site=site)
input <- rbind(input1,input2)

PL <- list()

PL$irr.att <- 6/100 # m #value from Best ModcostFit
PL$irr     <- 0.50   # d-1 #value from Best ModcostFit
PL$prod.att <- 13./100. 

PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)

output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                    Prod = PL$prod.grid$mid, atol = 1e-8)

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
sel.dat2 <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)& (df$replicate == 2)

plot(x=output$y,y=PL$grid$x.mid,ylim=c(0.20,0.),type="l",lwd=3,lty=2,col="black",cex=2,xlab=expression("DIC concentration (mmol m"^-3~")"),ylab="Depth (m)",
     xlim=c(min(df$value[sel.dat],na.rm=T),max(df$value[sel.dat],na.rm=T)),cex.axis =1.8, cex.lab =1.8,cex.main=1.8,main = "DIC-Model Fit")
points(x=df$value[sel.dat],y=df$depth[sel.dat]/100,pch=19,cex=2, col="darkblue")
legend("topright",bty='n',cex=1.5,legend=paste("Irr.att =",PL$irr.att,"Irr =",PL$irr))

BestFit <- output$y


# Irrigation vs Depth Plot HC


Depth_HC_Irr <- PL$grid$x.int
Irr_HC <- PL$irr.grid$int

#Plot Irrigation profile of HC

plot_HC_Irr = plot (y=Depth_HC_Irr,x=Irr_HC,
      ylim=c(0.20,0),xlim=c(0,0.7), pch=16,col="darkblue",
      xlab="Irrigation (flux/d)",ylab="Depth (m)",axes=T,cex=1.5, cex.axis =1.8,cex.lab =1.8,cex.main=1.8,
      main="Irrigation")

#-----------------------------------------------------------
# Gullmar Center
#-----------------------------------------------------------

fjord <- "Gullmar"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# DIC

spec  <- "DIC"

sel <- (Flux.list$Flux$code == "GC") & (Flux.list$Flux$var == spec) & (Flux.list$Flux$p.F < 0.05)  
R.int <- mean(Flux.list$Flux$flux[sel])
env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species="HCO3")[["HCO3"]]*3600*24 # [m2 d-1]

input1 <- prepare.input.function(input=df,var=spec,rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var=spec,rep=2,fjord=fjord,site=site)
input <- rbind(input1,input2)
input <- input[-30,]

PL <- list()

PL$irr.att <- 10./100 # m #value from Best ModcostFit
PL$irr     <- 0.4   # d-1 #value from Best ModcostFit
PL$prod.att <- 13./100. 

PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)

output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                    Prod = PL$prod.grid$mid, atol = 1e-8)

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
sel.dat2 <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)& (df$replicate == 2)

plot(x=output$y,y=PL$grid$x.mid,ylim=c(0.25,0.),type="l",lwd=3,lty=2,col="black",cex = 2, xlab=expression("DIC concentration (mmol m"^-3~")"),ylab="Depth (m)",
     xlim=c(min(df$value[sel.dat],na.rm=T),max(df$value[sel.dat],na.rm=T)), cex.axis =1.8,cex.lab =1.8,cex.main=1.8,main = "DIC-Model Fit")
points(x=df$value[sel.dat],y=df$depth[sel.dat]/100,pch=19,cex=2, col="darkblue")
legend("topright",bty='n',cex=1.5,legend=paste("Irr.att =",PL$irr.att,"Irr =",PL$irr))

# Irrigation vs Depth Plot GC

Depth_GC_Irr <- PL$grid$x.int
Irr_GC <- PL$irr.grid$int
plot (y=Depth_GC_Irr,x=Irr_GC,
      ylim=c(0.25,0),xlim=c(0,0.7), pch=16,col="darkblue",
      xlab="Irrigation (flux/d)",ylab="Depth (m)",axes=T,cex=1.5,cex.axis =1.8,cex.lab =1.8,cex.main=1.8,
      main="Irrigation")

#---------------------------------------------------
# Save Plot
#---------------------------------------------------

savePlot(filename=paste("Irrigation_GC_HC_DIC.pdf"),type=c("pdf"))
savePlot(filename=paste("Irrigation_GC_HC_DIC.png"),type=c("png"))

#---------------------------------------------------
# Save File
#---------------------------------------------------
save(data.analysis,file="BioIrrDerivation_Original.Rdata")

#-----------------------------------------------------------
#Sensitivity analysis
#-----------------------------------------------------------

#Prepare Plotting Area

x11(height=25,width=40)
par(mfrow=c(1,2))

# Example Hakefjord

fjord <- "Hake"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

# DIC

spec  <- "DIC"

sel <- (Flux.list$Flux$code == "HC") & (Flux.list$Flux$var == spec) & (Flux.list$Flux$p.F < 0.05)  
R.int <- mean(Flux.list$Flux$flux[sel])

env.parms$Dmol <- diffcoeff(S=env.parms$S,t=env.parms$TC,P=env.parms$P,species="HCO3")[["HCO3"]]*3600*24 # [m2 d-1]

input1 <- prepare.input.function(input=df,var=spec,rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var=spec,rep=2,fjord=fjord,site=site)
input <- rbind(input1,input2)

#Vary irr.att from 0.01-0.1 m
Sens_DIC_list <- list()
Par_irr.att_list <- list()
for (i in 1:10) {
  
PL <- list()

PL$irr.att <- 0.01*i
PL$irr     <- 0.60   # d-1
PL$prod.att <- 13./100. 

PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)

output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                    Prod = PL$prod.grid$mid, atol = 1e-8)

Sens_DIC_list[[i]] <- paste(output$y)  # Add a parameter to the list
Par_irr.att_list[[i]] <- paste(PL$irr.att)
}

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)

PL$irr.att2 <- as.matrix(Par_irr.att_list)
  
Irr.att_Plot<-plot(x=BestFit,y=PL$grid$x.mid,ylim=c(0.20,0.),type="l",lwd=4,lty=2,col="black",cex=1.5,xlab=expression("DIC concentration (mmol m"^-3~")"),ylab="Depth (m)",
     xlim=c(2000,8000),cex.axis =1.5, cex.lab =1.5,cex.main=1.5,main = "DIC-Model Fit_HC")
points(x=df$value[sel.dat],y=df$depth[sel.dat]/100,pch=19,cex=1.5, col="darkblue")
for (i in 1:10) {
  output$y <- as.matrix(Sens_DIC_list[[i]])
  lines(x=output$y,y=PL$grid$x.mid, col = i + 1,lty = 1, lwd = 2)
}
legend("topright",cex=0.8,legend=paste("Irr.att =",PL$irr.att2), lty = 1, lwd = 2, col = c(2:11),bg="white")

#Vary irr from 0-1
Sens_DIC_list <- list()
Par_irr_list <- list()
for (i in 1:10) {
  
  PL <- list()
  
  PL$irr.att <- 0.05
  PL$irr     <- 0+0.1*i   # d-1
  PL$prod.att <- 13./100. 
  
  PL <- set.parmlist.function(PL=PL,env.parms=env.parms,input=input,fjord=fjord,site=site,spec=spec,R.int=R.int)
  
  output <- steady.1D(y = rep(0, PL$grid$N), func = Basic.model, nspec = 1,parms = PL,names = "C", 
                      Prod = PL$prod.grid$mid, atol = 1e-8)
  
  Sens_DIC_list[[i]] <- paste(output$y)  # Add a parameter to the list
  Par_irr_list[[i]] <- paste(PL$irr)
}

sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)

PL$irr2 <- as.matrix(Par_irr_list)

Irr_Plot<-plot(x=BestFit,y=PL$grid$x.mid,ylim=c(0.20,0.),type="l",lwd=4,lty=2,col="black",cex=1.5,xlab=expression("DIC concentration (mmol m"^-3~")"),ylab="Depth (m)",
     xlim=c(2000,8000),cex.axis =1.5, cex.lab =1.5,cex.main=1.5,main = "DIC-Model Fit_HC")
points(x=df$value[sel.dat],y=df$depth[sel.dat]/100,pch=19,cex=1.5, col="darkblue")
for (i in 1:10) {
  output$y <- as.matrix(Sens_DIC_list[[i]])
  lines(x=output$y,y=PL$grid$x.mid, col = i + 1,lty = 1, lwd = 2)
}
legend("topright",cex=0.8,legend=paste("Irr =",PL$irr2), lty = 1, lwd = 2, col = c(2:11),,bg="white")

# Save Plots

savePlot(filename=paste("Sensitivity_Irrigation_HC_DIC.pdf"),type=c("pdf"))
savePlot(filename=paste("Sensitivity_Irrigation_HC_DIC.png"),type=c("png"))





#=====================================================
# Porewater data analysis - Mn and Fe
#=====================================================

rm(list = ls())

data.analysis <- list()


#=====================================================
# Load data
#=====================================================

setwd("C:/Users/cawittig/OneDrive - UGent/UGent-PC/Cawittig/Documents/PhD/Sweden/Final Scripts_CW")
library(readr)
source("user interface function.R")
source("FLIPPER_plotfunction.R")
source("profile v01.R")

load("~/PhD/Sweden/Final Scripts_CW/SwedishFjords_Data.Rdata")
load("~/PhD/Sweden/Final Scripts_CW/SwedishFjords_Fluxes.Rdata")


#=====================================================
# Aux functions
#=====================================================

prepare.input.function <- function(input,var,rep,fjord,site,factor=1.){
  sel     <- (df$fjord==fjord) & (df$site==site) & (df$var == var)& (df$rep == rep)
  sel.por <- (df$fjord==fjord) & (df$site==site) & (df$var == "porosity")
  
  var.temp <- as.data.frame(cbind(df$value[sel]*factor,df$depth[sel]*1e-2))
  var.temp <- var.temp[(var.temp[,2]>0.),]
  por.temp <- cbind(df$value[sel.por],df$depth[sel.por]*1e-2)
  
  max.length  <- nrow(var.temp)#max(nrow(var.temp),nrow(por.temp))
  output.temp <- as.data.frame(cbind(rep(NA,max.length),rep(NA,max.length),rep(NA,max.length)))
  for (i in 1:max.length){
    depth.i         <- var.temp[i,2]
    #por.i           <- por.temp[por.temp[,2]==depth.i,1]
    por.i           <-  por.temp[order(abs(por.temp[,2]-depth.i)),][1,1]
    #por.temp        <- por.temp[-c(por.temp[,2]==depth.i),]
    val.i           <- var.temp[var.temp[,2]==depth.i,1]
    output.temp[i,] <- c(val.i,depth.i,por.i)
  } 
  #  por.temp <- por.temp[,]
  #  if (nrow(por.temp)!=0){
  #  output.temp <- rbind(output.temp,cbind(rep(NA,nrow(por.temp)),por.temp[,2],por.temp[,1]))
  #}
  
  colnames(output.temp) <- c("C","x","por")
  
  output.temp$tort <- 1-2*log(as.numeric(output.temp$por))
  
  return(output.temp)
}

set.parmlist.function <- function(PL,env.parms,input,fjord,site,spec="DIC",R.int){
  
  sel.dat <- (df$fjord==fjord) & (df$site==site) & (df$var == spec)
  
  PL$N    <- 200
  PL$x.up <- min(input$x)
  PL$x.down <- max(input$x[!is.na(input$C)])
  
  PL$grid     <- setup.grid.1D(x.up = 0., x.down = PL$x.down, N = PL$N)
  PL$por.grid <- setup.prop.1D(xy=cbind(input$x,input$por),interpolate="linear", grid = PL$grid)
  PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr,y.inf=0,x.att=PL$irr.att) 
  PL$Ds       <- setup.prop.1D(xy=cbind(input$x,(env.parms$Dmol/input$tort)),interpolate="linear", grid = PL$grid)
  PL$v        <- 0.
  PL$flux.up  <- R.int
  PL$C.up     <- mean(df$value[sel.dat][which(df$depth[sel.dat]==min(df$depth[sel.dat]))])
  #mean(input$C[which(input$x == PL$x.up)])
  PL$C.down   <- mean(input$C[which(input$x == PL$x.down)])
  
  # R.int = - x.att*y.0* (exp(-x.down/x.att) - exp(0))
  # y.0 = R.int / (x.att * (exp(-x.down/x.att) -1.) )
  
  #Prod <- R.int/PL$grid$x.down # mmol m-3 d-1
  PL$prod.0   <-R.int/ (- PL$prod.att*(exp(-PL$grid$x.down/PL$prod.att) - exp(0)))
  PL$prod.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$prod.0,y.inf=0,x.att=PL$prod.att)
  #sum(PL$grid$dx*PL$prod.grid$mid)
  #plot(PL$prod.grid,PL$grid)
  
  return(PL)
}


#=====================================================
# Hakefjord center
#Hakefjord_Center
#=====================================================

fjord <- "Hake"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

Pi_HC <- 1 #Irrigation efficiency acquired by fitting modelled R.int_dMn to measured R.int_dMn

#Mn

input1 <- prepare.input.function(input=df,var="dMn",rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var="dMn",rep=2,fjord=fjord,site=site)
input_HC_Mn <- rbind(input1, input2)

# zones <- 2
discrete.parms <- list(LBC="conc.down", initial.zones=11, irr=0.50*Pi_HC, irr.att=0.06)
data.analysis$Mn$HakeCenter_1_2 <- FLIPPER.func(input=input_HC_Mn,species=c("Mn"),
                                                env.parms=env.parms,discrete.parms=discrete.parms,
                                                method="discrete")


#plot.FLIPPER(data.analysis$Mn$HakeCenter_1_2)
Fit_HC_Mn <- data.analysis$Mn$HakeCenter_1_2$output$fit

#Fe

input1 <- prepare.input.function(input=df,var="dFe",rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var="dFe",rep=2,fjord=fjord,site=site)
input_HC_Fe <- rbind(input1, input2)
input_HC_Fe <- input_HC_Fe[-6,]

# zones <- 5
discrete.parms <- list(LBC="conc.down", initial.zones=8, irr=0.50*Pi_HC, irr.att=0.06)
data.analysis$Fe$HakeCenter_1_2 <- FLIPPER.func(input=input_HC_Fe,species=c("Fe"),
                                                env.parms=env.parms,discrete.parms=discrete.parms,
                                                method="discrete")


#plot.FLIPPER(data.analysis$Fe$HakeCenter_1_2)
Fit_HC_Fe <- data.analysis$Fe$HakeCenter_1_2$output$fit

#=====================================================
# Gullmarfjord center
#Gullmarsfjord_Center
#=====================================================

fjord <- "Gullmar"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

Pi_GC <- 0.25 #Irrigation efficiency acquired by fitting modelled R.int_dMn to measured R.int_dMn

#Mn

input1 <- prepare.input.function(input=df,var="dMn",rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var="dMn",rep=2,fjord=fjord,site=site)
input_GC_Mn <- rbind(input1, input2)

# zones <- 2
discrete.parms <- list(LBC="conc.down", initial.zones=10, irr=0.4*Pi_GC, irr.att=0.1) #reducing irr_DIC = 0.4 to Irr_Mn = 0.105 due to oxidative behavior to fit the measured R.int of dMn
data.analysis$Mn$GullmarCenter_1_2 <- FLIPPER.func(input=input_GC_Mn,species=c("Mn"),
                                                   env.parms=env.parms,discrete.parms=discrete.parms,
                                                   method="discrete")

#plot.FLIPPER(data.analysis$Mn$GullmarCenter_1_2)
Fit_GC_Mn <- data.analysis$Mn$GullmarCenter_1_2$output$fit

#Fe

fjord <- "Gullmar"
site  <- "Center"

env.parms <- list("TC"=df$value[(df$var=="T") & (df$fjord==fjord) & (df$site==site)],
                  "S"=df$value[(df$var=="S") & (df$fjord==fjord) & (df$site==site)],
                  "P"=(df$value[(df$var=="water.depth") & (df$fjord==fjord) & (df$site==site)]/10. + 1.)*1.1013)

input1 <- prepare.input.function(input=df,var="dFe",rep=1,fjord=fjord,site=site)
input2 <- prepare.input.function(input=df,var="dFe",rep=2,fjord=fjord,site=site)
input_GC_Fe <- rbind(input1, input2)

# zones <- 4
discrete.parms <- list(LBC="no.flux", initial.zones=10, irr=0.4*Pi_GC, irr.att=0.1) #reducing irr_DIC = 0.4 to Irr_Mn = 0.105, Irr_Mn = Irr_Fe assumed due to lack of measured flux data
data.analysis$Fe$GullmarCenter_1_2 <- FLIPPER.func(input=input_GC_Fe,species=c("Fe"),
                                                   env.parms=env.parms,discrete.parms=discrete.parms,
                                                   method="discrete")

#plot.FLIPPER(data.analysis$Fe$GullmarCenter_1_2)
Fit_GC_Fe <- data.analysis$Fe$GullmarCenter_1_2$output$fit

#=====================================================
# Plot Mn and Fe Model Fits
#=====================================================

#Prepare Plotting Area

x11(height=25,width=40)
par(mfrow=c(2,2))

plot (x=input_HC_Mn$C,y=input_HC_Mn$x,
      ylim=c(0.20,0),xlim=c(0,100), pch=19,col="violetred",
      xlab=expression("Mn concentration (µmol m"^-3~")"),ylab="Depth (m)",axes=T,cex=2,cex.axis =1.8, cex.lab =1.8,cex.main=1.8,
      main="Mn-Model Fit")
lines (x=Fit_HC_Mn$C, y=Fit_HC_Mn$x, lwd=3, lty=2, col="black")

plot (x=input_HC_Fe$C,y=input_HC_Fe$x,
      ylim=c(0.20,0),xlim=c(0,100), pch=19,col="darkblue",
      xlab=expression("Fe concentration (µmol m"^-3~")"),ylab="Depth (m)",axes=T,cex=2,cex.axis =1.8, cex.lab =1.8,cex.main=1.8,
      main="Fe-Model Fit")
lines (x=Fit_HC_Fe$C, y=Fit_HC_Fe$x, lwd=3, lty=2, col="black")

plot (x=input_GC_Mn$C,y=input_GC_Mn$x,
      ylim=c(0.25,0),xlim=c(0,300), pch=19,col="violetred",
      xlab=expression("Mn concentration (µmol m"^-3~")"),ylab="Depth (m)",axes=T,cex=2,cex.axis =1.8, cex.lab =1.8,cex.main=1.8,
      main="Mn-Model Fit")
lines (x=Fit_GC_Mn$C, y=Fit_GC_Mn$x, lwd=3, lty=2, col="black")

plot (x=input_GC_Fe$C,y=input_GC_Fe$x,
      ylim=c(0.25,0),xlim=c(0,300), pch=19,col="darkblue",
      xlab=expression("Fe concentration (µmol m"^-3~")"),ylab="Depth (m)",axes=T,cex=2,cex.axis =1.8, cex.lab =1.8,cex.main=1.8,
      main="Fe-Model Fit")
lines (x=Fit_GC_Fe$C, y=Fit_GC_Fe$x, lwd=3, lty=2, col="black")

#Save Plots
savePlot(filename=paste("Irrigation_GC_HC_Mn_Fe.pdf"),type=c("pdf"))
savePlot(filename=paste("Irrigation_GC_HC_Mn_Fe.png"),type=c("png"))
