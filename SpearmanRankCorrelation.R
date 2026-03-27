<Spearman rank correlation>
    Copyright (C) <2026>  <Cathrin Wittig>
#=====================================================
# Load data and packages
#=====================================================
setwd("D:/Sweden_2026/R_Scripts/SpearmanRankCorrelation_Figure_7")

library(readr)
library(ggplot2)
library(gridExtra)

Sweden_df <- read_csv("Input_Sweden_Parameters_2026.csv")
View(Sweden_df)

#-----------------------------------------------------------------------------
# Set Parameter Data Frames
#-----------------------------------------------------------------------------

#Fjord Background Data
Depth <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="Depth")]) #Waterdepth in meters                    
Temp <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="Temp")]) #Temperature in degreeCelsius
BW_O2 <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="BW_O2")]) #Bottom water oxygen concentration in µM

#Sediment depositional environment
MAR <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="MAR")]) # Mass Accumulation Rate in g m-2 yr-1

#OM Quality
OC <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="OC")]) #Organic Carbon content in dry wet %
C.N <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="C.N")]) #Organic C/Total N ratio

#Macrofauna community
Abundance <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="Abundance")]) #Abundance in individuals/m2 
Biomass <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="Biomass")]) #Biomass in g/m2
B.A <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="B.A")]) #Biomass/Abundance ratio
BPc <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="BPc")]) #Bioturbation community potential
IPc <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="IPc")]) #Irrigation community potential
BPc_sc <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="BPc_sc")]) #Bioturbation community potential scaled to 1
IPc_sc <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="IPc_sc")]) #Irrigation community potential scaled to 1
BPc.IPc <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="BPc.IPc")]) #Scaled BPc/IPc ratio

#Bioturbation
Db <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="Db")]) #bioturbation rate or particle mixing coefficient in cm2 yr-1                    
xmix <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="xmix")]) #Biomixing depth in cm= mixed layer depth
Irr <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="IrrDIC")]) #Species specific Irrigation rate for DIC in d-1
xirr <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="xirr")]) #Bioirrigation depth in cm=depth attenuation constant for irrigation

#Geochemical Fluxes 

#Burial Fluxes
J_rFe_down <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_rFe_down")]) #reactive Fe burial rate in mmol m-2 d-1     
J_FeS2_down <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_FeS2_down")]) #FeS2 burial rate in mmol m-2 d-1 
J_rMn_down <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_rMn_down")]) #reactive Mn burial rate in mmol m-2 d-1 

#Dissolved upward Fluxes
#Note: Only available in Center stations ! n = 3
#Measured with Lander
J_DIC_up.L <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_DIC_up.L")]) #benthic efflux of dissolved DIC measured by the Lander
J_dMn_up.L <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_dMn_up.L")]) #benthic efflux of dissolved Mn measured by the Lander 
#Modelled with Reactive Transport MOdelling + Irrigation
J_dFe_up.M <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_dFe_up.M")]) #benthic efflux of dissolved Fe modelled
J_dMn_up.M <- as.numeric(Sweden_df$Value[(Sweden_df$Parameter=="J_dMn_up.M")]) #benthic efflux of dissolved Mn modelled
  
#-----------------------------------------------------------------------------
# Calculate Spearmans rank correlation
#-----------------------------------------------------------------------------
#Set Parameters to analyse

#Input Matrix
rows = 9 
cols = 17 
input_matrix <- matrix(c(Depth,Temp,BW_O2,MAR, OC, 
                         C.N,Abundance,Biomass,B.A,BPc,
                         IPc,BPc.IPc,Db,xmix,J_rFe_down,
                         J_FeS2_down,J_rMn_down), nrow=rows, ncol=cols)

#Dependables Matrix
#Statistical relationships tested
Dependables_matrix <- matrix(c(BPc,IPc,Db,xmix,BPc.IPc), nrow=rows, ncol=5)

#Output Matrix
rows = 5
cols = 17

Output_matrix <- matrix(nrow=rows, ncol=cols)
Output_matrix_p <- matrix(nrow=rows, ncol=cols)

#Calculate Spearman rank correlation of all Parameters at all Stations with Db

for (j in 1:rows) {
  x <- Dependables_matrix[,j]

for (i in 1:cols) {
  
  y <- input_matrix[,i]
  
  #Calculate Spearman rank correlation 
  
  spearman_corr <- cor(x, y, method = "spearman")
  correlation <- cor.test(x, y, method = "spearman", exact=FALSE) #Correlation Coefficient of +-0.4 significant, rho is not equal to 0
  p <- correlation$p.value # under-equal 0.05 significant

  #Save results in a matrix
  
  Output_matrix[j,i] <- spearman_corr
  Output_matrix_p[j,i] <- p
}
  
}

#Assign row and column names to output matrix

rownames(Output_matrix) <- c("BPc","IPc","Db","xmix","BPc.IPc")
colnames(Output_matrix) <- c("Depth","Temp","BW_O2","MAR", "OC", "C.N", 
                             "Abundance","Biomass","B.A","BPc","IPc","BPc.IPc","Db","xmix",
                             "J_rFe_down","J_FeS2_down","J_rMn_down")

rownames(Output_matrix_p) <- c("BPc","IPc","Db","xmix","BPc.IPc")
colnames(Output_matrix_p) <- c("Depth","Temp","BW_O2","MAR", "OC", "C.N", 
                               "Abundance","Biomass","B.A","BPc","IPc","BPc.IPc","Db","xmix",
                               "J_rFe_down","J_FeS2_down","J_rMn_down")
print(Output_matrix)
print(Output_matrix_p)


#-----------------------------------------------------------------------------
# Note
#-----------------------------------------------------------------------------
# 
# Strength of Monotonic Relationship
# rs: 0.00 to ±0.19;	±0.20 to ±0.39;	±0.40 to ±0.59;	±0.60 to ±0.79;	±0.80 to ±1.00
# Very weak or no correlation;	Weak correlation;	Moderate correlation;	Strong correlation;	Very strong correlation	
# 
# Direction of Monotonic Relationship
# 1=strong positive monotonic relationship											
# -1=strong negative monotonic relationship																	
# 0=little to no monotonic relationship																	
# 
# Statistical significance															
# p: <0.1, <0.05, <0.01	
# 
#Here: (p≤0.05; rs≥0.60)

#-----------------------------------------------------------------------------
# Save Data
#-----------------------------------------------------------------------------

#Save output in a .csv file

write.csv(Output_matrix, "Spearman_Rank_Corr_Sweden_2026.csv", row.names=TRUE)
write.csv(Output_matrix_p, "Spearman_Rank_p_Sweden_2026.csv", row.names=TRUE)

#-----------------------------------------------------------------------------
# Plot Data
#-----------------------------------------------------------------------------

x11(height=20,width=21)
par(mar = c(4, 3, 2, 3))
par(mfrow=c(3,2))

# Plot identifed negative and positive monotonic relationships (p≤0.05; rs≥0.60; n = 9).

plot_1 = plot (y=IPc,x=Db,
               ylim=c(0,6000),xlim=c(0,25), pch=16,col="black",
               xlab="Db",ylab="Indices",axes=T,cex=2.0, cex.axis =2.0,cex.lab =2.0)
points(x=Db,y=BPc,pch=3,cex=2.0, col="tomato4")
# legend("topright",cex=1.5,legend=c("IPc", "BPc"), pch=c(16,16),col = c("darkblue","darkgreen"),,bg="white")

plot_2 = plot (y=IPc,x=xmix,
               ylim=c(0,6000),xlim=c(0,20), pch=16,col="black",
               xlab="xmix",ylab="Bioturbation",axes=T,cex=2.0, cex.axis =2.0,cex.lab =2.0)# Multiply Db to make it bigger *1000 ?
points(x=xmix,y=Db*200,pch=3,cex=2.0, col="tomato4")
# legend("topright",cex=1.5,legend=c("Db", "IPc"), pch=c(16,16),col = c("darkgreen","darkblue"),,bg="white")

plot_3 = plot (y=BW_O2,x=BPc.IPc,
                   ylim=c(0,300),xlim=c(0,5), pch=16,col="black",
                   xlab="BPc/IPc",ylab="Environment",axes=T,cex=2.0, cex.axis =2.0,cex.lab =2.0)# Multiply OC *50 to make it bigger ?
points(x=BPc.IPc,y=OC*50,pch=3,cex=2.0, col="tomato4")
# legend("topright",cex=1.5,legend=c("BW O2", "OC"), pch=c(16,3),col = c("black","tomato4"),,bg="white")

plot_4 = plot (y=MAR,x=xmix,
               ylim=c(0,5000),xlim=c(0,20), pch=16,col="black",
               xlab="xmix",ylab="Environment",axes=T,cex=2.0, cex.axis =2.0,cex.lab =2.0)# Multiply C.N *100 to make it bigger ?
points(x=xmix,y=C.N*100,pch=3,cex=2.0, col="tomato4")
# legend("topright",cex=1.5,legend=c("MAR", "C/N"), pch=c(16,3),col = c("black","tomato4"),,bg="white")

plot_5 = plot (y=Biomass*10,x=Db,
               ylim=c(0,6000),xlim=c(0,25), pch=16,col="black",
               xlab="Db",ylab="Macrobenthos",axes=T,cex=2.0, cex.axis =2.0,cex.lab =2.0)# Multiply C.N *100 to make it bigger ?
points(x=Db,y=Abundance,pch=3,cex=2.0, col="tomato4")
# legend("topright",cex=1.5,legend=c("Biomass", "Abundance"), pch=c(16,3),col = c("black","tomato4"),,bg="white")

#-----------------------------------------------------------------------------
# Save Plot
#-----------------------------------------------------------------------------

savePlot(filename=paste("SpearmanRankCorrelation.pdf"),type=c("pdf"))
savePlot(filename=paste("SpearmanRankCorrelation.png"),type=c("png"))





