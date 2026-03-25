

require(ReacTran)
require(marelac)

#=============================================================================
# Model formulation
#=============================================================================

Pb.model <- function (t,state,parameters,full.output=FALSE) 
{
  with(as.list(c(parameters)),{
    
    # Initialisation of state variables 
    Chla  <- state[1:N]
    Phba  <- state[(N+1):(2*N)]
    rest  <- state[(2*N+1):(3*N)]
    Pb210 <- state[(3*N+1):(4*N)]
    
    # Transport terms
    tran.Chla  <- tran.1D(C=Chla, flux.up=F.Chla,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$dC
    tran.Phba  <- tran.1D(C=Phba, flux.up=F.Phba,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$dC
    tran.rest  <- tran.1D(C=rest, flux.up=F.rest,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$dC
    tran.Pb210 <- tran.1D(C=Pb210,flux.up=F.Pb210,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$dC
    
    # Reaction rates
    R.1      <- svf.grid$mid*k.1*Chla* (Chla>0.)
    R.2      <- svf.grid$mid*k.2*Phba * (Phba>0.)
    Pb.decay <- svf.grid$mid*(log(2)/hl.Pb210)*Pb210 * (Pb210>0.)
    
    # Reaction terms
    Chla.min   <- (-          R.1      )/svf.grid$mid
    Phba.min   <- (+     frac*R.1 - R.2)/svf.grid$mid
    rest.min   <- (+ (1-frac)*R.1 + R.2)/svf.grid$mid
    
    Pb210.min  <- (- Pb.decay)/svf.grid$mid
    
    # Rate of change
    ddt.Chla  <- tran.Chla  + Chla.min
    ddt.Phba  <- tran.Phba  + Phba.min
    ddt.rest  <- tran.rest  + rest.min
    ddt.Pb210 <- tran.Pb210 + Pb210.min
    
    # Mineralization [?mol cm-3 yr-1]
    
    return(list(c(ddt.Chla,ddt.Phba,ddt.rest,ddt.Pb210),
                
                #Reaction rates
                R.1=R.1,
                R.2=R.2,
                Pb.decay=Pb.decay*1e-6*6.022e23/svf.grid$mid/rho.sed*1e3/365.25/24/3600)) # convert to Bq kg-1
    
  })}

#=============================================================================
# auxiliary functions: 
#=============================================================================
#-----------------------------------------------------------------------------
# Function: element budgets 
#-----------------------------------------------------------------------------

IntegratedRate <- function(rate, depth = NULL)  {      # integrated rate for liquids
  if (is.null(depth))
    sum(rate* PL$por.grid$mid * PL$grid$dx)
  else
    sum(rate* PL$por.grid$mid * PL$grid$dx * (PL$grid$x.mid < depth))
} 
IntegratedRateSolid <- function(rate, depth = NULL) {  #                     solids
  if (is.null(depth))
    sum(rate* PL$svf.grid$mid * PL$grid$dx)
  else
    sum(rate* PL$svf.grid$mid * PL$grid$dx * (PL$grid$x.mid < depth))
} 

O2budget <- function(output) {
  
  flux <- output$O2.SWI.flux - output$O2.deep.flux
  
  cons <- IntegratedRate(output$AR) + IntegratedRate(output$CSO)
  
  return(list(Flux = flux, Cons = cons, 
              Delta = flux - cons))
}

#=============================================================================
# Units used in the program: 
#=============================================================================

# Mass =  ?mol
# Space = cm
# Time = yr 

#=============================================================================
# Assumptions
#=============================================================================

# (1) Porosity follows exponential profile with depth 
# (2) Fixed flux of organic matter at the sediment-water interface 
# (3) First order decay of organic matter 

PL <- list()

#=============================================================================
# Model domain and grid definition
#=============================================================================

PL$L    <- 20     # depth of sediment domain [cm]
PL$N    <- 400    # number of grid layers
PL$grid <- setup.grid.1D(x.up = 0,  L =  PL$L, N =  PL$N, dx.1 = PL$L/(10*PL$N))
PL$N.var <- 4 # number of state variables

#=============================================================================
# Model parameters: 
#=============================================================================

# Environmental parameters

PL$S  <- 11    # salinity 
PL$TC <- 6.3    # temperature [deg C]
PL$P  <- 171.*1.013 # pressure [bar]

# Porosity profile 

PL$por.grid     <- setup.prop.1D(value=0.8,grid=PL$grid)
PL$svf.grid     <- setup.prop.1D(value=(1.-0.8),grid=PL$grid)

# Initialisation tortuosity 
# No dissolved species -> not necessary

#PL$tort.grid     <- setup.prop.1D(value=0,grid=PL$grid)
#PL$tort.grid$mid <- (1-2*log(PL$por.grid$mid))
#PL$tort.grid$int <- (1-2*log(PL$por.grid$int))

# Transport parameters 

PL$SedFlux <- 0.1  # sedimentation flux [g cm-2 yr-1]
PL$rho.sed <- 2.6  # density solid sediment [g cm-3]

PL$u.inf <- PL$SedFlux/PL$rho.sed/(1-PL$por.inf)  # sedimentation velocity pore water [cm yr-1]
PL$v.inf <- PL$SedFlux/PL$rho.sed/(1-PL$por.inf)  # sedimentation velocity solids [cm yr-1]

PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v # advection velocity solid
PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u # advection velocity pore water

# Bioturbation profile

PL$x.L    <- 0
PL$x.att  <- 2
PL$y.inf  <- 0
PL$Db.0   <- 0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)

# Reaction parameters 

PL$k.1    <- 0.04*365.25 # d-1 -> yr-1 (Bianchi et al., 2000)
PL$k.2    <- 0.04*365.25
PL$frac   <- 1
PL$K.Chla <- 0.5

PL$hl.Pb210 <- 22.2 # ys

# Boundary conditions 

PL$F.Chla   <- 10*36.525 # flux of organic carbon [?mol C cm-2 yr-1]
PL$F.Phba   <- 0 # flux of organic carbon [?mol C cm-2 yr-1]
PL$F.rest   <- 0 # flux of organic carbon [?mol C cm-2 yr-1]
PL$F.Pb210  <- 4e-11 # flux of Pb210 [?mol C cm-2 yr-1]

# others

PL$conv.Pb210F <- (1e-6*6.022e23*(log(2)/PL$hl.Pb210)/365.25/24/3600*1e4) # convert Bq m-2 yr-1 to model units (moles Pb210 cm-2 yr-1)

#=============================================================================
# Initial conditions:
#=============================================================================

Chla.in  <- rep(0,length.out = PL$N)
Phba.in  <- rep(0,length.out = PL$N)
rest.in  <- rep(0,length.out = PL$N)
Pb210.in <- rep(0,length.out = PL$N)

# Initialization state variables vector 

state <- c(Chla.in,Phba.in,rest.in, Pb210.in)






