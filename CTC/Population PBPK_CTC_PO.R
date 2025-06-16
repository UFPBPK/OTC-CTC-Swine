## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(patchwork)   # combining multiple ggplot2 plots
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

## Load Model
source (file = "Swine_CTC_pbpk.R") # Loading the generic PBPK model code
mod <- mcode_cache("Swine_CTC_pbpk", Swine_CTC_PBPK.code)

#Individual Data
Pars_CTC <- c(
  BW          = 40,   
  Htc         = 0.412,     
  QCC         = 8.7,    
  QLCa        = 0.243,	                        
  QKCa        = 0.114,	                    
  QMCa        = 0.342,	                       
  QFCa        = 0.128,
  QRestCa     = (1-0.243-0.114-0.342-0.128),
  VLCa        = 0.0204,                                   
  VKCa        = 0.0037,                                  
  VFCa        = 0.1544,                                     
  VMCa        = 0.3632, 
  VbloodCa    = 0.0412,
  VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412),
  PB          = 0.55,
  Ka          = 0.044, #Oral administration via feed ka is 0.044, oral administration via water ka is 0.08
  KurineC     = 0.001,           
  Kint        = 0.65,
  KbileC      = 0.45,
  Kst         = 0.008,
  PL          = 4.30,
  PK          = 6.06, 
  PM          = 1.7,
  PF          = 1.2,
  PRest       = 1.8)

#######-Mean value for each parameter-############
##Physiologcial parameter
#Bodyweight
BW.mean           = 40

#Hematocrit
Htc.mean          = 0.412

#Cardiac output (L/h/kg)
QCC.mean          = 8.7

#Blood Flow (fraction of Cardiac output, unitless)
QLCa.mean         = 0.243	                        
QKCa.mean         = 0.114	                    
QMCa.mean         = 0.342	                       
QFCa.mean         = 0.128
QRestCa.mean      = (1-QLCa.mean-QKCa.mean -QMCa.mean-QFCa.mean )

#Volume of tissues (fraction of body weight, unitless)
VLCa.mean         = 0.0204                                   
VKCa.mean         = 0.0037                                  
VFCa.mean         = 0.1544                                     
VMCa.mean         = 0.3632 
VbloodCa.mean     = 0.0412
VPlas.mean        = (1-Htc.mean) * VbloodCa.mean 
VRestCa.mean      = (1-VLCa.mean -VKCa.mean-VFCa.mean-VMCa.mean-VbloodCa.mean )
#Chemical-specific parameters(mean)
#Fraction of protein binding drug
PB.mean=0.55

#Intestinal absorption rate constant (/h), oral administration with feed
Ka.mean            = 0.044

#Urinary elimination rate constant [L/(h*kg)]
KurineC.mean       = 0.001  

#Bile elimination rate constant [L/(h*kg)]
KbileC.mean        = 0.45

#Intestinal transit rate constant (/h)
Kint.mean          = 0.65

#Gastric empty rate constant (/h)
Kst.mean           = 0.008

#Partition coefficient
PL.mean            = 4.30
PK.mean            = 6.06 
PM.mean            = 1.7
PF.mean            = 1.2
PRest.mean         = 1.8

#######-CV value for each parameter-############
#Physiological parameters
BW.CV           = 0.2

#Hematocrit
Htc.CV          = 0.121

##Cardiac output (L/h/kg)
QCC.CV          = 0.186

#Blood Flow (Fraction of cardiac output)
QLCa.CV         = 0.3	                        
QKCa.CV         = 0.281	                    
QMCa.CV         = 0.895                       
QFCa.CV         = 0.3
QRestCa.CV      = 0.3

#Volume of tissues (Fraction of Bodyweight)
VLCa.CV         = 0.162                                  
VKCa.CV         = 0.297 
VMCa.CV         = 0.073 
VFCa.CV         = 0.172                                     
VbloodCa.CV     = 0.112
VRestCa.CV      = 0.3

#Fraction of protein binding 
PB.CV           = 0.1
## Coefficient of variation for partition coefficients sd/mean 
PC.CV           = 0.2   
## Coefficient of variation for kinetic parameters sd/mean 
K.CV            = 0.3  
## Coefficient of variation for physiological parameters sd/mean 
P.CV            = 0.3                      

########-SD values-###################
BW.sd = BW.CV*BW.mean

#Cardiac output (L/h/kg)
QCC.sd         = 1.6182

## Standard deviation for percentage of drug bound to plasma proteins 
PB.sd = PB.CV*PB.mean 

##Hematocrit
Htc.sd         = 0.05  

## Standard deviation for blood flow fraction
QLCa.sd = P.CV*QLCa.mean             ## Standard deviation for fractions of blood flow to liver 
QKCa.sd = 0.032                      ## Standard deviation for fractions of blood flow to Kidney 
QMCa.sd = 0.306                      ## Standard deviation for fractions of blood flow to Muscle 
QFCa.sd = P.CV*QFCa.mean             ## Standard deviation for fractions of blood flow to fat
QRestCa.sd = P.CV*QRestCa.mean       ## Standard deviation for fractions of blood flow to rest of the body

## Standard deviation for bodyweight 
VLCa.sd = 0.0033                     ## Standard deviation for fractions of volume of liver 
VKCa.sd = 0.0011                     ## Standard deviation for fractions of volume of kidney 
VMCa.sd = 0.0266                     ## Standard deviation for fractions of volume of muscle 
VFCa.sd = 0.0265                     ## Standard deviation for fractions of volume of fat 
VRestCa.sd = P.CV*VRestCa.mean       ## Standard deviation for fractions of volume of rest of the body
VbloodCa.sd= 0.0046                  ## Standard deviation for fractions of volume of blood

## Standard deviation for kinetic parameter
Kst.sd = K.CV*Kst.mean               ## Standard deviation for gastric emptying rate constant
Kint.sd = K.CV*Kint.mean             ## Standard deviation for intestinal transit rate constant
KurineC.sd = K.CV*KurineC.mean       ## Standard deviation for urinary excretion rate constant
KbileC.sd = K.CV*KbileC.mean         ## Standard deviation for bile excretion rate constant
Ka.sd = K.CV*Ka.mean                 ## Standard deviation for intestinal absorption rate constant


## Standard deviation for PC
PL.sd = PC.CV*PL.mean               ## Standard deviation for liver:Plasma Partition coefficient
PK.sd = PC.CV*PK.mean               ## Standard deviation for kidney:Plasma Partition coefficient
PM.sd =PC.CV*PM.mean                ## Standard deviation for muscle:Plasma Partition coefficient
PF.sd = PC.CV*PF.mean               ## Standard deviation for fat:Plasma Partition coefficient
PRest.sd = PC.CV*PRest.mean         ## Standard deviation for rest of the body:Plasma Partition coefficient

#######Calculate log-normal mean and SD#######

#blood flow fraction of Liver to blood
m.log.QLCa= log(QLCa.mean^2/(QLCa.sd^2+QLCa.mean^2)^0.5) 
sd.log.QLCa = (log(1+QLCa.sd^2/QLCa.mean^2))^0.5 

#blood flow fraction of Kidney to blood
m.log.QKCa= log(QKCa.mean^2/(QKCa.sd^2+QKCa.mean^2)^0.5) 
sd.log.QKCa = (log(1+QKCa.sd^2/QKCa.mean^2))^0.5 

#blood flow fraction of Muscle to blood
m.log.QMCa= log(QMCa.mean^2/(QMCa.sd^2+QMCa.mean^2)^0.5) 
sd.log.QMCa = (log(1+QMCa.sd^2/QMCa.mean^2))^0.5 

#blood flow fraction of Fat to blood
m.log.QFCa= log(QFCa.mean^2/(QFCa.sd^2+QFCa.mean^2)^0.5) 
sd.log.QFCa = (log(1+QFCa.sd^2/QFCa.mean^2))^0.5 

#blood flow fraction of Rest to blood
m.log.QRestCa= log(QRestCa.mean^2/(QRestCa.sd^2+QRestCa.mean^2)^0.5) 
sd.log.QRestCa = (log(1+QRestCa.sd^2/QRestCa.mean^2))^0.5 

#PC of liver to plasma
m.log.PL= log(PL.mean^2/(PL.sd^2+PL.mean^2)^0.5) 
sd.log.PL = (log(1+PL.sd^2/PL.mean^2))^0.5 

#PC of Kidney to plasma
m.log.PK = log(PK.mean^2/(PK.sd^2+PK.mean^2)^0.5) 
sd.log.PK = (log(1+PK.sd^2/PK.mean^2))^0.5 

#PC of Muscle to plasma
m.log.PM = log(PM.mean^2/(PM.sd^2+PM.mean^2)^0.5)
sd.log.PM = (log(1+PM.sd^2/PM.mean^2))^0.5

#PC of fat to plasma
m.log.PF= log(PF.mean^2/(PF.sd^2+PF.mean^2)^0.5)
sd.log.PF = (log(1+PF.sd^2/PF.mean^2))^0.5

#PC of rest of body to plasma
m.log.PRest= log(PRest.mean^2/(PRest.sd^2+PRest.mean^2)^0.5)
sd.log.PRest = (log(1+PRest.sd^2/PRest.mean^2))^0.5

#Gastric empty rate constant (/h)
m.log.Kst= log(Kst.mean^2/(Kst.sd^2+Kst.mean^2)^0.5)
sd.log.Kst = (log(1+Kst.sd^2/Kst.mean^2))^0.5

#Intestinal transit rate constant (/h)
m.log.Kint= log(Kint.mean^2/(Kint.sd^2+Kint.mean^2)^0.5)
sd.log.Kint = (log(1+Kint.sd^2/Kint.mean^2))^0.5

#Urinary elimination rate constant [L/(h*kg)]
m.log.KurineC= log(KurineC.mean^2/(KurineC.sd^2+KurineC.mean^2)^0.5) 
sd.log.KurineC = (log(1+KurineC.sd^2/KurineC.mean^2))^0.5

#Bile elimination rate constant [L/(h*kg)]
m.log.KbileC= log(KbileC.mean^2/(KbileC.sd^2+KbileC.mean^2)^0.5)
sd.log.KbileC = (log(1+KbileC.sd^2/KbileC.mean^2))^0.5

#Intestinal absorption rate constant (/h)
m.log.Ka= log(Ka.mean^2/(Ka.sd^2+Ka.mean^2)^0.5) 
sd.log.Ka = (log(1+Ka.sd^2/Ka.mean^2))^0.5 

#Fraction of protein binding 
m.log.PB= log(PB.mean^2/(PB.sd^2+PB.mean^2)^0.5) 
sd.log.PB = (log(1+PB.sd^2/PB.mean^2))^0.5 

##############-Set the distribution-###############
N         <- 1000                 ## Number of iteration
set.seed(42) 

idata_S_CTC<- data_frame(ID=1:N) %>% 
  mutate(
    #Bodyweight  
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = BW.mean, sd = BW.sd), 
                    b = qnorm(0.975, mean = BW.mean, sd = BW.sd), 
                    mean = BW.mean, sd = BW.sd), 
    # Cardiac Output
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QCC.mean, sd = QCC.sd),
      b = qnorm(0.975, mean = QCC.mean, sd = QCC.sd),
      mean = QCC.mean,sd = QCC.sd), 
    # Fraction of protein binding 
    PB = rlnormTrunc(
      N,
      meanlog = m.log.PB,
      sdlog = sd.log.PB,
      min = qlnorm(0.025, meanlog = m.log.PB, sdlog = sd.log.PB),
      max = qlnorm(0.975, meanlog = m.log.PB, sdlog = sd.log.PB)),
    # Hematocrit
    Htc= rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = Htc.mean, sd = Htc.sd), 
                    b = qnorm(0.975, mean = Htc.mean, sd = Htc.sd), 
                    mean = Htc.mean, sd = Htc.sd), 
      
    # Fraction of blood flow to liver 
    QLCa = rlnormTrunc(
      N,
      meanlog=m.log.QLCa,
      sdlog=sd.log.QLCa,
      min = qlnorm(0.025, meanlog=m.log.QLCa, sdlog=sd.log.QLCa),
      max = qlnorm(0.975, meanlog=m.log.QLCa, sdlog=sd.log.QLCa)),
    # Fraction of blood flow to the kidney
    QKCa = rlnormTrunc(
      N,
      meanlog=m.log.QKCa,
      sdlog=sd.log.QKCa,
      min = qlnorm(0.025, meanlog=m.log.QKCa, sdlog=sd.log.QKCa),
      max = qlnorm(0.975, meanlog=m.log.QKCa, sdlog=sd.log.QKCa)),
    # Fraction of blood flow to the muscle 
    QMCa = rlnormTrunc(
      N,
      meanlog=m.log.QMCa,
      sdlog=sd.log.QMCa,
      min = qlnorm(0.025, meanlog=m.log.QMCa, sdlog=sd.log.QMCa),
      max = qlnorm(0.975, meanlog=m.log.QMCa, sdlog=sd.log.QMCa)),
    # Fraction of blood flow to the fat 
    QFCa = rlnormTrunc(
      N,
      meanlog=m.log.QFCa,
      sdlog=sd.log.QFCa,
      min = qlnorm(0.025, meanlog=m.log.QFCa, sdlog=sd.log.QFCa),
      max = qlnorm(0.975, meanlog=m.log.QFCa, sdlog=sd.log.QFCa)),
    # Fraction of blood flow to the rest of body 
    QRestCa = rlnormTrunc(
      N,
      meanlog=m.log.QRestCa,
      sdlog=sd.log.QRestCa,
      min = qlnorm(0.025, meanlog=m.log.QRestCa, sdlog=sd.log.QRestCa),
      max = qlnorm(0.975, meanlog=m.log.QRestCa, sdlog=sd.log.QRestCa)),

    # Fractional organ tissue volumes
    # Fractional liver tissue 
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VLCa.mean, sd =  VLCa.sd),
      b = qnorm(0.975, mean = VLCa.mean, sd =  VLCa.sd) ,
      mean = VLCa.mean,sd =  VLCa.sd),
    # Fractional kidney tissue 
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VKCa.mean, sd =  VKCa.sd),
      b = qnorm(0.975, mean = VKCa.mean, sd =  VKCa.sd),
      mean = VKCa.mean, sd =  VKCa.sd),
    # Fractional fat tissue 
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VFCa.mean, sd =  VFCa.sd),
      b = qnorm(0.975, mean = VFCa.mean, sd =  VFCa.sd),
      mean = VFCa.mean,sd =  VFCa.sd),
    # Fractional muscle tissue 
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VMCa.mean, sd =  VMCa.sd),
      b = qnorm(0.975, mean = VMCa.mean, sd =  VMCa.sd),
      mean = VMCa.mean,sd =  VMCa.sd),
    ##Fraction of blood
    VbloodCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VbloodCa.mean, sd =  VbloodCa.sd),
      b = qnorm(0.975, mean = VbloodCa.mean, sd =  VbloodCa.sd),
      mean = VbloodCa.mean,
      sd =  VbloodCa.sd),
    # Fractional rest of body 
    VRestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VRestCa.mean, sd =  VRestCa.sd),
      b = qnorm(0.975, mean = VRestCa.mean, sd =  VRestCa.sd),
      mean = VRestCa.mean,sd =  VRestCa.sd),
    ###{PC values}
    #Liver: Plasma coefficient
    PL = rlnormTrunc(
      N,
      meanlog = m.log.PL,
      sdlog = sd.log.PL,
      min = qlnorm(0.025, meanlog = m.log.PL, sdlog = sd.log.PL),
      max = qlnorm(0.975, meanlog = m.log.PL, sdlog = sd.log.PL)
    ),
    #Kidney:plasma coefficient
    PK = rlnormTrunc(
      N,
      meanlog = m.log.PK,
      sdlog = sd.log.PK,
      min = qlnorm(0.025, meanlog = m.log.PK, sdlog = sd.log.PK),
      max = qlnorm(0.975, meanlog = m.log.PK, sdlog = sd.log.PK)
    ),
    #Muscle:plasma coefficient
    PM = rlnormTrunc(
      N,
      meanlog = m.log.PM,
      sdlog = sd.log.PM,
      min = qlnorm(0.025, meanlog = m.log.PM, sdlog = sd.log.PM),
      max = qlnorm(0.975, meanlog = m.log.PM, sdlog = sd.log.PM)
    ),
    #Fat:plamsa coefficient
    PF = rlnormTrunc(
      N,
      meanlog = m.log.PF,
      sdlog = sd.log.PF,
      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF),
      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)
    ),
    #Rest:plasma coefficient
    PRest = rlnormTrunc(
      N,
      meanlog = m.log.PRest,
      sdlog = sd.log.PRest,
      min = qlnorm(0.025, meanlog = m.log.PRest, sdlog = sd.log.PRest),
      max = qlnorm(0.975, meanlog = m.log.PRest, sdlog = sd.log.PRest)
    ),
    #{Kinetic Constants}
    ## Gastric emptying rate constant
    Kst = rlnormTrunc(
      N,
      meanlog = m.log.Kst,
      sdlog = sd.log.Kst,
      min = qlnorm(0.025, meanlog = m.log.Kst, sdlog = sd.log.Kst),
      max = qlnorm(0.975, meanlog = m.log.Kst, sdlog = sd.log.Kst)
    ),
    ## Intestinal transit rate constant
    Kint = rlnormTrunc(
      N,
      meanlog = m.log.Kint,
      sdlog = sd.log.Kint,
      min = qlnorm(0.025, meanlog = m.log.Kint, sdlog = sd.log.Kint),
      max = qlnorm(0.975, meanlog = m.log.Kint, sdlog = sd.log.Kint)
    ),
    ## Intestinal absorption rate constant
    Ka = rlnormTrunc(
      N,
      meanlog = m.log.Ka,
      sdlog = sd.log.Ka,
      min = qlnorm(0.025, meanlog = m.log.Ka, sdlog = sd.log.Ka),
      max = qlnorm(0.975, meanlog = m.log.Ka, sdlog = sd.log.Ka)
    ),
    # Bile eliminated Rate constant
    KbileC = rlnormTrunc(
      N,
      meanlog = m.log.KbileC,
      sdlog = sd.log.KbileC,
      min = qlnorm(0.025, meanlog = m.log.KbileC, sdlog = sd.log.KbileC),
      max = qlnorm(0.975, meanlog = m.log.KbileC, sdlog = sd.log.KbileC)
    ),
    # Urinary excreted Rate constants
    KurineC = rlnormTrunc(
      N,
      meanlog = m.log.KurineC,
      sdlog = sd.log.KurineC,
      min = qlnorm(0.025, meanlog = m.log.KurineC, sdlog = sd.log.KurineC),
      max = qlnorm(0.975, meanlog = m.log.KurineC, sdlog = sd.log.KurineC)
    )
  )

## Define the prediction function
pred <- function (pars, drug, idata, tinterval = 24, Dose, Dtimes, route) {
  
  ## Exposure scenarios
  if (is.null(idata$BW) == TRUE) {BW = rep(pars["BW"], N)} else {BW =idata$BW}
  
  
  tinterval   = tinterval
  TDOSE       = Dtimes
  MW          = 478.88
  DOSE        = Dose*BW/MW  
  
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii = tinterval, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    
    ex <- ev_1 }
  
  if (route == "oral") {
    ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii =  tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1:N, amt  = DOSE, ii =  tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2}
  
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*60, 0.5) 
  
  ## Simulation
  out <- mod %>% param(pars)%>%
    update(atol = 1E-10, rtol = 1E-5, maxsteps = 50000) %>%
    mrgsim (idata = idata, ev = ex, tgrid = tsamp)
  
  
  outdf = cbind.data.frame( ID    = out$ID,
                            Time  = out$time,
                            CP    = out$Plasma*MW,
                            CL    = out$Liver*MW,
                            CK    = out$Kidney*MW,
                            CM    = out$Muscle*MW,
                            CF    = out$Fat*MW)
  
  
  return (outdf)
  
}


##--------------Define the plot function-------------------#######

##Plot Theme information    
ptheme<-theme (
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border = element_rect(colour = "black", fill=NA, size=2.5),
  axis.line.x = element_line(size = 2, linetype = "solid", colour = "black"),
  axis.line.y = element_line(size = 2, linetype = "solid", colour = "black"),
  panel.background        = element_rect (fill="White"),
  panel.grid.major        = element_blank(),
  panel.grid.minor        = element_blank(),
  ggh4x.axis.ticks.length.minor = rel(0.5),
  axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
  axis.title              = element_text (size   = 18, colour = "black", face = "bold"),
  axis.ticks.length.x = unit(.25, "cm"),
  axis.ticks = element_line(size = 1),# label of axes
  legend.position='none') 


MCplot <- function(data, target, TOL, tdoses, tinterval = 12,
                   color = "lightsteelblue1", Y.limit = 1e-3) {
  
  # Grouping and summarizing the data
  PDat <- data %>% 
    group_by(Time) %>%
    reframe(
      Time = Time / 24,
      Median = switch(target,
                      "Plasma" = median(CP),
                      "Liver" = median(CL),
                      "Kidney" = median(CK),
                      "Muscle" = median(CM),
                      "Fat" = median(CF)),
      
      Ub = switch(target,
                  "Plasma" = quantile(CP, probs = 0.99),
                  "Liver" = quantile(CL, probs = 0.99),
                  "Kidney" = quantile(CK, probs = 0.99),
                  "Muscle" = quantile(CM, probs = 0.99),
                  "Fat" = quantile(CF, probs = 0.99)),
      
      Lb = switch(target,
                  "Plasma" = quantile(CP, probs = 0.01),
                  "Liver" = quantile(CL, probs = 0.01),
                  "Kidney" = quantile(CK, probs = 0.01),
                  "Muscle" = quantile(CM, probs = 0.01),
                  "Fat" = quantile(CF, probs = 0.01)),
      .groups = 'drop'
    )
  
  # Calculate the x-intercept based on the tolerance and dosing interval
  x.intercept <- PDat %>%filter(Time >= (((tdoses-1)*tinterval/24)))%>% filter (Ub <= TOL) %>%select(Time)%>%min()
  endtime <- x.intercept + (tdoses * tinterval / 24 )+5
  
  # Create the plot
  p <- ggplot(PDat, aes(Time, Median)) + 
    geom_ribbon(aes(ymin = Lb, ymax = Ub), fill = color) +
    geom_line(aes(x = Time, y = Ub), colour = "black", lwd = 1, linetype = 2) +
    geom_line(aes(x = Time, y = Lb), colour = "black", lwd = 1, linetype = 2) +
    geom_line(lwd = 1, colour = "black") +
    scale_x_continuous(
      breaks = seq((tdoses - 1) * tinterval / 24, endtime, by = 2),
      labels = function(x) x - ((tdoses - 1) * tinterval / 24),
      limits = c(0, endtime)
    )+
    scale_y_log10(
      labels = function(x) format(x, scientific = TRUE),
      limits = c(Y.limit, NA)
    ) +
    annotation_logticks(short = unit(2, "mm"),
                        mid = unit(3, "mm"),
                        long = unit(3, "mm"), size = 0.6,
                        sides = "l") +
    labs(x = "", y = "", face = "bold")
  
  # Adding the vertical line for intercept and the tolerance line
  p1 <- p +ptheme+
    geom_vline(aes(xintercept = x.intercept),
               size = 1, color = "red", linetype = 2, show.legend = FALSE) +
    geom_line(aes(y = TOL), color = 'black', size = 1, linetype = 'twodash',
              show.legend = FALSE)
  
  return(p1)
}
MCplot_95 <- function(data, target, TOL, tdoses, tinterval = 12,
                      color = "lightsteelblue1", Y.limit = 1e-3) {
  
  # Grouping and summarizing the data
  PDat <- data %>% 
    group_by(Time) %>%
    reframe(
      Time = Time / 24,
      Median = switch(target,
                      "Plasma" = median(CP),
                      "Liver" = median(CL),
                      "Kidney" = median(CK),
                      "Muscle" = median(CM),
                      "Fat" = median(CF)),
      
      Ub = switch(target,
                  "Plasma" = quantile(CP, probs = 0.95),
                  "Liver" = quantile(CL, probs = 0.95),
                  "Kidney" = quantile(CK, probs = 0.95),
                  "Muscle" = quantile(CM, probs = 0.95),
                  "Fat" = quantile(CF, probs = 0.95)),
      
      Lb = switch(target,
                  "Plasma" = quantile(CP, probs = 0.01),
                  "Liver" = quantile(CL, probs = 0.01),
                  "Kidney" = quantile(CK, probs = 0.01),
                  "Muscle" = quantile(CM, probs = 0.01),
                  "Fat" = quantile(CF, probs = 0.01)),
      .groups = 'drop'
    )
  
  # Calculate the x-intercept based on the tolerance and dosing interval
  x.intercept <- PDat %>%filter(Time >= (((tdoses-1)*tinterval/24)))%>% filter (Ub <= TOL) %>%select(Time)%>%min()
  endtime <- x.intercept + (tdoses * tinterval / 24 )+5
  
  # Create the plot
  p <- ggplot(PDat, aes(Time, Median)) + 
    geom_ribbon(aes(ymin = Lb, ymax = Ub), fill = color) +
    geom_line(aes(x = Time, y = Ub), colour = "black", lwd = 1, linetype = 2) +
    geom_line(aes(x = Time, y = Lb), colour = "black", lwd = 1, linetype = 2) +
    geom_line(lwd = 1, colour = "black") +
    scale_x_continuous(
      breaks = seq((tdoses - 1) * tinterval / 24, endtime, by = 2),
      labels = function(x) x - ((tdoses - 1) * tinterval / 24),
      limits = c(0, endtime)
    )+
    scale_y_log10(
      labels = function(x) format(x, scientific = TRUE),
      limits = c(Y.limit, NA)
    ) +
    annotation_logticks(short = unit(2, "mm"),
                        mid = unit(3, "mm"),
                        long = unit(3, "mm"), size = 0.6,
                        sides = "l") +
    labs(x = "", y = "", face = "bold")
  
  # Adding the vertical line for intercept and the tolerance line
  p1 <- p +ptheme+
    geom_vline(aes(xintercept = x.intercept),
               size = 1, color = "red", linetype = 2, show.legend = FALSE) +
    geom_line(aes(y = TOL), color = 'black', size = 1, linetype = 'twodash',
              show.legend = FALSE)
  
  return(p1)
}


## Estimated the Withdraw intervals
##-------------------Scenario A : Dose 22.05 mg/kg/12h 7 days via oral with feed -----------------#####
Sim_S_A_CTC <-pred (pars = Pars_CTC, idata = idata_S_CTC, 
                    tinterval = 12, Dose = 11.025, 
                    Dtimes = 14, route = 'oral')

CTC_pop_99<-Sim_S_A_CTC %>% group_by (Time) %>% summarise (
  M99 = quantile(CM, probs  = 0.99),
  L99 = quantile(CL, probs  = 0.99),
  K99 = quantile(CK, probs  = 0.99),
  F99= quantile(CF, probs  = 0.99))

CTC_pop_95<-Sim_S_A_CTC %>% group_by (Time) %>% summarise (
  M95 = quantile(CM, probs  = 0.95),
  L95 = quantile(CL, probs  = 0.95),
  K95 = quantile(CK, probs  = 0.95),
  F95= quantile(CF, probs  = 0.95))

## Withdrawal interval determination
#L,K,M,F,USA
TOL_USA<-c(6,12,2,12)

WDIs_CTC_L_USA <- CTC_pop_99 %>% filter(round(L99,3) <= TOL_USA[1] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_USA <- CTC_pop_99 %>% filter(round(K99,3) <= TOL_USA[2] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_USA <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_USA[3] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_F_USA <- CTC_pop_99 %>% filter(round(F99,3) <= TOL_USA[4] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()

#L,K,M,Canada
TOL_Ca<-c(0.6,1.2,0.2)
WDIs_CTC_L_Ca <- CTC_pop_99 %>% filter(round(L99,3) <= TOL_Ca[1] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_Ca <- CTC_pop_99 %>% filter(round(K99,3) <= TOL_Ca[2] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_Ca <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_Ca[3] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()

#L,K,M,F,Japanese
TOL_Ja<-c(0.6,1,0.2,0.2)
WDIs_CTC_L_Ja <- CTC_pop_99 %>% filter(round(L99,3) <= TOL_Ja[1] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_Ja <- CTC_pop_99 %>% filter(round(K99,3) <= TOL_Ja[2] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_Ja <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_Ja[3] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_F_Ja <- CTC_pop_99 %>% filter(round(F99,3) <= TOL_Ja[4] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()

#L,K,M,EU
TOL_EU<-c(0.3,0.6,0.1)
WDIs_CTC_L_EU <- CTC_pop_95 %>% filter(round(L95,3) <= TOL_EU[1] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_EU <- CTC_pop_95 %>% filter(round(K95,3) <= TOL_EU[2] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_EU <- CTC_pop_95 %>% filter(round(M95,3) <= TOL_EU[3] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()

#L,K,M,China
TOL_CH<-c(0.6,1.2,0.2)
WDIs_CTC_L_CH <- CTC_pop_95 %>% filter(round(L95,3) <= TOL_CH[1] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_CH <- CTC_pop_95 %>% filter(round(K95,3) <= TOL_CH[2] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_CH <- CTC_pop_95 %>% filter(round(M95,3) <= TOL_CH[3] & Time >=24*6.5)%>%mutate(WDIs = (Time/24)-6.5)%>%select(WDIs)%>%min()

# Summary
WDIS_TOTAL_USA<-c(WDIs_CTC_L_USA,WDIs_CTC_K_USA,WDIs_CTC_M_USA,WDIs_CTC_F_USA)
WDIS_TOTAL_Ca<-c(WDIs_CTC_L_Ca,WDIs_CTC_K_Ca,WDIs_CTC_M_Ca)
WDIS_TOTAL_Ja<-c(WDIs_CTC_L_Ja,WDIs_CTC_K_Ja,WDIs_CTC_M_Ja,WDIs_CTC_F_Ja)
WDIS_TOTAL_EU<-c(WDIs_CTC_L_EU,WDIs_CTC_K_EU,WDIs_CTC_M_EU)
WDIS_TOTAL_CH<-c(WDIs_CTC_L_CH,WDIs_CTC_K_CH,WDIs_CTC_M_CH)

#Dataframe
data <- ceiling(data.frame(
  liver = c(WDIS_TOTAL_USA[1], WDIS_TOTAL_Ca[1], WDIS_TOTAL_Ja[1], WDIS_TOTAL_EU[1], WDIS_TOTAL_CH[1]),
  kidney = c(WDIS_TOTAL_USA[2], WDIS_TOTAL_Ca[2], WDIS_TOTAL_Ja[2], WDIS_TOTAL_EU[2], WDIS_TOTAL_CH[2]),
  muscle = c(WDIS_TOTAL_USA[3], WDIS_TOTAL_Ca[3], WDIS_TOTAL_Ja[3], WDIS_TOTAL_EU[3], WDIS_TOTAL_CH[3]),
  fat = c(WDIS_TOTAL_USA[4], NA, WDIS_TOTAL_Ja[4], NA, NA)  ))

rownames(data) <- c("USA", "Canada", "Japan", "EU", "China")

print(data)


###Figure for residue depletion for USA
p_CTC_S_A_L_USA<- MCplot(Sim_S_A_CTC , target = "Liver", 
                         TOL =TOL_USA[1] , tdoses = 14, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_USA<- MCplot(Sim_S_A_CTC , target = "Kidney",
                         TOL =TOL_USA[2] , tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_USA<- MCplot(Sim_S_A_CTC , target = "Muscle",
                         TOL =TOL_USA[3], tdoses =14, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
p_CTC_S_A_F_USA<- MCplot(Sim_S_A_CTC , target = "Fat",
                         TOL =TOL_USA[4] , tdoses = 14,tinterval=12, Y.limit = 1e-5, color = 'red')
## Combine the figures 
p1_CTC_USA <-p_CTC_S_A_L_USA+p_CTC_S_A_K_USA+p_CTC_S_A_M_USA +p_CTC_S_A_F_USA

# # ## Save figures for USA
# ggsave("Fig WDT_USA_A.tiff",scale = 1.5,
#        plot = p1_CTC_USA,
#        width = 25, height = 15, units = "cm", dpi=320)

###Figure for depletion for Canada
p_CTC_S_A_L_Ca<- MCplot(Sim_S_A_CTC , target = "Liver", 
                        TOL =TOL_Ca[1] , tdoses = 14, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_Ca<- MCplot(Sim_S_A_CTC , target = "Kidney",
                        TOL =TOL_Ca[2] , tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_Ca<- MCplot(Sim_S_A_CTC , target = "Muscle",
                        TOL =TOL_Ca[3], tdoses =14, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
## Combine the figures
p1_CTC_Ca <-p_CTC_S_A_L_Ca+p_CTC_S_A_K_Ca+p_CTC_S_A_M_Ca

# # ## Save figures
# ggsave("Fig WDT_Canada_A.tiff",scale = 1.5,
#        plot = p1_CTC_Ca,
#        width = 25, height = 15, units = "cm", dpi=320)
# 

###Figure for depletion for Ja 
p_CTC_S_A_L_Ja <- MCplot(Sim_S_A_CTC , target = "Liver", 
                            TOL =TOL_Ja [1] , tdoses = 14, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_Ja <- MCplot(Sim_S_A_CTC , target = "Kidney",
                            TOL =TOL_Ja [2] , tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_Ja <- MCplot(Sim_S_A_CTC , target = "Muscle",
                            TOL =TOL_Ja [3], tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
p_CTC_S_A_F_Ja<- MCplot(Sim_S_A_CTC , target = "Fat",
                           TOL =TOL_Ja[4] , tdoses = 14,tinterval=12, Y.limit = 1e-5, color = 'red')
## Combine the figures
p1_CTC_Ja <-p_CTC_S_A_L_Ja+p_CTC_S_A_K_Ja+p_CTC_S_A_M_Ja+p_CTC_S_A_F_Ja

# # ## Save figures
# ggsave("Fig WDT_Japanese_A.tiff",scale = 1.5,
#        plot = p1_CTC_Ja,
#        width = 25, height = 15, units = "cm", dpi=320)


###Figure for residue depletion for EU 
p_CTC_S_A_L_EU <- MCplot_95(Sim_S_A_CTC , target = "Liver", 
                            TOL =TOL_EU[1] , tdoses = 14, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_EU <- MCplot_95(Sim_S_A_CTC , target = "Kidney",
                            TOL =TOL_EU[2] , tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_EU <- MCplot_95(Sim_S_A_CTC , target = "Muscle",
                            TOL =TOL_EU[3], tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'yellow')  

## Combine the figures
p1_CTC_EU <-p_CTC_S_A_L_EU+p_CTC_S_A_K_EU+p_CTC_S_A_M_EU

# # ## Save figures
# ggsave("Fig WDT_EU_A.tiff",scale = 1.5,
#        plot = p1_CTC_EU,
#        width = 25, height = 15, units = "cm", dpi=320)

###Figure for residue depletion for China 
p_CTC_S_A_L_CH <- MCplot_95(Sim_S_A_CTC , target = "Liver", 
                            TOL =TOL_CH[1] , tdoses = 14, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_CH <- MCplot_95(Sim_S_A_CTC , target = "Kidney",
                            TOL =TOL_CH[2] , tdoses = 14, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_CH <- MCplot_95(Sim_S_A_CTC , target = "Muscle",
                            TOL =TOL_CH[3], tdoses =14, tinterval=12,Y.limit = 1e-5, color = 'yellow')  

p1_CTC_CH <-p_CTC_S_A_L_CH+p_CTC_S_A_K_CH+p_CTC_S_A_M_CH

# # ## Save figures
# ggsave("Fig WDT_CH_A.tiff",scale = 1.5,
#        plot = p1_CTC_CH,
#        width = 25, height = 15, units = "cm", dpi=320)

P_1_SUMMARY<-p_CTC_S_A_L_Ca+p_CTC_S_A_K_Ca+p_CTC_S_A_M_Ca+
            p_CTC_S_A_L_CH+p_CTC_S_A_K_CH+p_CTC_S_A_M_CH+
            p_CTC_S_A_L_EU+p_CTC_S_A_K_EU+p_CTC_S_A_M_EU+
            p_CTC_S_A_L_Ja+p_CTC_S_A_K_Ja+p_CTC_S_A_M_Ja+p_CTC_S_A_F_Ja+
            p_CTC_S_A_L_USA+p_CTC_S_A_K_USA+p_CTC_S_A_M_USA +p_CTC_S_A_F_USA

# ggsave("Figure 22 mg_kg for 7 days via PO-CTC.tiff",scale = 1.5,
#        plot = P_1_SUMMARY,
#        width = 40,height = 20,units = "cm", dpi=320)

## Estimated the Withdraw intervals
##-------------------Scenario A : Dose 22.05 mg/kg/12h 14 days-----------------#####
Sim_S_A_CTC <-pred (pars = Pars_CTC, idata = idata_S_CTC, 
                    tinterval = 12, Dose = 11.025, 
                    Dtimes = 28, route = 'oral')

CTC_pop_99<-Sim_S_A_CTC %>% group_by (Time) %>% summarise (
  M99 = quantile(CM, probs  = 0.99),
  L99 = quantile(CL, probs  = 0.99),
  K99 = quantile(CK, probs  = 0.99),
  F99= quantile(CF, probs  = 0.99))

CTC_pop_95<-Sim_S_A_CTC %>% group_by (Time) %>% summarise (
  M95 = quantile(CM, probs  = 0.95),
  L95 = quantile(CL, probs  = 0.95),
  K95 = quantile(CK, probs  = 0.95),
  F95= quantile(CF, probs  = 0.95))

## Withdrawal interval determination
#L,K,M,F,USA
TOL_USA<-c(6,12,2,12)

WDIs_CTC_L_USA <- CTC_pop_99%>% filter(round(L99,3) <= TOL_USA[1] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_USA<- CTC_pop_99 %>% filter(round(K99,3) <= TOL_USA[2] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_USA <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_USA[3]& Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_F_USA<- CTC_pop_99 %>% filter(round(F99,3) <= TOL_USA[4] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()

#L,K,M,F,Canada
TOL_Ca<-c(0.6,1.2,0.2)
WDIs_CTC_L_Ca <- CTC_pop_99%>% filter(round(L99,3) <= TOL_Ca[1] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_Ca<- CTC_pop_99 %>% filter(round(K99,3) <= TOL_Ca[2] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_Ca <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_Ca[3]& Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()

#L,K,M,F,Japanese
TOL_Ja<-c(0.6,1,0.2,0.2)
WDIs_CTC_L_Ja <- CTC_pop_99%>% filter(round(L99,3) <= TOL_Ja[1] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_Ja<- CTC_pop_99%>% filter(round(K99,3) <= TOL_Ja[2] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_Ja <- CTC_pop_99 %>% filter(round(M99,3) <= TOL_Ja[3]& Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_F_Ja<- CTC_pop_99 %>% filter(round(F99,3) <= TOL_Ja[4] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()

#L,K,M,F,EU
TOL_EU<-c(0.3,0.6,0.1)
WDIs_CTC_L_EU <- CTC_pop_95%>% filter(round(L95,3) <= TOL_EU[1] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_EU<- CTC_pop_95 %>% filter(round(K95,3) <= TOL_EU[2] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_EU <- CTC_pop_95 %>% filter(round(M95,3) <= TOL_EU[3]& Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()

#L,K,M,F,China
TOL_CH<-c(0.6,1.2,0.2)
WDIs_CTC_L_CH <- CTC_pop_95%>% filter(round(L95,3) <= TOL_CH[1] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_K_CH<- CTC_pop_95 %>% filter(round(K95,3) <= TOL_CH[2] & Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()
WDIs_CTC_M_CH <- CTC_pop_95 %>% filter(round(M95,3) <= TOL_CH[3]& Time >=24*13.5)%>%mutate(WDIs = (Time/24)-13.5)%>%select(WDIs)%>%min()

# Summary
WDIS_TOTAL_USA<-c(WDIs_CTC_L_USA,WDIs_CTC_K_USA,WDIs_CTC_M_USA,WDIs_CTC_F_USA)
WDIS_TOTAL_Ca<-c(WDIs_CTC_L_Ca,WDIs_CTC_K_Ca,WDIs_CTC_M_Ca)
WDIS_TOTAL_Ja<-c(WDIs_CTC_L_Ja,WDIs_CTC_K_Ja,WDIs_CTC_M_Ja,WDIs_CTC_F_Ja)
WDIS_TOTAL_EU<-c(WDIs_CTC_L_EU,WDIs_CTC_K_EU,WDIs_CTC_M_EU)
WDIS_TOTAL_CH<-c(WDIs_CTC_L_CH,WDIs_CTC_K_CH,WDIs_CTC_M_CH)

#Dataframe
data <- ceiling(data.frame(
  liver = c(WDIS_TOTAL_USA[1], WDIS_TOTAL_Ca[1], WDIS_TOTAL_Ja[1], WDIS_TOTAL_EU[1], WDIS_TOTAL_CH[1]),
  kidney = c(WDIS_TOTAL_USA[2], WDIS_TOTAL_Ca[2], WDIS_TOTAL_Ja[2], WDIS_TOTAL_EU[2], WDIS_TOTAL_CH[2]),
  muscle = c(WDIS_TOTAL_USA[3], WDIS_TOTAL_Ca[3], WDIS_TOTAL_Ja[3], WDIS_TOTAL_EU[3], WDIS_TOTAL_CH[3]),
  fat = c(WDIS_TOTAL_USA[4], NA, WDIS_TOTAL_Ja[4], NA, NA)  ))

rownames(data) <- c("USA", "Canada", "Japan", "EU", "China")

print(data)


###Figure for residue depletion for USA
p_CTC_S_A_L_USA<- MCplot(Sim_S_A_CTC , target = "Liver", 
                         TOL =TOL_USA[1] , tdoses = 28, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_USA<- MCplot(Sim_S_A_CTC , target = "Kidney",
                         TOL =TOL_USA[2] , tdoses = 28, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_USA<- MCplot(Sim_S_A_CTC , target = "Muscle",
                         TOL =TOL_USA[3], tdoses =28, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
p_CTC_S_A_F_USA<- MCplot(Sim_S_A_CTC , target = "Fat",
                         TOL =TOL_USA[4] , tdoses = 28,tinterval=12, Y.limit = 1e-5, color = 'red')
## Combine the figures 
p1_CTC_USA <-p_CTC_S_A_L_USA+p_CTC_S_A_K_USA+p_CTC_S_A_M_USA +p_CTC_S_A_F_USA

# # ## Save figures for USA
# ggsave("Fig WDT_USA_A.tiff",scale = 1.5,
#        plot = p1_CTC_USA,
#        width = 25, height = 15, units = "cm", dpi=320)

###Figure for depletion for Canada
p_CTC_S_A_L_Ca<- MCplot(Sim_S_A_CTC , target = "Liver", 
                        TOL =TOL_Ca[1] , tdoses = 28, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_Ca<- MCplot(Sim_S_A_CTC , target = "Kidney",
                        TOL =TOL_Ca[2] , tdoses = 28, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_Ca<- MCplot(Sim_S_A_CTC , target = "Muscle",
                        TOL =TOL_Ca[3], tdoses =28, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
## Combine the figures
p1_CTC_Ca <-p_CTC_S_A_L_Ca+p_CTC_S_A_K_Ca+p_CTC_S_A_M_Ca

# # ## Save figures
# ggsave("Fig WDT_Canada_A.tiff",scale = 1.5,
#        plot = p1_CTC_Ca,
#        width = 25, height = 15, units = "cm", dpi=320)
# 

###Figure for depletion for Ja 
p_CTC_S_A_L_Ja <- MCplot(Sim_S_A_CTC , target = "Liver", 
                         TOL =TOL_Ja [1] , tdoses = 28, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_Ja <- MCplot(Sim_S_A_CTC , target = "Kidney",
                         TOL =TOL_Ja [2] , tdoses = 28, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_Ja <- MCplot(Sim_S_A_CTC , target = "Muscle",
                         TOL =TOL_Ja [3], tdoses =28, tinterval=12,Y.limit = 1e-5, color = 'yellow')  
p_CTC_S_A_F_Ja<- MCplot(Sim_S_A_CTC , target = "Fat",
                        TOL =TOL_Ja[4] , tdoses = 28,tinterval=12, Y.limit = 1e-5, color = 'red')
## Combine the figures
p1_CTC_Ja <-p_CTC_S_A_L_Ja+p_CTC_S_A_K_Ja+p_CTC_S_A_M_Ja+p_CTC_S_A_F_Ja

# # ## Save figures
# ggsave("Fig WDT_Japanese_A.tiff",scale = 1.5,
#        plot = p1_CTC_Ja,
#        width = 25, height = 15, units = "cm", dpi=320)


###Figure for residue depletion for EU 
p_CTC_S_A_L_EU <- MCplot_95(Sim_S_A_CTC , target = "Liver", 
                            TOL =TOL_EU[1] , tdoses = 28, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_EU <- MCplot_95(Sim_S_A_CTC , target = "Kidney",
                            TOL =TOL_EU[2] , tdoses = 28, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_EU <- MCplot_95(Sim_S_A_CTC , target = "Muscle",
                            TOL =TOL_EU[3], tdoses =28, tinterval=12,Y.limit = 1e-5, color = 'yellow')  

## Combine the figures
p1_CTC_EU <-p_CTC_S_A_L_EU+p_CTC_S_A_K_EU+p_CTC_S_A_M_EU

# # ## Save figures
# ggsave("Fig WDT_EU_A.tiff",scale = 1.5,
#        plot = p1_CTC_EU,
#        width = 25, height = 15, units = "cm", dpi=320)

###Figure for residue depletion for China 
p_CTC_S_A_L_CH <- MCplot_95(Sim_S_A_CTC , target = "Liver", 
                            TOL =TOL_CH[1] , tdoses = 28, tinterval=12, Y.limit = 1e-5) 
p_CTC_S_A_K_CH <- MCplot_95(Sim_S_A_CTC , target = "Kidney",
                            TOL =TOL_CH[2] , tdoses = 28, tinterval=12,Y.limit = 1e-5, color = 'bisque3')  
p_CTC_S_A_M_CH <- MCplot_95(Sim_S_A_CTC , target = "Muscle",
                            TOL =TOL_CH[3], tdoses =28, tinterval=12,Y.limit = 1e-5, color = 'yellow')  

p1_CTC_CH <-p_CTC_S_A_L_CH+p_CTC_S_A_K_CH+p_CTC_S_A_M_CH

# # ## Save figures
# ggsave("Fig WDT_CH_A.tiff",scale = 1.5,
#        plot = p1_CTC_CH,
#        width = 25, height = 15, units = "cm", dpi=320)

P_1_SUMMARY<-p_CTC_S_A_L_Ca+p_CTC_S_A_K_Ca+p_CTC_S_A_M_Ca+
  p_CTC_S_A_L_CH+p_CTC_S_A_K_CH+p_CTC_S_A_M_CH+
  p_CTC_S_A_L_EU+p_CTC_S_A_K_EU+p_CTC_S_A_M_EU+
  p_CTC_S_A_L_Ja+p_CTC_S_A_K_Ja+p_CTC_S_A_M_Ja+p_CTC_S_A_F_Ja+
  p_CTC_S_A_L_USA+p_CTC_S_A_K_USA+p_CTC_S_A_M_USA +p_CTC_S_A_F_USA

ggsave("Fig p_A.tiff",scale = 1.5,
       plot = P_1_SUMMARY,
       width = 40,height = 20,units = "cm", dpi=320)