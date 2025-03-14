## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)
library(EnvStats)
library(readxl)     
library(ggprism)   
library(patchwork) 
library(ggprism) 
library(DescTools)
library(ggExtra)
library(gridExtra)


# Input the PBPK model
source (file = "Swine_OTC_pbpk.R") # Loading the PBPK model code

## Load Model
mod <- mcode_cache("Swine_OTC_PBPK", Swine_OTC_PBPK.code) #refer to mcode function in mrgsolve user guide 3.1.2 Inline


## Input the swine parameters
## Physiological parameters (Please refer to Table 2 for the full name of parameters)
Phy_pars_Swine <- c(
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
  VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412))

## Chemical-specific parameters (Please refer to Table 3 for the full name of parameters)
Chme_pars_OTC <- c(
  PB           =0.709,
  Ka           =0.015,
  KurineC      =0.001,           
  Kint         =0.3,
  KbileC       =0.24,
  Kst          =0.015)   

PC_Swine_OTC <- c(PL = 1.89,
                  PK = 4.75, 
                  PM = 0.851,
                  PF = 0.085,
                  PRest = 4.9)

## Defined the parameters
VarPars_Swine <- c(Chme_pars_OTC, PC_Swine_OTC)   ## Variable parameters 
FixPars_Swine <- c(Phy_pars_Swine, VarPars_Swine) ## Fixed parameters
##---------------------------------------------------------------------
## Model fitting
Data_OTC <- read.csv(file = "Data_OTC.csv")
head(Data_OTC)

## Read the dataset and later used in model calibration
## Study 1: Nielsen et al., 1996; IV; Matrix = serum/ plasma; Dose: 10(iv); Repeat: 1,1
Obs_A1    <- Data_OTC %>% filter(Study == 1)
Obs_A1_IV_P <- Obs_A1 %>% filter(Route == "IV")%>% filter(Matrix == 'P')%>% select(Time = Time, CP = Conc)


## Study 2: FDA document (NADA 008-804 ); PO; K,L,M,F; Dose: 22.05 mg/kg/day; Repeat: 14
Obs_A2    <- Data_OTC %>% filter(Study ==2)
Obs_A2_PO_K  <- Obs_A2 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A2_PO_L  <- Obs_A2 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc)
Obs_A2_PO_F  <- Obs_A2 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc)
Obs_A2_PO_M  <- Obs_A2 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)

### Study 3 :Mercer, 1984;IV;P,M,K;Dose:11 mg/kg/day; Repeat: 1
Obs_A3    <- Data_OTC %>% filter(Study ==3)
Obs_A3_IV_P  <- Obs_A3 %>% filter(Matrix == "P")%>%select(Time = Time, CP = Conc)
Obs_A3_IV_M  <- Obs_A3 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A3_IV_K  <- Obs_A3 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)

### Study 4 :Jacobson et al., 1984; PO; P;Dose: 22mg/kg/day; Repeat: 8
Obs_A4    <- Data_OTC %>% filter(Study ==4)
Obs_A4_PO_P  <- Obs_A4 %>% filter(Matrix == "P")%>%select(Time = Time, CP = Conc)

## Study 5 : Mercer et al., 1978;IV;P,K,L,M,F; Dose: 11 mg/kg/day; Repeat: 1
Obs_A5    <- Data_OTC %>% filter(Study ==5)
Obs_A5_IV_P  <- Obs_A5 %>% filter(Matrix == "P")%>%select(Time = Time, CP = Conc)
Obs_A5_IV_M  <- Obs_A5 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A5_IV_K  <- Obs_A5 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A5_IV_F  <- Obs_A5 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc)
Obs_A5_IV_L  <- Obs_A5 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc)

### Study 6 :FDA document(NADA 130-435); POW; K;Dose: 22 mg/kg/day; Repeat: 5*2; 
Obs_A6    <- Data_OTC %>% filter(Study ==6)
Obs_A6_POW_K  <- Obs_A6 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)

### Study 7 :FDA document(NADA 008-622); POW; K;Dose: 22 mg/kg/day; Repeat: 14*2
Obs_A7   <- Data_OTC %>% filter(Study ==7)
Obs_A7_POW_K  <- Obs_A7 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A7_POW_M  <- Obs_A7 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A7_POW_F  <- Obs_A7 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc)
Obs_A7_POW_L  <- Obs_A7 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc)


pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
  ## Get out of log domain
  parsinput <- exp(Vpars) 
  fixpars["BW"] <- BW
  
  ## Exposure scenarios
  BW          = BW
  tinterval   = tinterval
  TDOSE       = Dtimes
  MW          = 460.439
  DOSE        = Dose*BW/MW  
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_1+ev_2 }
  
  if (route == "oral") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2}
  
  if (route == "ow") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii =  tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2}
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*20, 0.1) 
  
  ## Simulation
  out <- mod %>% param (fixpars) %>% param (parsinput) %>% 
    update(atol = 1E-10, rtol = 1E-5, maxsteps = 50000) %>%
    mrgsim_d (data = ex, tgrid = tsamp)
  
  outdf = cbind.data.frame( Time  = out$time,
                            CP    = out$Plasma*MW,
                            CL    = out$Liver*MW,
                            CK    = out$Kidney*MW,
                            CM    = out$Muscle*MW,
                            CF    = out$Fat*MW)
  return (outdf)
}

## Create the Cost function
## If you run the cost function, you will find the Warning message listed below:
## In regularize.values(x, y, ties, missing(ties)) :collapsing to unique 'x' values
## It casued by the function approx in later R 3.6.0 version; But it doesn't influence the results
Cost_Swine <- function (pars, w) {
  # Prediction
  out_A1_IV  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =30, Dose = 10, Dtimes = 1, route = "iv")
  out_A2_PO  <-pred (fixpars = FixPars_Swine, Vpars = pars, BW =40, Dose = 11.025, Dtimes = 28, tinterval =12,route = "oral")
  out_A5_IV  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =40, Dose = 11, Dtimes = 1, route = "iv")
  out_A4_PO   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =25, Dose = 25, Dtimes = 8, tinterval =12, route = "oral")
  out_A7_POW <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =40, Dose = 11.025, Dtimes = 28,tinterval =12,route = "ow")
  cost<- modCost  (model = out_A1_IV, obs = Obs_A1_IV_P, x ="Time", weight = w)
  cost<- modCost  (model = out_A2_PO, obs = Obs_A2_PO_F, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A2_PO, obs = Obs_A2_PO_K, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A2_PO, obs = Obs_A2_PO_L, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A2_PO, obs = Obs_A2_PO_M, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A4_PO, obs = Obs_A4_PO_P, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A5_IV, obs = Obs_A5_IV_F, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A5_IV, obs = Obs_A5_IV_K, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A5_IV, obs = Obs_A5_IV_L, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A5_IV, obs = Obs_A5_IV_M, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A5_IV, obs = Obs_A5_IV_P, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A7_POW, obs = Obs_A7_POW_F, x ="Time",cost=cost, weight = w)
  cost<- modCost  (model = out_A7_POW, obs = Obs_A7_POW_K, x ="Time",cost=cost, weight = w)
  cost<- modCost  (model = out_A7_POW, obs = Obs_A7_POW_L, x ="Time", cost=cost, weight = w)
  cost<- modCost  (model = out_A7_POW, obs = Obs_A7_POW_M, x ="Time", cost=cost, weight = w)
  return(cost)}

## Use initial parameter to predict concentration
Cost_Swine(log(VarPars_Swine), w="mean")

## Optimize chemical specific parameters
theta_init <- c(
  PB           = 0.709,
  Ka           = 0.02,
  KurineC      = 0.005,           
  Kint         = 0.38,
  KbileC       = 0.25,
  Kst          = 0.015,
  PL           = 1.89,
  PK           = 4.97, 
  PM           = 0.851,
  PF           = 0.2,
  PRest        = 3.0)
#------------------------------------------------------------------------------
# ## Sensitivity function (FME) 
# ## Check the sensitive parameters in the model
# Sns_Swine <- sensFun(func = Cost_Swine, w = "mean",
#                      parms = log(theta_init), varscale = 1)
#  
# Sen_1 <- summary(Sns_Swine)
# # plot(summary(Sns_Swine))
# # ## Selected sensitive parameters;
# theta_Swine <- theta_init[abs(Sen_1$Mean) >1.2*mean(abs(Sen_1$Mean))]
# theta_Swine
# 
theta<- c(
  # PB           = 0.709,
  # Ka           = 0.02,
  # KurineC      = 0.006,
  # Kint       = 0.38,
   KbileC       =0.28,
  # Kst          =0.015,
  # PL = 1.89,
  PK = 4.9736874
  # PM = 0.851,
  # PF = 0.1932379,
  # PRest =2.90
)

# ## Selected sensitive parameters;
Fit_Swine_OTC <- modFit(f = Cost_Swine,
                        p = log(theta),
                        w = "mean",
                        method ="Marq",
                        lower = -Inf, upper = log(theta*100000),
                        control = nls.lm.control(nprint = 1))

summary(Fit_Swine_OTC)

exp(Fit_Swine_OTC$par)

#------------------------------------------------------------------------------

### make a prediction with the optimized parameters 

CSwine_OTC<-Cost_Swine(log(theta_init), w="mean")

## Make a plot for evaluation data
PDat_Swine <- cbind.data.frame (
  OBS = CSwine_OTC$residuals$obs,
  PRE = CSwine_OTC$residuals$mod,
  RES = CSwine_OTC$residuals$res)

PDat_Swine <- PDat_Swine %>% filter(OBS>0)%>%mutate (Log.OBS = log(OBS,10), 
                                                     Log.PRE = log(PRE,10))
# Linear regression analysis
fit_Swine <- lm(Log.PRE ~ Log.OBS, data = PDat_Swine)
summary(fit_Swine)

## Form the calibrated result dataset
PlotDat_Swine <- PDat_Swine %>% mutate(prediction = predict(fit_Swine), 
                                       OPR = PRE/OBS,Class='Calibration')
#Plot
p1 <-
  ggplot(PlotDat_Swine, aes(Log.PRE, Log.OBS)) +
  geom_point  (colour = "red", size = 4)  +
  geom_abline (intercept = 0,
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-3,1), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-3,1),labels = scales::math_format(10^.x))
p1
## Save the fitting results
# saveRDS(CSwine_OTC, file = 'CSwine_OTC.rds')

#///////////////////////////////////////////////////////////////////////////////
## Model evaluation
### Make cost function
Cost_Swine_eva <- function (pars, w) {
  out_A3_IV  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =40, Dose = 11, Dtimes = 1, route = "iv")
  out_A6_POW  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =40, Dose = 11.025, Dtimes = 10,tinterval = 12, route = "ow")
  cost<- modCost  (model = out_A3_IV, obs = Obs_A3_IV_P, x ="Time", weight = w)
  cost<- modCost  (model = out_A3_IV, obs = Obs_A3_IV_K, x ="Time", cost=cost,weight = w)
  cost<- modCost  (model = out_A3_IV, obs = Obs_A3_IV_M, x ="Time", cost=cost,weight = w)
  cost<- modCost(model = out_A6_POW, obs = Obs_A6_POW_K, x ="Time",cost=cost,weight = w)
  return(cost)
}
## The final parameters for model 
Pars_OTC <- c(
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
  PB          = 0.709,
  Ka          = 0.02,
  KurineC     = 0.005,           
  Kint        = 0.38,
  KbileC      = 0.25,
  Kst         = 0.015,
  PL          = 1.89,
  PK          = 4.97, 
  PM          = 0.851,
  PF          = 0.2,
  PRest       = 3.0)

## Predict concentration for evaluation datasets
Cost_Swine_eva (log(Pars_OTC), w="mean")

## Make the data for the plot
CSwine_OTC_eva  <-Cost_Swine_eva(log(Pars_OTC), w="mean")

## Make a plot
PDat_eva   <- cbind.data.frame (OBS = CSwine_OTC_eva$residuals$obs,
                                PRE = CSwine_OTC_eva$residuals$mod,
                                RES = CSwine_OTC_eva$residuals$res)

PDat_eva <- PDat_eva %>% mutate (Log.OBS = log(OBS, 10), 
                                 Log.PRE = log(PRE, 10),)

## linear regression analysis 
fit_eva <- lm(Log.OBS ~ Log.PRE, data = PDat_eva)
summary(fit_eva)

##Plot Linear regression
PlotDat_eva <- PDat_eva %>% mutate(prediction = predict(fit_eva), 
                                   OPR = PRE/OBS,Class='Evaluation')

p2 <-
  ggplot(PlotDat_eva , aes(Log.PRE, Log.OBS)) +
  geom_point  (colour = "blue", size = 4)  +
  geom_abline (intercept = 0,
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-3,1), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-3,1),labels = scales::math_format(10^.x))
p2


# Make a dataframe for total data
Total_output<-rbind(PlotDat_Swine,PlotDat_eva)
# write.csv(Total_output,"Evaluation.csv")

### Estimate the percentage within 2-fold error
n   <- Total_output%>%summarise (count = n())
n_2 <- Total_output%>%filter(OPR>=0.5 & OPR<=2)%>% summarise (count = n())
n_3 <- Total_output%>%filter(OPR>=0.33 & OPR<=3)%>% summarise (count = n())

N2 <- (n_2$count/n$count)
N3 <- (n_3$count/n$count)

print(N2)
print(N3)
# # Save the fitting results
# saveRDS(Pars_OTC,file="parameter.rds")
# saveRDS(CSwine_OTC_eva, file = 'CSwine_OTC_eva.rds')

#Comparision of observed and predicted data
## Form a dataset for observed data
Data_Cal <-rbind.data.frame(
  Obs_A1_1  <- Obs_A1 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_1'),
  Obs_A2_1  <- Obs_A2 %>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_1'),
  Obs_A2_2  <- Obs_A2 %>% filter(Matrix == "L")%>%select(Time = Time,Conc = Conc)%>%mutate(Study = 'A2_2'),
  Obs_A2_3  <- Obs_A2 %>% filter(Matrix == "F")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_3'),
  Obs_A2_4  <- Obs_A2 %>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_4'),
  Obs_A3_1  <- Obs_A3 %>% filter(Matrix == "P")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_1'),
  Obs_A3_2  <- Obs_A3 %>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_2'),
  Obs_A3_3  <- Obs_A3 %>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_3'),
  Obs_A4_1  <- Obs_A4 %>% filter(Matrix == "P")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_1'),
  Obs_A5_1  <- Obs_A5 %>% filter(Matrix == "P")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_1'),
  Obs_A5_2  <- Obs_A5 %>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_2'),
  Obs_A5_3  <- Obs_A5 %>% filter(Matrix == "K")%>% select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_3'),
  Obs_A5_4  <- Obs_A5 %>% filter(Matrix == "F")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_4'),
  Obs_A5_5  <- Obs_A5 %>% filter(Matrix == "L")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_5'),
  Obs_A6_1  <- Obs_A6 %>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_1'),
  Obs_A7_1  <- Obs_A7 %>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_1'),
  Obs_A7_2  <- Obs_A7 %>% filter(Matrix == "L")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_2'),
  Obs_A7_3  <- Obs_A7 %>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_3'),
  Obs_A7_4  <- Obs_A7 %>% filter(Matrix == "F")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_4')
)

# Dose function
preda <- function (Vpars, BW, tinterval = 24, Dose, Dtimes, route,endtime) {
  parsinput<-Vpars
  parsinput["BW"]<-BW
  
  ## Exposure scenarios
  BW          = BW
  tinterval   = tinterval
  TDOSE       = Dtimes
  MW          = 460.439
  DOSE        = Dose*BW/MW  
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_1+ev_2 }
  
  if (route == "oral") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii =tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2}
  
  if (route == "ow") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2}
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, 0.1) 
  
  ## Simulation
  out <- mod %>%param (parsinput) %>% 
    update(atol = 1E-10, rtol = 1E-5, maxsteps = 50000) %>%
    mrgsim_d (data = ex, tgrid = tsamp)
  
  outdf = cbind.data.frame( Time  = out$time,
                            CP    = out$Plasma*MW,
                            CL    = out$Liver*MW,
                            CK    = out$Kidney*MW,
                            CM    = out$Muscle*MW,
                            CF    = out$Fat*MW)
  return (outdf)
}

## Form a dataset of predicted data 
Out_Cal     <- rbind.data.frame(
  out_A1_1    <- preda (Vpars = Pars_OTC, BW = 30, Dose = 10, Dtimes = 1, route = "iv",endtime=24)%>%
    select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_1'),
  out_A2_1    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28, route = "oral",tinterval=12,endtime=500)%>%
    select(Time = Time, Conc = CK)%>%mutate(Study = 'A2_1'),
  out_A2_2    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes =28, route = "oral",tinterval=12,endtime=500)%>%
    select(Time = Time, Conc = CL)%>%mutate(Study = 'A2_2'),
  out_A2_3    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28, route = "oral",tinterval=12,endtime=500)%>%
    select(Time = Time, Conc = CF)%>%mutate(Study = 'A2_3'),
  out_A2_4    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28, route = "oral",tinterval=12,endtime=500)%>%
    select(Time = Time, Conc = CM)%>%mutate(Study = 'A2_4'),
  out_A3_1    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1, route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CP)%>%mutate(Study = 'A3_1'),  
  out_A3_2    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1, route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CM)%>%mutate(Study = 'A3_2'),  
  out_A3_3    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1, route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CK)%>%mutate(Study = 'A3_3'), 
  out_A4_1    <- preda (Vpars = Pars_OTC, BW = 25, Dose = 11, Dtimes = 8, tinterval = 12, route = "oral",endtime=120)%>%
    select(Time = Time, Conc = CP)%>%mutate(Study = 'A4_1'), 
  out_A5_1    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1,  route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CP)%>%mutate(Study = 'A5_1'), 
  out_A5_2    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1,  route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CM)%>%mutate(Study = 'A5_2'), 
  out_A5_3    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1,  route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CK)%>%mutate(Study = 'A5_3'), 
  out_A5_4    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1,  route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CF)%>%mutate(Study = 'A5_4'), 
  out_A5_5    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11, Dtimes = 1,  route = "iv",endtime=32)%>%
    select(Time = Time, Conc = CL)%>%mutate(Study = 'A5_5'),
  out_A6_1    <- preda (Vpars = Pars_OTC, BW = 40,Dose = 11.025, Dtimes = 10,tinterval=12,route = "ow",endtime=216)%>%
    select(Time = Time, Conc = CK)%>%mutate(Study = 'A6_1'),
  out_A7_1    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28,tinterval=12, route = "ow",endtime=500)%>%
    select(Time = Time, Conc = CK)%>%mutate(Study = 'A7_1'),
  out_A7_2    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28, tinterval=12,route = "ow",endtime=500)%>%
    select(Time = Time, Conc = CL)%>%mutate(Study = 'A7_2'),
  out_A7_3    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28,tinterval=12, route = "ow",endtime=500)%>%
    select(Time = Time, Conc = CM)%>%mutate(Study = 'A7_3'),
  out_A7_4    <- preda (Vpars = Pars_OTC, BW = 40, Dose = 11.025, Dtimes = 28,tinterval=12,route = "ow",endtime=500)%>%
    select(Time = Time, Conc = CF)%>%mutate(Study = 'A7_4')
)

###/////////////////////////////////////////////////////////////////////////////
## Plot label
Levels <- c('A1_1', 'A2_1', 'A2_2', 'A2_3',
            'A2_4', 'A3_1', 'A3_2', 'A3_3', 'A4_1', 'A5_1',
            'A5_2', 'A5_3', 'A5_4', 'A5_5', 'A6_1', 'A7_1',
            'A7_2', 'A7_3', 'A7_4')

# Export the datasets 
Data_Cal$Study <- factor(Data_Cal$Study, levels=Levels)
Out_Cal$Study <- factor(Out_Cal$Study, levels=Levels)
# write.csv(Data_Cal, "Observation.csv", row.names = FALSE)
# write.csv(Out_Cal, "Simulation.csv", row.names = FALSE)

#Raw data
df_combined=merge.data.frame(Out_Cal,Data_Cal,by=c("Time","Study"),all = TRUE)
write.csv(df_combined, "Raw data-Comparision Observed and Predicted.csv", row.names = FALSE)

# Determine the Figure 3
p_3<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
  geom_point(shape = 21, colour = "red", fill = "white", size = 1.5, stroke = 1.2) + 
  geom_line(data = data.frame(Out_Cal), aes(x = Time, y = Conc), 
            size = 0.6, colour = "black") +
  
  scale_y_log10 (labels = function(x) format(x, scientific = TRUE))  + 
  facet_wrap(~Study, scales = "free",ncol=6) +
  theme_bw() +
  theme (
    panel.background = element_rect(fill = "white"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    strip.background        = element_blank(),
    strip.text              = element_blank(),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 10, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 10, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "",  y = "")

# ## Export the figure
# ggsave("Fig_2.tiff",scale = 1.5,
#        plot = p_3,
#        width = 25, height = 15, units = "cm", dpi=320)

##Model calibration
## Plot theme
ptheme<-theme (
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect (colour = "black", fill=NA, size=2),
  panel.background        = element_rect (fill="White"),
  panel.grid.major        = element_blank(),
  panel.grid.minor        = element_blank(),
  axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
  axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
  legend.position='none') 

## Plot for model calibration of OTC
p_4 <- 
  ggplot(Total_output, aes(Log.PRE, Log.OBS)) + 
  geom_point  (aes(colour = factor(Class),), size = 4)+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  scale_shape_manual(values = c(15, 16, 17)) +
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", size = 1, alpha = 0.8, linetype = 2) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-3,5), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-3,5),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")

p_4 

#R2 value
summary(fit_Swine)#Calibrated data
summary(fit_eva) # Evaluated data

#Figure for the ratio of predicted value to observed value 

p5 <-
  ggplot(Total_output, aes(Log.PRE, log(OPR,10))) +
  geom_hline(yintercept = log10(2),linetype = 3,color   = "black", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3,color   = "black", size =1) +
  geom_point(aes(colour = factor(Class),),size=2.5) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  scale_shape_manual(values = c(15, 16, 17)) +
  geom_smooth(se = FALSE) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-3,5), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-3,5),labels = scales::math_format(10^.x)) +
  ptheme + labs (x = "", y = "")

ggMarginal(p5, type = "histogram", margins = "y",
           yparams = list(binwidth = 0.1, fill = "#FFA488"))

print(N2)#Within 2 fold error
print(N3)#within 3 fold error
