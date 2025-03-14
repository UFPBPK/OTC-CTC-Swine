## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)

# Input the PBPK model
source (file = "Swine_CTC_pbpk.R") # Loading the PBPK model code

## Load Model
mod <- mcode_cache("Swine_CTC_PBPK", Swine_CTC_PBPK.code) #refer to mcode function in mrgsolve user guide 3.1.2 Inline

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
  VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412),
  PB          = 0.55)

## Chemical-specific parameters (Please refer to Table 3 for the full name of parameters)
Chme_pars_CTC <- c(
  Ka           = 0.08,
  KurineC      = 0.001,           
  Kint         = 0.65,
  KbileC       = 0.45,
  Kst          = 0.008)   

PC_Swine_CTC <- c(  PL = 4.30,
                    PK = 6.06, 
                    PM = 1.7,
                    PF = 1.2,
                    PRest = 1.8)

## Defined the parameters
VarPars_Swine <- c(Chme_pars_CTC, PC_Swine_CTC)   ## Variable parameters 
FixPars_Swine <- c(Phy_pars_Swine, VarPars_Swine) ## Fixed parameters
##---------------------------------------------------------------------
## Model fitting
Data_CTC <- read.csv(file = "Data_CTC.csv")
head(Data_CTC)

### Study 4 :FDA document (NADA 65-480 );POW;K,L,M,F;Dose:22.05 mg/kg/day; Repeat: 45
Obs_A4    <- Data_CTC %>% filter(Study ==4)
Obs_A4_POW_F  <- Obs_A4 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc)
Obs_A4_POW_M  <- Obs_A4 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A4_POW_K  <- Obs_A4 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A4_POW_L  <- Obs_A4 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc)

pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
  ## Get out of log domain
  parsinput <- exp(Vpars) 
  fixpars["BW"] <- BW
  
  ## Exposure scenarios
  BW          = BW
  tinterval   = tinterval
  TDOSE       = Dtimes
  MW          = 478.88
  DOSE        = Dose*BW/MW  
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_1+ev_2 }
  
  if (route == "ow") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
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
  out_A4_POW   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW =20, Dose = 11.025, Dtimes = 90, tinterval =12, route = "ow")
  cost<- modCost  (model = out_A4_POW, obs = Obs_A4_POW_F, x ="Time", weight = w)
  cost<- modCost  (model = out_A4_POW, obs = Obs_A4_POW_K, x ="Time", cost=cost,weight = w)
  cost<- modCost  (model = out_A4_POW, obs = Obs_A4_POW_L, x ="Time", cost=cost,weight = w)
  cost<- modCost  (model = out_A4_POW, obs = Obs_A4_POW_M, x ="Time", cost=cost,weight = w)
  return(cost)}

CSwine_CTC_OW<-Cost_Swine(log(VarPars_Swine), w="mean")

#----------------------------------------------------------------
## plot data of swine for CTC
PDat_S_CTC <- cbind.data.frame (OBS = CSwine_CTC_OW$residuals$obs, 
                                PRE = CSwine_CTC_OW$residuals$mod,
                                RES = CSwine_CTC_OW$residuals$res)

PDat_CTC <- rbind(PDat_S_CTC)
PDat_CTC <- PDat_CTC %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))
fit_CTC <- lm(Log.PRE ~ Log.OBS, data = PDat_CTC)
summary(fit_CTC)

# ## Add the observed-to-prediction ratios column
PlotDat_CTC  <- PDat_CTC %>% mutate(prediction = predict(fit_CTC), 
                                    OPR = PRE/OBS,Class='Calibration')

# Ensure both data frames have columns in the same ordeR
Plotdat1<-rbind(PlotDat_CTC)
Plotdat_2fold<-rbind(Plotdat1)

### Estimate the percentage within 2-fold error
n   <- Plotdat_2fold%>%summarise (count = n())
n_2 <- Plotdat_2fold %>%filter(OPR>=0.5 & OPR<=2)%>% summarise (count = n())
n_3 <- Plotdat_2fold %>%filter(OPR>=0.33 & OPR<=3)%>% summarise (count = n())

N2 <- (n_2$count/n$count)
N3 <- (n_3$count/n$count)

print(N2)
print(N3)

write.csv(PlotDat_CTC, "Evaluation_OW.csv", row.names = TRUE)

a<-read.csv("Evaluation.csv")
b<-read.csv("Evaluation_OW.csv")

Eval_total<-rbind(a,b)

write.csv(Eval_total, "Raw data- GOF-CTC.csv", row.names = FALSE)
#---------------------------------------------------------------------
#Final Parameters
Pars_S_CTC<- c(
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
  Ka          = 0.08,
  KurineC     = 0.001,           
  Kint        = 0.65,
  KbileC      = 0.45,
  Kst         = 0.008,
  PL          = 4.30,
  PK          = 6.06, 
  PM          = 1.70,
  PF          = 1.20,
  PRest       = 1.8)
#---------------------------------------------------------------------
# Adminstration fuction 
preda <- function (Vpars, BW, tinterval = 24, Dose, Dtimes, route,endtime) {
  parsinput<-Vpars
  parsinput["BW"]<-BW
  
  ## Exposure scenarios
  BW          = BW
  tinterval   = tinterval
  TDOSE       = Dtimes
  MW          = 478.88
  DOSE        = Dose*BW/MW  
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = 24, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_1+ev_2 }
  
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

#Import dataset
Data_PG <- read.csv(file = "Data_CTC.csv")
head(Data_PG)
Obs_A4    <- Data_PG %>% filter(Study == 4)
#Observed data
Data_Cal <-rbind.data.frame(
  Obs_A4_1  <- Obs_A4 %>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_1'),
  Obs_A4_2  <- Obs_A4 %>% filter(Matrix == "F")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_2'),
  Obs_A4_3 <- Obs_A4 %>% filter(Matrix == "L")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_3'),
  Obs_A4_4 <- Obs_A4 %>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_4')
)
#Model predicted data
Out_Cal     <- rbind.data.frame(
  out_A4_1    <- preda (Vpars = Pars_S_CTC, BW = 20, Dose = 11.025, Dtimes = 90,tinterval = 12, route = "ow",endtime=1300)%>%
                select(Time = Time, Conc = CK)%>%mutate(Study = 'A4_1'),
  out_A4_2    <- preda (Vpars = Pars_S_CTC, BW = 20, Dose = 11.025, Dtimes = 90,tinterval = 12, route = "ow",endtime=1300)%>%
              select(Time = Time, Conc = CF)%>%mutate(Study = 'A4_2'),
  out_A4_3    <- preda (Vpars = Pars_S_CTC, BW = 20, Dose = 11.025, Dtimes = 90,tinterval = 12,route = "ow",endtime=1300)%>%
              select(Time = Time, Conc = CL)%>%mutate(Study = 'A4_3'),
  out_A4_4    <- preda (Vpars = Pars_S_CTC, BW = 20, Dose = 11.025, Dtimes = 90,tinterval = 12,route = "ow",endtime=1300)%>%
                select(Time = Time, Conc = CM)%>%mutate(Study = 'A4_4')
)
##Figure label 
Levels <- c('A4_1','A4_2', 'A4_3', 'A4_4')

#Graphpad prism data
Data_Cal$Study <- factor(Data_Cal$Study, levels=Levels)# Observed data
Out_Cal$Study <- factor(Out_Cal$Study, levels=Levels)# Predicted data

df_combined_ow=merge.data.frame(Out_Cal,Data_Cal,by=c("Time","Study"),all = TRUE)
write.csv(df_combined_ow, "Raw data-Comparision of observed and predicted-CTC-POW.csv", row.names = FALSE)


#Comparision of observed and simulated data
p<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
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

# # # ## Export the figure
# ggsave("Fig 2-CTC.tiff",scale = 1.5,
#        plot = p,
#        width = 25, height = 15, units = "cm", dpi=320)





