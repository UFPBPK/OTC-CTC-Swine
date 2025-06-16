
##---------------------------------------------------------------------------------------------------
# This mrgsolve PBPK model code is developed for the generic model 
# to simulate the concentrations of Oxytetracycline
# for swine
#----------------------------------------------------------------------------------------------------

Swine_OTC_PBPK.code <- '
$PROB @annotated
# PBPK model for OTC
- Author    : Kun Mi
- Adivisor  : Zhoumeng Lin
- Date      : June, 2024
- Strucutre : GI tract, Muscle, rest of body, kidneys, liver, fat, venous and arterial plasma
- Default physiological parameters value obtained from Lin et al. (2020)

$PARAM @annotated
BW                  : 40       : kg,                  Body weight
Htc                 : 0.412    : Unitless             Hematocrit for swine
QCC                 : 8.7      : L/h/kg,              Cardiac output (L/h/kg) (Lin et al., 2020, Table 36)
QLCa                : 0.243    : Unitless,            Fraction of blood flow to the liver (Lin et al., 2020, Table 28), hepatic artery + portal vein
QKCa                : 0.114    : Unitless,            Fraction of blood flow to the kidneys (Lin et al., 2020., Table 28)
QMCa                : 0.342    : Unitless,            Fraction of blood flow to the muscle (Lin et al., 2020, Table 28)
QFCa                : 0.128    : Unitless,            Fraction of blood flow to the fat (Li et al., 2017, Table 2)
QRestCa             : 0.173    :                      Fraction of blood flow to the rest of body, Value from (total sum equals to 1) (1-QLC-QKC-QFC-QMC)
VLCa                : 0.0204   : L/kg BW,             Fractional liver tissue (Lin et al., 2020, Table 9)
VKCa                : 0.0037   : L/kg BW,             Fractional kidney tissue (Lin et al., 2020, Table 9)
VMCa                : 0.3632   : L/kg BW,             Fractional muscle tissue (Lin et al., 2020, Table 9)
VFCa                : 0.1544   : L/kg BW,             Fractional fat  (Lin et al., 2020, Table 9)
VbloodCa            : 0.041    : L/kg BW,             Fractional blood (Lin et al., 2020, Table 9)
VRestCa             : 0.3654   : L/kg BW,             Fractional rest of body, Value calculated from 1-VLC-VKC-VFC-VMC-VbloodC
PL                  : 1.8      : Unitless,            Liver:plasma PC 
PK                  : 4.5      : Unitless,            Kidney:plasma PC 
PM                  : 0.85     : Unitless,            Muscle:plasma PC 
PF                  : 0.085    : Unitless,            Fat:plasma PC 
PRest               : 0.851    : Unitless,            Rest:plasma PC 
PB                  : 0.709    : Unitless,            Fraction of protein binding drug in plasma 
Ka                  : 0.02     : 1/h/kg ,             Rate of absorption of drug in small intestines
KurineC             : 0.001    : 1/h/kg,              Urinary elimination rate constant
KbileC              : 0.02     : 1/h/kg,              Bile elimination rate constant
Kint                : 0.38     : 1/h/kg,              Intestine transit rate constant
Kst                 : 0.015    : 1/h,                 Gastric emptying rate constant
MW                  : 460.439  : g/mol,               Molecular mass MW as default value 

$MAIN
double sumQ        = QLCa + QKCa + QMCa + QFCa + QRestCa;             // Sum up cardiac output fraction
double sumV        = VLCa + VKCa + VMCa + VbloodCa + VFCa + VRestCa;  // Sum up the tissue volumes
double QLC         = QLCa/sumQ;                                       // Adjusted blood flow rate fraction to liver
double QKC         = QKCa/sumQ;                                       // Adjusted blood flow rate fraction to kidney
double QMC         = QMCa/sumQ;                                       // Adjusted blood flow rate fraction to muscle
double QFC         = QFCa/sumQ;                                       // Adjusted blood flow rate fraction to fat
double QRestC      = 1-QLC-QKC-QMC-QFC;                               // Adjusted blood flow rate fraction to rest of body
double VLC         = VLCa/sumV;                                       // Adjusted fraction of tissue volume of liver
double VKC         = VKCa/sumV;                                       // Adjusted fraction of tissue volume of kidney
double VMC         = VMCa/sumV;                                       // Adjusted fraction of tissue volume of muscle
double VFC         = VFCa/sumV;                                       // Adjusted fraction of tissue volume of fat
double VbloodC     = VbloodCa/sumV;                                   // Adjusted fraction of tissue volume of blood
double VRestC      = 1-VLC-VbloodC-VKC-VMC-VFC;                       // Adjusted fraction of tissue volume of rest of body
double QC          = QCC*BW;                                          // Cardiac output
double QL          = QLC*QC;                                          // Blood flow in Liver
double QK          = QKC*QC;                                          // Blood flow in Kidney
double QM          = QMC*QC;                                          // Blood flow in Muscle
double QF          = QFC*QC;                                          // Blood flow in Fat
double QRest       = QRestC*QC;                                       // Blood flow in Rest of body 
double VL          = VLC*BW;                                          // Tissue volume of Liver
double VK          = VKC*BW;                                          // Tissue volume of Kidney
double VM          = VMC*BW;                                          // Tissue volume of Muscle
double VF          = VFC*BW;                                          // Tissue volume of Fat
double VRest       = VRestC*BW;                                       // Tissue volume of Rest of body 
double Vblood      = VbloodC*BW;                                      // Volume of Blood
double VPlas       = Vblood*(1-Htc);                                  // Volume of Plasma
double Kurine      = KurineC*BW;                                      // Urine elimination rate 
double Kbile       = KbileC*BW;                                       // Bile elimination rate 
double Free        = 1-PB;                                            // Fraction of free drug in plasma
//---------------------------------------------------------------------------------------
$INIT @annotated
ADOSE           : 0   : mg, Amount of input dose; virtual compartment
APlas_free      : 0   : mg, Amount of unbound drug in the plasma compartment  
AST             : 0   : mg, Amount of drug in the stomach
AabsI           : 0   : mg, Amount of drug absorbed by small intestine 
Afeces          : 0   : mg, Amount of drug excretion by feces  
ASI             : 0   : mg, Amount of drug in the small intestine
Abile           : 0   : mg, Amount of drug through biliary secretion
AL              : 0   : mg, Amount of drug in the liver compartment    
AM              : 0   : mg, Amount of drug in the muscle compartment  
AK              : 0   : mg, Amount of drug in the kidney compartment   
Aurine          : 0   : mg, Amount of drug excreted by the urine   
AF              : 0   : mg, Amount of drug in the fat compartment
ARest           : 0   : mg, Amount of drug in the rest of body compartment
AUCCP           : 0   : mg/L*hr, Area under curve of drug in the plasma
AUCCL           : 0   : mg/L*hr, Area under curve of drug in the liver
AUCCK           : 0   : mg/L*hr, Area under curve of drug in the kidney
AUCCM           : 0   : mg/L*hr, Area under curve of drug in the muscle
AUCCF           : 0   : mg/L*hr, Area under curve of drug in the fat
//----------------------------------------------------------------------------------------
$ODE
// Concentrations in the tissues and in the venous plasma leaving each of the tissues
double CPlas_free            = APlas_free/VPlas;      
double CPlas                 = CPlas_free/Free;                           
double CL                    = AL/VL;                                   
double CVL                   = CL/PL;                                  
double CF                    = AF/VF;  
double CVF                   = CF/PF;  
double CK                    = AK/VK;
double CVK                   = CK/PK;
double CM                    = AM/VM;
double CVM                   = CM/PM;
double CRest                 = ARest/VRest;                          
double CVRest                = CRest/PRest;
double CV                    = (QL*CVL + QK*CVK + QF*CVF + QM*CVM + QRest*CVRest)/QC; 

//-- The changing rate of the amount of OTC in plasma 
double RPlas_free   = QC*(CV-CPlas)*Free ;
//-- The amount and AUC of OTC in plasma  
dxdt_APlas_free     = RPlas_free;

// OTC Concentration in gut lumen
double RST                   = -Kst*AST;   
double RabsI                 = Ka*ASI;   
double Rfeces                = Kint*ASI;   
double RSI                   = Kst*AST - Ka*ASI +Kbile*CVL- Kint*ASI;    
double Rbile                 = Kbile*CVL;

//-- The amount and AUC of OTC in lumen
dxdt_AST              = RST;
dxdt_AabsI            = RabsI; 
dxdt_Afeces           = Rfeces;  
dxdt_ASI              = RSI;
dxdt_Abile            = Rbile;


// OTC Concentration in Liver 
double RL             = QL*(CPlas-CVL)*Free-Kbile*CVL+RabsI;
//-- The amount and AUC of OTC in Liver
dxdt_AL               = RL;

// OTC Concentration in Muscle
double RM             = QM*(CPlas- CVM)*Free;
//-- The amount and AUC of OTC in Muscle
dxdt_AM               = RM;

// OTC Concentration in Kidney
double RK             = QK*(CPlas - CVK)*Free - Kurine*CVK;
double Rurine         = Kurine*CVK;
//-- The amount and AUC of OTC in Kidney
dxdt_AK               = RK;
dxdt_Aurine           = Rurine;

// OTC Concentration in Fat
double RF             = QF*(CPlas - CVF)*Free;
//-- The amount and AUC of OTC in Fat
dxdt_AF               = RF;

// OTC Concentration in the rest of body
double Rrest          = QRest*(CPlas- CVRest)*Free;  
//-- The amount and AUC of OTC in the rest of body
dxdt_ARest            = Rrest;

// // {AUC equation}

dxdt_AUCCP            = CPlas;
dxdt_AUCCL            = CL;
dxdt_AUCCK            = CK;
dxdt_AUCCM            = CM;
dxdt_AUCCF            = CF;

// #+ Virtural compartment; input dose for mass balance
dxdt_ADOSE            = 0;

// Mass balance equations
double Qbal           = QC - QL - QK - QM - QRest-QF;
double Vbal           = BW - Vblood - VL - VM - VRest - VK - VF;
double Tmass          = AF + AL + AM + ARest + AK ;
double Loss           = Aurine+Abile;
double ATotal         = Tmass + Loss;

$TABLE
capture Plasma        = CPlas;
capture Liver         = CL;
capture Kidney        = CK;
capture Muscle        = CM;
capture Fat           = CF;
capture AUC_CP        = AUCCP;
capture AUC_CL        = AUCCL;
capture AUC_CK        = AUCCK;
capture AUC_CM        = AUCCM;
capture AUC_CF        = AUCCF;

'
