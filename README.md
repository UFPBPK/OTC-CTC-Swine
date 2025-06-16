# OTC-CTC-Swine

This repository contains the average and population PBPK model code files used in our oxytetracycline (OTC) and chlortetracycline (CTC) PBPK study in swine.

The manuscript title of this study is: “Withdrawal interval estimations of oxytetracycline and chlortetracycline in swine using physiologically based pharmacokinetic models: a global trade perspective”. This manuscript is under review.

The 'AASV_OTC&CTC_PBPK' folder is organized into two subfolders, one for each drug:

'OTC' – containing files for oxytetracycline. 'CTC' – containing files for chlortetracycline.

OTC Folder

The 'OTC' folder contains: (1) General PBPK model code: 'Swine_OTC_pbpk.R' (2) Calibration and evaluation function: 'Modfit_OTC_Swine.R' (3) Population PBPK model for withdrawal interval calculations: 'Population_PBPK_OTC.R' (4) Experimental pharmacokinetic raw data: 'Data_OTC.csv' (used for model calibration and evaluation)

CTC Folder

The 'CTC' folder contains: (1) General PBPK model code: 'Swine_CTC_pbpk.R' (2) Calibration and evaluation functions: 'Modfit_OTC_Swine_PO.R' (for oral administration via feed) 'Modfit_OTC_Swine_POW.R' (for oral administration via water) (3) Population PBPK models: 'Population_PBPK_CTC_PO.R' (for feed-based administration) 'Population_PBPK_CTC_POW.R' (for water-based administration) (4) Experimental pharmacokinetic raw data: 'Data_CTC.csv' (used for model calibration and evaluation)
