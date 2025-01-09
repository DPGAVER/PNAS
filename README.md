# PNAS
Ventilation Work Repository
Gaver_PNAS_Ventilation.m

**Description:**
This MATLAB script analyzes pressure-volume (PV) loops from respiratory data to calculate components of ventilation work. 

**Key Features:** 
•	Loads experimental data from specified directories. 
•	Calculates input work, total work, tissue work, and recruitment work. 
•	Generates visualizations of PV loops and recruitment regions. 
•	Saves calculated work values to Excel files. 

**Usage:** 
1. Prepare Input Files: 
•	PNAS_DataFiles.txt: A text file listing the directories containing the experimental data in the root directory. 
•	Data files are in the directory “PNAS data”
o	Ensure that each directory contains: 
o	`header.txt`: A file with experimental parameters (Tlow, Phigh). 
o	`files.xls`: An Excel file listing the filenames of respiratory data files. 
o	Respiratory data files in the specified format. 

2. Run the Script: 
•	Execute the `Gaver_PNAS_Ventilation.m` script in MATLAB. 
•	The script will process the data, generate figures, and save the results to Excel files. 

Calculations:
1.	Volume and Pressure Calculations:
o	Computes transpulmonary pressure (Ptp) and volume (V) using the given flow (Q) and airway pressure (Pairway).
o	Calculates the tissue pressure (P_tissue) from the P_airway, flow Q, and resistance, R.
2.	Work Calculations:
o	Calculates input work, total work, tissue work, and recruitment work from PV loops.

Output Files: 
•	PNAS_Ventilation_Work.xlsx: Contains work values for each file. 
•	PNAS_Ventilation_Work_Combined.xlsx: Contains work values for all files in a combined format. 
•	Figures saved in directories within each data directory:
o	`PNAS_Figures_PairwayV_Recruit`  
o	`PNAS_Figures_PV` 

Dependencies: 
•	MATLAB with necessary toolboxes 

Parameters: 
Some parameters can be adjusted within the script, such as: 
•	Filter cutoff frequency 
•	Inflection point detection criteria 

Notes:
•	The accuracy of the results depends heavily on the quality of the input data. 
•	The definition and calculation of recruitment work may vary depending on the specific research context. 
•	This script provides a basic framework. Adaptations may be necessary based on specific data characteristics and research needs. 
•	Sample datafiles provided are:
o	OD+RD+: "Pig 22 calibrated data files"
o	OD+RD-: "Pig 45 calibrated data files"
o	OD-RD+: "Pig 30 calibrated data files"
o	OD-RD-: "Pig 27 calibrated data files"

Contact:
Donald P. Gaver 
dpg@tulane.edu

