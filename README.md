**PNAS Data and Script Repository**
<ul>
  Ventilation Work Repository
</ul>

**Normalized Compliance Data**
<ul>
  Normalized Compliance.xlsx
</ul>

**PF Ratio Data**
<ul>
  PF Ratio Data.xlsx
</ul>

**MATLAB Script**
<ul>
  Gaver_PNAS_Ventilation.m
</ul>

**Description:**
This MATLAB script analyzes pressure-volume (PV) loops from respiratory data to calculate components of ventilation work. 

**Key Features:** 
<ul>
<li>Loads experimental data from specified directories. </li>
<li>Calculates input work, total work, tissue work, and recruitment work. </li>
<li>Generates visualizations of PV loops and recruitment regions. </li>
<li>Saves calculated work values to Excel files. </li>
</ul>

**Usage:** 
1. Prepare Input Files: 
<ul>
  <li>	PNAS_DataFiles.txt: A text file listing the directories containing the experimental data in the root directory. </li>
<li>	Data files are in the directory “PNAS data”</li>
<ul>
  <li>	Ensure that each directory contains: </li>
<li>	`header.txt`: A file with experimental parameters (Tlow, Phigh). </li>
<li>	`files.xls`: An Excel file listing the filenames of respiratory data files. </li>
<li>	Respiratory data files in the specified format.</li>
</ul>  
</ul>

2. Run the Script: 
<ul>
  <li>Execute the `Gaver_PNAS_Ventilation.m` script in MATLAB. </li>
<li>The script will process the data, generate figures, and save the results to Excel files. </li>
</ul>
  
**Calculations:**
1.	Volume and Pressure Calculations:
<ul>
  <li>Computes transpulmonary pressure (Ptp) and volume (V) using the given flow (Q) and airway pressure (Pairway).</li>
<li>Calculates the tissue pressure (P_tissue) from the P_airway, flow Q, and resistance, R. </li>
</ul>
2.	Work Calculations:
<ul>
  <li>Calculates input work, total work, tissue work, and recruitment work from PV loops. </li>
</ul>

**Output Files:** 
<ul>
  <li>NAS_Ventilation_Work.xlsx: Contains work values for each file.</li>
 <li>PNAS_Ventilation_Work_Combined.xlsx: Contains work values for all files in a combined format.</li>
 <li>Figures saved in directories within each data directory:</li>
<ul>
 <li>`PNAS_Figures_PairwayV_Recruit` </li>
<li>`PNAS_Figures_PV` </li>
</ul>
</ul>

**Dependencies:**
<ul>
  <li>	MATLAB with necessary toolboxes </li>
</ul>

**Parameters:**

Some parameters can be adjusted within the script, such as: 
<ul>
  <li>	Filter cutoff frequency  </li>
 <li>	Inflection point detection criteria</li>
</ul>
 

**Notes:**
<ul>
  <li>	The accuracy of the results depends heavily on the quality of the input data.  </li>
  <li>The definition and calculation of recruitment work may vary depending on the specific research context. </li>
  <li>Sample datafiles provided are:</li>
<ul>
     <li>OD+RD+: "Pig 22 calibrated data files"</li>
     <li>OD+RD-: "Pig 45 calibrated data files"</li>
     <li>OD-RD+: "Pig 30 calibrated data files"</li>  
    <li>OD-RD-: "Pig 27 calibrated data files"</li>  
</ul>
  </ul>


**Contact:**
<ul>
Donald P. Gaver (dpg@tulane.edu)
  </ul>

