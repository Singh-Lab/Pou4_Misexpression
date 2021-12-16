# README
#### LM_inputs directory

This directory contains all the code and lists pertaining to Linear Model analyses.
It is self-contained, in that you can run the code within this directory. 

## **Code**
(1) "final_nmf_pipeline.m" is the function used in:
(2) "final_run_nmf_pipeline.m"

## **Input to "final_run_nmf_pipeline.m"**
(1) "oe*_data"
- "oe*_data" consists of expression vectors of OE cells within a given OE Subcluster.  

(2) "t_av_all_wt" csv file
- "t_av_all_wt" contains average expression vectors of all Wildtype Cell Types (4 Sensory Cell Types and WT Epidermis). 

### **Notes**
- These files  were generated from "final_pou4_analysis.pipeline.Rmd".
-  **Change the value of the "oe_csv" variable in "final_run_nmf_pipeline.m" to generate various outputs.**
- Outputs include:
	- Histogram of all Cell-Type Specific Solved Coefficients ("oe*_LM.png").
	- Stacked barplots of different combinations of BTN and other Sensory Cell Type Solved Coefficients ("oe*_BTN/*.png").
	-  Statistics of average Cell-type Specific Coefficients of each Sensory Cell Type and the differences between them and average BTN-specific Solved Coefiicent ("_oe*_stats.txt").

## **Outputs of "final_run_nmf_pipeline.m"**
- Present in **"subclusters" subdirectory**.
- **"subclusters" subdirectory** consists of various **"oe directories"**, which correspond to results from the respective oe subcluster. 

## **Example Usage**
i.e. Running **"final_run_nmf_pipeline.m"** with **"t_av_all_wt.csv"** and **"oe6_data.csv"** generated outputs that can be found in **"subclusters/oe6"** directory.

### Note
Results are for solved coefficients with respect to 4 Sensory Cell Types. WT Epidermis is excluded from these analyses.

