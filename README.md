# Bioinfo_Project
1. given csv files (from Vital):
    1.1 "diffacto_weightedsum.csv" is the protein intensity determined by OpenMS, committed by Vital
    1.2 "goterm.csv" contains protein annotation
    1.3 "category_loc.csv" contains additional localization of protein: whether it is cytoplasmic or mb-embedded
    
 2. created csv file (by me)
    2.1 "pvalue_weightedsum.csv" 
    2.2 "qvalue_weightedsum.csv"
    2.3 "median.csv"
    2.4 "localization.csv"
    2.5 "analysis.csv" (by PANTHER) from GO Enrichment analysis
    
 3. python file
    3.1 "f_oneway.py" is the script to calculate distribution, p-value and if the p-value is significant, based on one-way ANOVA. It takes:
        - input: "diffacto_weightedsum.csv"
        - output: "pvalue_weightedsum.csv"
    3.2 "qvalue.py" corrects the p-values from "f_oneway.py" based on multiple hypothesis
        - input: "diffacto_weigtedsum.csv"
        - output: "qvalue_weightedsum.csv"
    3.3 "median.py" calculates the median of 5 intensity levels of X60, X100, X200, X300, X1000
        - input: "qvalue_weightedsum.csv"
        - output: "median.csv"
    3.4 "localization.py" is the extension of "median.py". It first calculates median of 5 intensity, then matches the protein name with "category_loc.csv" and list the localization of each protein on the last column.
        - input1: "qvalue_weightedsum.csv"
        - input2: "category_loc.csv"
        - output: "localization.py"
    3.5 "VoverA.py" calculates the volume over surface area in 5 intensity levels and plot a xy graph to represent the change.
        - input: "localizatoin.csv"
        - output: a figure
        - Volume estimation: intensity sum of (median of "Cytoplasmic" proteins)
        - Area estimation: intensity sum of (median of "CytoplasmicMembrane" + "Periplasmic" + "OuterMembrane" proteins
    3.6 "Vocalno.py" plots the log2(FDR) over fold change of GO analysis results
        - input: "analysis.csv"
        - output: a figure
        
 4. Facts
    4.1 Alltogether there are 1394 proteins
    4.2 ANOVA p < 0.05 = 802 proteins
    4.3 ANOVA q < 0.05 = 724 proteins
    4.4 Altogether there are five kinds of localization by running np.unique = ["Cytoplasmic", "CytoplasmicMembrane", "OuterMembrane", "Periplasmic", "Extracellular", "Unknown"]
    4.5 V to A of 5 intensity = VtoA = [ 4.11469476,  3.93226573,  3.65069334,  3.90560723,  3.78239169]
'
    
    
     
    
