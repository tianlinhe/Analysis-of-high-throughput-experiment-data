# Bioinfo_Project
1. given csv files (from Vital):
    1.1 "diffacto_weightedsum.csv" is the protein intensity determined by OpenMS, committed by Vital
    1.2 "goterm.csv" contains protein annotation
    1.4 "category_loc.csv" contains additional localization of protein: whether it is cytoplasmic or mb-embedded
    
 2. created csv file (by me)
    2.1 "pvalue_weightedsum.csv" 
    2.2 "qvalue_weightedsum.csv"
   
    
 2. python file
    2.1 "f_oneway.py" is the script to calculate distribution, p-value and if the p-value is significant. It takes:
        - input: "diffacto_weightedsum.csv"
        - output: "pvalue_weightedsum.csv"
    2.2 "qvalue.py" corrects the p-values from "f_oneway.py" based on multiple hypothesis
        - input: "diffacto_weigtedsum.csv"
        - output: "qvalue_weightedsum.csv"
    
