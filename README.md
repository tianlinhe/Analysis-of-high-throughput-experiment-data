# Bioinfo_Project
1. csv files:
    1.1 "diffacto_weightedsum.csv" is the protein intensity determined by OpenMS, committed by Vital
    1.2 "goterm.csv" contains protein annotation
    1.3 "p_value_weightedsum.csv" description goes to 2.1
    1.4 "category_loc.csv" contains additional localization of protein: whether it is cytoplasmic or mb-embedded
    
 2. python file
    2.1 "f_oneway.py" is the script to calculate distribution, p-value and if the p-value is significant. It takes:
        - input: "diffacto_weightedsum.csv"
        - output: "pvalue_weightedsum.csv"
    
