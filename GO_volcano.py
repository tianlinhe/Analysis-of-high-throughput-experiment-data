import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

#a function to log base 2
def f(nummer):
    return((-1)*math.log(nummer,2))

def volcano(input1, input2):

    df = pd.read_csv(input1)
    df=df.sort_values(['fold_change'], ascending = False)

    foldchange = np.array(df.loc[:,'fold_change'])
    penrichment = np.array(df.loc[:,'p_value_enrichment'])
    pdepletion = np.array(df.loc[:,'p_value_depletion'])
    #print (fc1)

    plt.figure(figsize=(6,6), dpi=100) #set figure size, it must precedes any plt term

    f1 = np.vectorize(f, otypes=[np.float])
    enrichment = 0
    for item in foldchange:
        if item>1:
            enrichment +=1
#to split the lists of p-enrichment and p-depletion according to fold change at 1
    fc = np.split(foldchange,[enrichment,len(foldchange)])
    p_enrichment=np.split(penrichment,[enrichment,len(penrichment)])
    p_depletion=np.split(pdepletion,[enrichment,len(pdepletion)])
    p_value = np.append (p_enrichment[0],p_depletion[1])

    Over = plt.scatter(-f1(fc[0]),f1(p_enrichment[0]),color='gold')
    Under = plt.scatter(-f1(fc[1]),f1(p_depletion[1]),color='cornflowerblue')
    df_go = pd.read_csv(input2)
    func_enrichment = []
    for i in range (len(fc[0])):
        if p_enrichment[0][i]<0.05:
            #plt.annotate(str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process'][-13:]),(-f1(fc2[0][i]),f1(fdr2[0][i])))
            for index, row in df_go.iterrows():
                if df.iloc[i,0]== row[4]:
                    print ('enrichment:' + str(df.iloc[i,0]) + '  '+str(p_enrichment[0][i]) + '  '+str(row[3]) )
                    break



    for i in range (len(fc[1])):
        if p_depletion[1][i]<0.05:
            #plt.annotate(str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process'][-13:]),(-f1(fc2[1][i]),f1(fdr2[1][i])))
            for index, row in df_go.iterrows():
                if df.iloc[i,0]== row[4]:
                    print ('depletion:' + str(df.iloc[i,0]) + '  '+str(p_depletion[1][i]) + '  ' +str(row[3]) )
                    break


    plt.legend((Over,Under),('GO Enrichment','GO Depletion'),loc='upper center')
    #y=4.321928, it correspond to q-value = log2(0.05)

    plt.axhline(y=4.321928,color='dimgrey')
    plt.text(-0.8,4.5,'p-value = 0.05', fontsize=12, color='black')

    plt.xlabel('Log2 (Fold Change)')
    plt.ylabel('- Log2 (p-value)')
    plt.title("GO Enrichment Analysis")

    plt.show()


volcano('enriched_GO.csv','goterm.csv')

#enrichment:GO:0015935  0.0089680440864  small ribosomal subunit
#enrichment:GO:0006950  0.0292899641767  response to stress
#enrichment:GO:0051287  0.0115624845586  NAD or NADH binding
#enrichment:GO:0006457  0.0200395283129  protein folding
#enrichment:GO:0030089  0.0304576206735  phycobilisome
#enrichment:GO:0030170  0.0397454882266  pyridoxal phosphate binding
#enrichment:GO:0015979  0.0468998351316  photosynthesis
#enrichment:GO:0005524  0.0433097909722  ATP binding
#depletion:GO:0008443  0.0124761022915  phosphofructokinase activity
#depletion:GO:0016805  0.0243068984316  dipeptidase activity
#depletion:GO:0050821  0.0128443927037  protein stabilization
