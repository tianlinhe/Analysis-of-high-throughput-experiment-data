import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
def f(fdr):
    return((-1)*math.log(fdr,2))

def volcano(input):

    df = pd.read_csv(input)
    df.drop(df.index[[0,1,2,3]], inplace=True)
    df.columns = df.iloc[0]
    df.drop(df.index[0],inplace = True)
    fdr = np.array(df.loc[:,'Client Text Box Input (FDR)'],dtype=float)
    foldchange = np.array(df.loc[:,'Client Text Box Input (fold Enrichment)'])


#fjdjsfkjdkjfkjfkjsdk



    f1 = np.vectorize(f, otypes=[np.float])
    #log2fdr = f(fdr)
    print (f1(fdr))


    count = 0
    for item in foldchange:
        if item.startswith(' <'):
            count += 1

    #print (count)
    small = np.random.uniform(low=0, high = 0.01, size=count)
    #print (small)
    big = np.array(foldchange[0:len(foldchange)-count],dtype=float)
    corrected_foldchange = np.concatenate([big,small])

    print (corrected_foldchange)
    plt.scatter(corrected_foldchange,f1(fdr))
    plt.xlabel('Fold Change')
    plt.ylabel('Log2 (FDR)')
    plt.title("GO Enrichment Analysis")
    plt.show()

volcano("analysis.csv")
