#volcano plot from panther
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
def f(fdr):
    return((-1)*math.log(fdr,2))
def volcano(input, output):

    df = pd.read_csv(input)
    df.drop(df.index[[0,1,2,3]], inplace=True)
    df.columns = df.iloc[0]
    df.drop(df.index[0],inplace = True)
    #fdr = np.array(df.loc[:,'Client Text Box Input (FDR)'],dtype=float)
    foldchange = np.array(df.loc[:,'Client Text Box Input (fold Enrichment)'])

#as foldchange < 0.01 means the GO does not exist, it should be excluded in the plot
    count = 0
    for item in foldchange:
        if item.startswith(' <'):
            count += 1
    #print (count)
    df= df[:count-1]

#sort the dataframe according to float(P-value)
    Pvalue = []
    for row_num in range (len(df)):
        p = float(df.iloc[row_num]['Client Text Box Input (raw P-value)'])
        Pvalue.append(p)
    df['Client Text Box Input (raw P-value)'] = Pvalue
    #print (df)
    sorted_df=df.sort_values(['Client Text Box Input (raw P-value)'], ascending = False)

    Qvalue = [1]


#calculate q-value and append it as last column
    for row_num in range (len(sorted_df)):
    #print (sorted_df.iloc[row_num]['Pvalue'])
#important, q-value is the minimum FDR, so it always compares it with the preceding term, her the last term in the Qvalue list.
        q = min(len(sorted_df)*sorted_df.iloc[row_num]['Client Text Box Input (raw P-value)']/(len(sorted_df)-row_num),Qvalue[-1])
        Qvalue.append(q)
#to remove the first term (1) from the list of q-value
    Qvalue.remove(1)
    sorted_df['Qvalue'] = Qvalue
    #print (sorted_df)


    sorted_df.to_csv(output)
    sorted_df= sorted_df.sort_values(['Client Text Box Input (fold Enrichment)'], ascending = False)
    #print (sorted_df)

    fdr1 = np.array(sorted_df.loc[:,'Qvalue'],dtype=float)
    fc1 = np.array(sorted_df.loc[:,'Client Text Box Input (fold Enrichment)'], dtype=float)
    #print (fc1)
#scatterplot
    plt.figure(figsize=(6,6), dpi=100)
    f1 = np.vectorize(f, otypes=[np.float])
    enrichment = 0
    for item in fc1:
        if item>1:
            enrichment +=1
    fc2=np.split(fc1,[enrichment,len(fc1)])
    fdr2=np.split(fdr1,[enrichment,len(fc1)])
    Over = plt.scatter(-f1(fc2[0]),f1(fdr2[0]),color='gold')

    for i in range (len(fc2[0])):
        if f1(fdr2[0][i])>4.321928:
            #plt.annotate(str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process'][-13:]),(-f1(fc2[0][i]),f1(fdr2[0][i])))
            print (str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process']),(-f1(fc2[0][i]),f1(fdr2[0][i])))


    Under = plt.scatter(-f1(fc2[1]),f1(fdr2[1]),color='cornflowerblue')
    for i in range (len(fc2[1])):
        if f1(fdr2[1][i])>4.321928:
            #plt.annotate(str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process'][-13:]),(-f1(fc2[1][i]),f1(fdr2[1][i])))
            print (str(sorted_df.loc[i,'PANTHER GO-Slim Biological Process']),(-f1(fc2[1][i]),f1(fdr2[1][i])))

    plt.legend((Over,Under),('Over-representation','Under-representation'),loc='upper right')
    #y=4.321928, it correspond to q-value = log2(0.05)

    plt.axhline(y=4.321928,color='dimgrey')
    plt.text(0.5,4.5,'q-value = 0.05', fontsize=12, color='black')

    plt.xlabel('Log2 (Fold Change)')
    plt.ylabel('- Log2 (q-value)')
    plt.title("GO Enrichment Analysis")
    #plt.figure(figsize=(4, 5), dpi=100)
    plt.show()


volcano('analysis.csv','analysis_plus.csv')

#polysaccharide metabolic process (GO:0005976) (0.13750352374993502, array(5.43097562343034))
#regulation of translation (GO:0006417) (-0.62148837674627011, array(5.3751223886968225))
#transcription, DNA-dependent (GO:0006351) (-0.71311885221183846, array(5.14908219034523))
#cellular protein modification process (GO:0006464) (-0.73696559416620622, array(5.43097562343034))
#monosaccharide metabolic process (GO:0005996) (-1.4739311883324124, array(8.508978135431613))
