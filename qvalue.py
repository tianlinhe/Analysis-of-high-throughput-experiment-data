import pandas as pd
from scipy import stats
#import the csv file
df = pd.read_csv('diffacto_weightedsum.csv')

#to format it in a more comprehensible way, remove the first column
df.drop(df.columns[[0]], axis=1, inplace=True)

F_distribution = []
Pvalue = []

for row_num in range (len(df.index)):

    f, p = stats.f_oneway(df.iloc[row_num][1:5],df.iloc[row_num][5:9],df.iloc[row_num][9:13],df.iloc[row_num][13:17],df.iloc[row_num][17:])
    F_distribution.append(f)
    Pvalue.append(p)

#add F_distribution, Pvalue and results as columns at the end of the table

df['F_distribution'] = F_distribution
df['Pvalue'] = Pvalue

sorted_Pvalue = Pvalue.sort(reverse=True)
print (sorted_Pvalue)

#print (my_df)
#export the results as a new csv file
df.to_csv('pvalue_weightedsum.csv')

#sort df according to P-value, so that the highest P-value is displayed on the top of the list
sorted_df=df.sort_values(['Pvalue'], ascending = False)

#Formula of Qvalue refer to the excel file in computerlab1
Result = []
Qvalue = [1]

count = 0
for row_num in range (len(sorted_df)):
    #print (sorted_df.iloc[row_num]['Pvalue'])
#important, q-value is the minimum FDR, so it always compares it with the preceding term, her the last term in the Qvalue list.
    q = min(len(sorted_df)*sorted_df.iloc[row_num]['Pvalue']/(len(sorted_df)-row_num),Qvalue[-1])
    if q < 0.05:
        Result.append('yes')
        count += 1
    else:
        Result.append('no')
    Qvalue.append(q)
#print (len(Result))
#print (len(Qvalue))
print (count)
#to remove the first term (1) from the list of q-value
Qvalue.remove(1)
sorted_df['Qvalue'] = Qvalue
sorted_df['Result'] = Result
#print (Qvalue)
#export the results as a new csv file
sorted_df.to_csv('qvalue_weightedsum.csv')

