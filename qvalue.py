import pandas as pd
from scipy import stats
#import the csv file
df = pd.read_csv('diffacto_weightedsum.csv')

#to format it in a more comprehensible way, remove the first column and set the protein name as row name, the new csv named "my_df"
df.drop(df.columns[[0]], axis=1, inplace=True)
my_df= df.set_index('protein')

F_distribution = []
Pvalue = []

#for row_num in range (0,3):
for row_num in range (len(my_df.index)):

    f, p = stats.f_oneway(my_df.iloc[row_num][0:4],my_df.iloc[row_num][4:8],my_df.iloc[row_num][8:12],my_df.iloc[row_num][12:16],my_df.iloc[row_num][16:])
    F_distribution.append(f)
    Pvalue.append(p)

#add F_distribution, Pvalue and results as columns at the end of the table
sorted_Pvalue = Pvalue.sort(reverse=True)

my_df['F_distribution'] = F_distribution
my_df['Pvalue'] = Pvalue

#sorted my_df according to P-value, so that the highest P-value is displayed on the top of the list
sorted_df=my_df.sort(['Pvalue'], ascending = False)
#print (sorted_df)
#Formula of Qvalue refer to the excel file in computerlab1
Result = []
Qvalue = [1]
for row_num in range (len(sorted_df)):
#important, q-value is the minimum FDR, so it always compares it with the preceding term, her the last term in the Qvalue list.
    q = min(len(sorted_df)*Pvalue[row_num]/(len(sorted_df)-row_num),Qvalue[-1])
    if q < 0.05:
        Result.append('yes')
    else:
        Result.append('no')
    Qvalue.append(q)
#print (len(Result))
#print (len(Qvalue))

#to remove the first term (1) from the list of q-value
Qvalue.remove(1)
sorted_df['Qvalue'] = Qvalue
sorted_df['Result'] = Result
#print (Qvalue)
#export the results as a new csv file
sorted_df.to_csv('qvalue_weightedsum.csv')
