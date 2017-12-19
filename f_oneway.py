import pandas as pd
from scipy import stats
#import the csv file
df = pd.read_csv('diffacto_weightedsum.csv')
#print (df)
#to format it in a more comprehensible way, remove the first column
df.drop(df.columns[[0]], axis=1, inplace=True)

print (df.index)
F_distribution = []
Pvalue = []
Result = []
count = 0
for row_num in range (len(df.index)):

    f, p = stats.f_oneway(df.iloc[row_num][1:5],df.iloc[row_num][5:9],df.iloc[row_num][9:13],df.iloc[row_num][13:17],df.iloc[row_num][17:])
    F_distribution.append(f)
    Pvalue.append(p)
#here we set the p-value threshold as 0.05
    if p < 0.05:
        Result.append('yes')
        count += 1
    else:
        Result.append('no')

#add F_distribution, Pvalue and results as columns at the end of the table

df['F_distribution'] = F_distribution
df['Pvalue'] = Pvalue
df['Result'] = Result
print (count)

#export the results as a new csv file
df.to_csv('pvalue_weightedsum.csv')



