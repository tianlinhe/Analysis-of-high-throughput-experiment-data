import pandas as pd
from scipy import stats
#import the csv file
df = pd.read_csv('diffacto_weightedsum.csv')

#to format it in a more comprehensible way, remove the first column and set the protein name as row name, the new csv named "my_df"
df.drop(df.columns[[0]], axis=1, inplace=True)
my_df= df.set_index('protein')

F_distribution = []
Pvalue = []
Result = []
for row_num in range (len(my_df.index)):

    f, p = stats.f_oneway(my_df.iloc[row_num][0:4],my_df.iloc[row_num][4:8],my_df.iloc[row_num][8:12],my_df.iloc[row_num][12:16],my_df.iloc[row_num][16:])
    F_distribution.append(f)
    Pvalue.append(p)
#here we set the p-value threshold as 0.05
    if p < 0.05:
        Result.append('yes')
    else:
        Result.append('no')

#add F_distribution, Pvalue and results as columns at the end of the table

my_df['F_distribution'] = F_distribution
my_df['Pvalue'] = Pvalue
my_df['Result'] = Result
#print (my_df)
#export the results as a new csv file
my_df.to_csv('pvalue_weightedsum.csv')
