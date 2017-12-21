import pandas as pd
import numpy as np
def localization (df_id,output):
    df = pd.read_csv(df_id)

    m = np.array([['protein'],['X60'],['X100'],['X200'],['X300'],['X1000']])

    #for row_num in range (0,3):
        #print (df_id.iloc[row_num][1])
    for row_num in range (len(df.index)):
        #print (np.array(df.iloc[row_num][0:5]))

        protein = np.array(df.iloc[row_num][1])
        m_60 = np.median(np.array(df.iloc[row_num][2:6]))
        m_100 = np.median(np.array(df.iloc[row_num][6:10]))
        m_200 = np.median(np.array(df.iloc[row_num][10:14]))
        m_300 = np.median(np.array(df.iloc[row_num][14:18]))
        m_1000 = np.median(np.array(df.iloc[row_num][18:22]))
        m_m = np.array([protein,m_60,m_100,m_200,m_300,m_1000])
        m = np.insert(m,1, m_m, axis=1)


        #print (m0)
    print (m)
    data_cal = pd.DataFrame(data=m)
    df = pd.DataFrame.transpose(data_cal)
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(0))
    print (df)

    df.to_csv(output)
    #print (df)

localization("qvalue_weightedsum.csv","median.csv")
