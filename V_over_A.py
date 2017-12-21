import pandas as pd
import numpy as np
def VoverA (df_loc):
    df = pd.read_csv(df_loc)
    df_cytoplasmic = df. groupby(['Localization']).get_group('Cytoplasmic')
    df_cytoMB = df. groupby(['Localization']).get_group('CytoplasmicMembrane')
    df_peri = df. groupby(['Localization']).get_group('Periplasmic')
    df_outerMB = df. groupby(['Localization']).get_group('OuterMembrane')
    #print (df)
    #print (df_cytoplasmic)
    #print (df_cytoMB)
    area = [[],[],[]]
    volume = []
    for column_num in range (2,7):
        area[0].append(df_cytoMB.iloc[:,column_num].sum())
        area[1].append(df_peri.iloc[:,column_num].sum())
        area[2].append(df_outerMB.iloc[:,column_num].sum())
        volume.append(df_cytoplasmic.iloc[:,column_num].sum())

    area_sum = np.sum(area, axis = 0)
    p = np.divide(volume,area_sum)


    print (volume)
    print (area_sum)
    print (p)

#VoverA("localization.csv")
VtoA = [ 4.11469476,  3.93226573,  3.65069334,  3.90560723,  3.78239169]

import matplotlib.pyplot as plt
x = [1,2,3,4,5]
x_axis = ['X60','X100','X200','X300','X1000']
VtoA = [ 4.11469476,  3.93226573,  3.65069334,  3.90560723,  3.78239169]
plt.xticks(x, x_axis)
#plt.plot([1,2,3,4,5],[ 4.11469476,  3.93226573,  3.65069334,  3.90560723,  3.78239169])
plt.plot(x,VtoA)
plt.ylim((0,5))
plt.xlabel('Light Intensity')
plt.ylabel('Volume/ Surface Area')
plt.grid(True)
plt.show()

