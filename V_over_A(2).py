import numpy as np
import matplotlib.pyplot as plt
def to_plot(input):
    x = [1,2,3,4,5]
    x_axis = ['X60','X100','X200','X300','X1000']
#VtoA = [ 4.11469476,  3.93226573,  3.65069334,  3.90560723,  3.78239169]
    #plt.xticks(x, x_axis)


    volume =[1630018405.5000002, 1506606691.3999999, 1571388973.7999997, 1528810441.3, 1508798214.7999997]
    cytoMB = [310636253.9, 301277282.0, 331563189.2, 302710007.6, 318035204.09999996]
    periplasmic = [65291998.0, 58280215.0, 66842115.5, 57548402.0, 51930649.0]
    extracellular= [15188366.0, 11863796.0, 18789647.0, 17996252.0, 15669269.2]
    outermembrane =[20217392.0, 23582096.0, 32030458.0, 31181474.0, 28934721.0]

    fig, ax1 = plt.subplots()

    VtoA= [x / input[0] for x in input]

    #VtoA = input/float(input[0])
    ax1.plot(x, VtoA, 'gold', label='V to A', marker ='o')
    ax1.set_xlabel('Light Intensity')
# Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('Normalized Volume/ Surface Area')
    #ax1.tick_params('y', colors='b')
    ax1.legend(loc='upper left')
    plt.xticks(x, x_axis)


    ax2 = ax1.twinx()

    #line_extra = ax2.plot(x, extracellular, 'black')
    ax2.plot(x, extracellular, 'cornflowerblue', label='Extracellular',marker ='o')

    #line_peri=ax2.plot(x, periplasmic, 'r')
    ax2.plot(x, periplasmic, 'springgreen', label="Periplasmic",marker ='o')

    #line_cytoMB=ax2.plot (x,cytoMB,'yellow')
    ax2.plot (x,cytoMB,'plum', label='Cytoplasmic MB',marker ='o')
    ax2.plot (x,volume,'pink', label='Cytoplasmic (V)',marker ='o')
    ax2.plot (x,outermembrane,'lightcoral', label='Outer MB (A)',marker ='o')
    ax2.set_ylabel('Protein Intensity')

    ax2.legend(loc=1)

    ax1.set_ylim(0, 1.2)
    ax2.set_ylim(0,2200000000)
    fig.tight_layout()
    plt.show()
to_plot([ 80.62456352,  63.88773464,  49.05921026,  49.02944746,  52.14490282])
