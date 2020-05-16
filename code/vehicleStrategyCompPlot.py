import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
from collections import OrderedDict

vehicleDeaths = [2214.068    ,
2192.622    ,
3919.11    ,
4517.784    ,
2212.224    ,
2216.636    ,
2209.59    ,
2216.702    ,
3915.034    ,
4534.06    ,
2215.954    ,
2203.738    ,
2200.542    ,
2205.292    ,
3933.934    ,
4512.684    ,
2184.056    ,
2237    ,
2224.052    ,
2199.668    ,
3930.096    ,
4515.772    ,
2192.944    ,
2227.086    

]
vehicleCosts = [313014.9771    ,
374496.8655    ,
181579.4847    ,
128945.2233    ,
362522.6514    ,
280978.5904    ,
315931.4308    ,
375791.0601    ,
181325.5517    ,
128489.2997    ,
362801.0587    ,
280735.0569    ,
312711.468    ,
375372.3096    ,
181577.0834    ,
128953.0488    ,
361751.0765    ,
280785.9296    ,
317410.0787    ,
375490.0878    ,
182129.24    ,
132263.7822    ,
362289.5282    ,
281332.4821    

]
droneDeaths = [2311.222    ,
2202.99    ,
3894.6    ,
4554.916    ,
2204.144    ,
2280.22    ,
2320.602    ,
2241.928    ,
3978.798    ,
4604.582    ,
2252.42    ,
2311    ,
2315.574    ,
2221.538    ,
3981.45    ,
4615.562    ,
2229.454    ,
2326.39    ,
2343.078    ,
2206.656    ,
3829.262    ,
4569.526    ,
2216.094    ,
2320.364    

]
droneCosts = [325608.2457    ,
376257.5823    ,
196580.9891    ,
140088.956    ,
364478.8668    ,
286526.2327    ,
325669.1014    ,
378343.3817    ,
189165.9792    ,
136884.5801    ,
363924.8648    ,
285437.5875    ,
324546.7406    ,
377632.458    ,
192239.6676    ,
140762.4134    ,
361833.35    ,
287540.6439    ,
283793.9643    ,
350187.9498    ,
177802.2967    ,
116032.4541    ,
343144.9742    ,
277131.2676    

]
combinedDeaths = [2214.864    ,
2212.57    ,
3727.96    ,
4513.37    ,
2219.874    ,
2240.742    ,
2189.976    ,
2205.678    ,
3733.558    ,
4510.436    ,
2212.84    ,
2215.014    ,
2207.916    ,
2180.942    ,
3696.532    ,
4518.576    ,
2218.054    ,
2223.034    ,
2199.36    ,
2210.522    ,
3704.068    ,
4499.896    ,
2206.79    ,
2215.524    

]
combinedCosts = [333433.6021    ,
375683.3796    ,
208919.5262    ,
137884.355    ,
364369.8701    ,
281764.3967    ,
333093.4509    ,
375210.4359    ,
209046.7259    ,
138408.6166    ,
363960.7833    ,
281648.5167    ,
333461.1425    ,
374411.9176    ,
210234.6492    ,
138582.7272    ,
364034.7738    ,
282258.7427    ,
323589.1402    ,
374845.1569    ,
202595.2289    ,
129720.1975    ,
361128.8083    ,
280084.8688    

]

teamStrat = ['S','I','I/N','Spread','N','EPE','S','I','I/N','Spread','N','EPE',
             'S','I','I/N','Spread','N','EPE','S','I','I/N','Spread','N','EPE']

teamStrat3 = ['S','S','S','I','I','I','I/N','I/N','I/N','Spread','Spread','Spread','N','N','N','EPE','EPE','EPE',
             'S','S','S','I','I','I','I/N','I/N','I/N','Spread','Spread','Spread','N','N','N','EPE','EPE','EPE',
             'S','S','S','I','I','I','I/N','I/N','I/N','Spread','Spread','Spread','N','N','N','EPE','EPE','EPE',
             'S','S','S','I','I','I','I/N','I/N','I/N','Spread','Spread','Spread','N','N','N','EPE','EPE','EPE',]


fig, ax =plt.subplots()
ax.set_xlim(100000,400000)
ax.set_ylim(2000,5000)
ax.set_title("Comparison of vaccine delivery methods\n(500 simulations each)")
ax.set_ylabel("Average number of deaths")
ax.set_xlabel("Average total delivery and vaccine cost ($)")
#line = plt.scatter(0, 6504.49, s = 50, color="k", label="No vaccinations")
#line.set_clip_on(False)

scatters = []

for i in range(0, 24):
    labl= '$' + teamStrat[i] +'$'
    if teamStrat[i] == 'I/N':
        labl = "$R$"
    elif teamStrat[i] == 'Spread':
        labl = "+"
    elif teamStrat[i] == 'EPE':
        labl = "$E$"
   
    scatters.append(ax.scatter(combinedCosts[i], combinedDeaths[i], label = "1 UAV + 1 vehicle", s=25,c='steelblue', marker=labl))
    scatters.append(ax.scatter(droneCosts[i], droneDeaths[i], label = "2 UAVs", s=25, c='gold',marker=labl))
    scatters.append(ax.scatter(vehicleCosts[i],vehicleDeaths[i], label = "2 vehicles", s=25,c='mediumseagreen',marker=labl))


ax.scatter(0,0, label = "1 UAV + 1 vehicle", s=12,c='steelblue')
ax.scatter(0,0, label = "2 UAVs", s=12, c='gold')
ax.scatter(0,0, label = "2 vehicles", s=12,c='mediumseagreen')



handles, labels = plt.gca().get_legend_handles_labels()
by_label =OrderedDict(zip(labels,handles))
leg1 = ax.legend(by_label.values(), by_label.keys(), title="Delivery method (colour):")
plt.gca().add_artist(leg1)

handles, labels = scatters, teamStrat3
by_label = OrderedDict(zip(labels,handles))
leg2 =ax.legend(by_label.values(), by_label.keys(), loc ='lower left', title="Team strategy (letter):")


#plt.legend()
plt.savefig("delivMethodScatter.pdf",bbox_inches='tight')
plt.show()