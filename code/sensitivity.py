import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

#===============================================================================
# Input
#===============================================================================
costs = np.array([360162.6906    ,
362318.5066    ,
371654.1126    ,
404835.3583    ,
465630.4849    ,
498536.7714    ,
658834.4542    ,
816269.6129    

])
deaths = np.array([2174.776    ,
2205.338    ,
2197.034    ,
2215.448    ,
2218.01    ,
2232.202    ,
2311.344    ,
2421.418    

])
x = np.array([0    ,
1.168    ,
3.31    ,
5.758    ,
7.916    ,
8.824    ,
12.64    ,
16.428    

])
#title = "Comparison of differing levels of road quality on results."
xlabel = "Average number of impassable roads in the network."

xt = np.array([0    ,
1,3,6,9,12,15 
])


#===============================================================================
# Plotting
#===============================================================================
fig, ax1 = plt.subplots()
#plt.title(title)
color='k'
ax1.set_ylim(0,1000000)
ax1.set_xlabel(xlabel)
ax1.set_ylabel('Average total delivery and vaccine cost ($)', color=color)
ax1.plot(x, costs, color=color, marker='o')
#ax1.scatter(x[int((len(x))/2)], costs[int((len(x))/2)],color=color,marker='s',s=50)
ax1.scatter(x[1], costs[1],color=color,marker='s',s=50)
ax1.tick_params(axis='y', labelcolor=color)
#x_enumerated = range(min(x), int(np.ceil(max(x)))+1)
plt.xticks(xt)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'darkred'
ax2.set_ylim(0,6500)
ax2.set_ylabel("Average number of deaths", color=color)  # we already handled the x-label with ax1
ax2.plot(x, deaths, color=color, marker='o')
#ax2.scatter(x[int((len(x))/2)], deaths[int((len(x))/2)],color=color,marker='s',s=50)
ax2.scatter(x[1], deaths[1],color=color,marker='s',s=50)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

fig.savefig("Sensitivity/roads.pdf",bbox_inches='tight')
plt.show()