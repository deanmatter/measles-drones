import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

#===============================================================================
# Input
#===============================================================================
givens = np.array([23,
22,
25,
25,
40,
60,
85,
120,
123,
159,
150,
135,
17,
8,
19,
0
])

Iavgs = np.array([25.7	,
27.9	,
34.7	,
43.7	,
54.9	,
68.3	,
82.8	,
99.2	,
112.8	,
131.4	,
130.8	,
124.3	,
114.4	,
102.4	,
90.3	,
77.2	

])
x = np.array([44	,
45	,
46	,
47	,
48	,
49	,
50	,
51	,
52	,
53	,
54	,
55	,
56	,
57	,
58	,
59	
])

histIn = []
for i in range(0,len(x)):
	for numTimes in range(0,givens[i]):
		histIn.append(x[i])
print(histIn)



title = "Validation of simulation using Matadi 2005-6 outbreak data."
xlabel = "Week of year 2005-2006."

xt = [44,46,48,50,52,54,56,58]
xtLabels = [44,46,48,50,52,2,4,6]

#===============================================================================
# Plotting
#===============================================================================
fig, ax1 = plt.subplots()
plt.title(title)
color='k'
ax1.set_ylim(0,200)
ax1.set_xlabel(xlabel)
ax1.set_ylabel('Actual reported measles cases per week', color=color)
ax1.hist(histIn, bins=x,rwidth=0.9,color='k')
#ax1.scatter(x[int((len(x))/2)], costs[int((len(x))/2)],color=color,marker='s',s=50)
#ax1.scatter(x[-1], costs[-1],color=color,marker='s',s=50)
ax1.tick_params(axis='y', labelcolor=color)
#x_enumerated = range(min(x), int(np.ceil(max(x)))+1)
plt.xticks(xt, xtLabels)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylim(0,200)
color = 'darkred'
ax2.set_ylabel("Simulated average infections per week", color=color)  # we already handled the x-label with ax1
ax2.plot(x+0.5, Iavgs, color=color, marker='o')
#ax2.scatter(x[int((len(x))/2)], deaths[int((len(x))/2)],color=color,marker='s',s=50)
#ax2.scatter(x[-1], deaths[-1],color=color,marker='s',s=50)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

fig.savefig("Matadi_updated.pdf",bbox_inches='tight')
plt.show()