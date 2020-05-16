import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150


costs = np.array([215685.2172,
332809.2184,
348866.0459,
357534.0787,
361837.6713,
366478.6284,
368842.1877,
374289.6082,
372477.8014
])
deaths = np.array([3810.158,
2655.72,
2434.688,
2254.99,
2197.862,
2171.658,
2042.54,
2006.842,
1842.806
])
x = np.array([5,10,12,14,15,16,18,20,25])



fig, ax1 = plt.subplots()

plt.title("Comparison of differing numbers of teams on costs and deaths")
color='k'
ax1.set_ylim(0,400000)
ax1.set_xlabel('Number of vaccination teams')
ax1.set_ylabel('Average total delivery and vaccine cost ($)', color=color)
ax1.plot(x, costs, color=color, marker='o')
ax1.tick_params(axis='y', labelcolor=color)
#x_enumerated = range(min(x), int(np.ceil(max(x)))+1)
plt.xticks(x)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'darkred'
ax2.set_ylim(0,6500)
ax2.set_ylabel("Average number of deaths", color=color)  # we already handled the x-label with ax1
ax2.plot(x, deaths, color=color, marker='o')
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

#fig.savefig("test.pdf",bbox_inches='tight')
plt.show()