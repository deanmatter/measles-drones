import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

with open("results/strategies/strategy_comparisons.csv") as f:


plt.xlim()
plt.ylim()
plt.title()
plt.ylabel()
plt.xlabel()
line = plt.scatter(0, 0, s = 50, color="k", label="No vaccinations.")
line.set_clip_on(False)

#fitting trendline
linear_reg = LinearRegression()
X = vehicleCosts
X.extend(droneCosts)
X.extend(combinedCosts)
X = np.array(X).reshape(-1,1)
Y = vehicleDeaths
Y.extend(droneDeaths)
Y.extend(combinedDeaths)
Y = np.array(Y).reshape(-1,1)
linear_reg.fit(X, Y)
y_pred = linear_reg.predict(X)
#plt.plot(X, y_pred, color="k", linewidth=0.1)

plt.legend()
plt.show()

