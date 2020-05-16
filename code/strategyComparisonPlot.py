import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

vehicleDeaths = [2214.068,
2192.622,
3919.11,
4517.784,
2212.224,
2209.59,
2216.702,
3915.034,
4534.06,
2215.954,
2200.542,
2205.292,
3933.934,
4512.684,
2184.056,
2224.052,
2199.668,
3930.096,
4515.772,
2192.944
]
vehicleCosts = [313014.9771,
374496.8655,
181579.4847,
128945.2233,
362522.6514,
315931.4308,
375791.0601,
181325.5517,
128489.2997,
362801.0587,
312711.468,
375372.3096,
181577.0834,
128953.0488,
361751.0765,
317410.0787,
375490.0878,
182129.24,
132263.7822,
362289.5282
]
droneDeaths = [2311.222,
2202.99,
3894.6,
4554.916,
2204.144,
2320.602,
2241.928,
3978.798,
4604.582,
2252.42,
2315.574,
2221.538,
3981.45,
4615.562,
2229.454,
2343.078,
2206.656,
3829.262,
4569.526,
2216.094
]
droneCosts = [325608.2457,
376257.5823,
196580.9891,
140088.956,
364478.8668,
325669.1014,
378343.3817,
189165.9792,
136884.5801,
363924.8648,
324546.7406,
377632.458,
192239.6676,
140762.4134,
361833.35,
283793.9643,
350187.9498,
177802.2967,
116032.4541,
343144.9742
]
combinedDeaths = [2214.864,
2212.57,
3727.96,
4513.37,
2219.874,
2189.976,
2205.678,
3733.558,
4510.436,
2212.84,
2207.916,
2180.942,
3696.532,
4518.576,
2218.054,
2199.36,
2210.522,
3704.068,
4499.896,
2206.79
]
combinedCosts = [333433.6021,
375683.3796,
208919.5262,
137884.355,
364369.8701,
333093.4509,
375210.4359,
209046.7259,
138408.6166,
363960.7833,
333461.1425,
374411.9176,
210234.6492,
138582.7272,
364034.7738,
323589.1402,
374845.1569,
202595.2289,
129720.1975,
361128.8083
]

plt.xlim(250000,400000)
plt.ylim(2000,2500)
plt.title("Results comparison for various delivery and team allocation strategies\n(500 simulations each)")
plt.ylabel("Average number of deaths")
plt.xlabel("Average total delivery and vaccine cost ($)")
line = plt.scatter(0, 6504.49, s = 50, color="k", label="No vaccinations.")
line.set_clip_on(False)
plt.scatter(vehicleCosts[:],vehicleDeaths[:], label = "2 vehicles only.", s=15)
plt.scatter(combinedCosts[:], combinedDeaths[:], label = "1 drone + 1 vehicle.", s=15)
plt.scatter(droneCosts[:], droneDeaths[:], label = "2 drones only.", s=15)

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

