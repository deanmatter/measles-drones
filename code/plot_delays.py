import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150


x_values = np.array([7,10,12,13,15,16,18,19,22])

percentages = {
    'Monocentric': np.array([0.854345439754455,	0.916405953301491,	0.95203634937124,	0.968506034407494,	1,	1.01367600407117,	1.03878468838065,	1.04931264786954,	1.07776870210947]),
    'Polycentric': np.array([0.966544115651693,	0.979421349089093,	0.989001096182039,	0.991732225864327,	1,	1.00381685642894,	1.01074828482433,	1.01334449366972,	1.02280917240084]),
    'City': np.array([0.944933384320279,	0.968082263471711,	0.982015840708799,	0.988272543589161,	1,	1.00522349433512,	1.0147662850173,	1.01904372530567,	1.02910014333318]),
    'Rural': np.array([0.761443867295545,	0.862791817004226,	0.915356868593013,	0.943001649056541,	1,	1.02542630676965,	1.0759008206153,	1.09992874775588,	1.17036508168657])
}

reg_coeffs = {      #"network": "(intercept, slope)"
    'Monocentric': (0.767424053773312,0.014878934358377),
    'Polycentric': (0.94235342605484, 0.003759369391041),
    'City': (0.912042535613728, 0.005629203481498),
    'Rural': (0.587199979550772, 0.027041025324392),
}

colors = {
    'Monocentric': 'C0',
    'Polycentric': 'C1',
    'City': 'C2',
    'Rural': 'C3',
}

plt.title("Effect of intervention delay on cases, for various network structures")
plt.xlabel("Days between epidemic detection and intervention")
plt.ylabel("Percentage change in number of cases (%)")
for key in percentages:
    plt.plot(x_values,(percentages.get(key)-1)*100,marker='.',label=key+" network",color=colors.get(key))

    reg_intercept = reg_coeffs.get(key)[0]
    reg_slope = reg_coeffs.get(key)[1]
    plt.plot(x_values,np.array([((reg_intercept + reg_slope*x)-1)*100 for x in x_values]),linestyle='dotted',color=colors.get(key))

plt.legend()

plt.savefig(f"results/sensitivity/int_delay_plot.pdf",bbox_inches='tight')
plt.show()