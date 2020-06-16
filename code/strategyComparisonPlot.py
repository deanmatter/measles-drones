import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
import pandas as pd

df = pd.read_csv('results/strategies/strategy_comparisons.csv')
print(df)

#plt.xlim()
#plt.ylim()
plt.title("Comparison of team and vaccine resource allocation strategies")
plt.ylabel("")
plt.xlabel("")

team_colours = {
    'I':'',
    'S':'',
    'N':'',
    'EPE':'',
    'I/N':'',
    'spread':'',
}
vacc_letters = {
    'I':'I',
    'S':'S',
    'N':'N',
    'EPE':'E',
    'absI':'i',
    'absS':'s',
    'absN':'n',
}


for ts in ['I','S','N','EPE','I/N','spread']:
    ts_df = df[df['Team strategy']==ts]
    plt.scatter(ts_df['Vaccinations'], ts_df['Cases'], s = 50, label=ts_df['Delivery strategy'])

plt.legend()
plt.show()

