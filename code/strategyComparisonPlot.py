import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
import pandas as pd
from collections import OrderedDict

def scatter_markers(vs, cs, s, col, ms, l):
    '''Arguments:
        vs: list of vaccination floats
        cs: list of cases floats
        s: size of scatter point
        col: colour of scatter point
        ms: list of scatter markers
        ls: list of scatter labels
    '''    
    for i in range(len(vs)):
        v = vs.iloc[i]
        c = cs.iloc[i]
        m = ms[i]
        print(v,c,s,col,m,l)
        scatters.append(ax.scatter(v,c,s=50,color=col,marker=m,label=l))

# Read in input from CSV
df = pd.read_csv('results/strategies/strategy_comparisons.csv')
df_capped = df[df['Delivery strategy'] != 'uncapped']
print(df_capped)

# Prepare the figure for scatters
fig, ax = plt.subplots()
scatters = []
ax.set_title("Comparison of team and vaccine resource allocation strategies")
ax.set_ylabel("")
ax.set_xlabel("")

# Dicts for information
ts_names = {
    'I':'Infections',
    'S':'Susceptible',
    'N':'Population',
    'EPE':'EPE Team Method',
    'I/N':'Infection Ratio',
    'spread':'Spread'
}
team_colours = {
    'I':'gold',
    'S':'crimson',
    'N':'k',
    'EPE':'forestgreen',
    'I/N':'darkorange',
    'spread':'cornflowerblue'
}
vacc_markers = {
    'I':'$I$',
    'S':'$S$',
    'N':'$N$',
    'EPE':'$E$',
    'absI':'$i$',
    'absS':'$s$',
    'absN':'$n$',
}

# Iterate through input to scatter points
for ts in ['I','S','N','EPE','I/N','spread']:
    ts_df = df_capped[df_capped['Team strategy']==ts]
    label_keys = ts_df['Delivery strategy'].values
    label_values = [vacc_markers[key] for key in label_keys]
    scatter_markers(
        vs=ts_df['Vaccinations'], 
        cs=ts_df['Cases'], 
        s=50,
        col=team_colours[ts], 
        ms=label_values,
        l=ts_names[ts]
    )
    

plt.legend()
plt.show()
quit()



handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels,handles))
leg1 = ax.legend(by_label.values(), by_label.keys(), title="Delivery strategy (colour):")
plt.gca().add_artist(leg1)

handles, labels = scatters, teamStrat3
by_label = OrderedDict(zip(labels,handles))
leg2 =ax.legend(by_label.values(), by_label.keys(), loc ='lower left', title="Team strategy (letter):")

plt.legend()
plt.show()

