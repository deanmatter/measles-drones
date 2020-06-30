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
        v = vs.iloc[i]                                                      # num vaccines for this point
        c = cs.iloc[i]                                                      # num cases for this point
        m = ms[i]                                                           # the marker: like '$I$'
        scatters.append(ax.scatter(v,c,s=50,color=col,marker=m,label=l))    # scatter a single point
        sc_letters.append(vacc_strat_names[m])                              # append vacc strat name

# Read in input from CSV
df = pd.read_csv('results/strategies/strategy_comparisons_targeted_monocentric.csv')
df_capped = df[df['Delivery strategy'] != 'uncapped']

# Prepare the figure for scatters
fig, ax = plt.subplots()
scatters = []
sc_letters = []
ax.set_title("Comparison of team and vaccine resource allocation strategies")
ax.set_ylabel("Total number of cases")
ax.set_xlabel("Total number of vaccines given")
ax.set_xlim(300000,650000)

# Dicts for information
ts_names = {                        # Key: team allocation strategy shortname
    'I':'Infections',               # Value: team allocation strategy name
    'S':'Susceptible',
    'N':'Population',
    'EPE':'EPE Team Allocations',
    'I/N':'Infection Ratio',
    'spread':'Spread'
}
team_colours = {                    # Key: team allocation strategy shortname
    'I':'gold',                     # Value: scatter point colour
    'S':'crimson',
    'N':'k',
    'EPE':'forestgreen',
    'I/N':'darkorange',
    'spread':'cornflowerblue'
}
vacc_markers = {                    # Key: vaccine delivery strategy shortname
    'I':'x',                      # Value: scatter point marker
    'S':'+',  
    'N':'$P$',
    'EPE':'$E$',
    'absI':'$I$',
    'absS':'$S$',
    'absN':'$N$',
}
vacc_strat_names = {                        # Key: scatter point marker
    'x':"Infections / minute",       # Value: vaccine delivery strategy name
    '+':"Susceptible / minute",
    '$P$':"Population / minute",
    '$E$':"EPE Vaccine Delivery",
    '$I$':"Infections",
    '$S$':"Susceptible",
    '$N$':"Population"
}

# Iterate through input to scatter points
for ts in ['I','S','N','EPE','I/N','spread']:
    ts_df = df_capped[df_capped['Team strategy']==ts]           # list of team alc strat shortnames
    label_keys = ts_df['Delivery strategy'].values              # list of vacc del strat shortnames
    label_values = [vacc_markers[key] for key in label_keys]    # list of scatter point markers
    scatter_markers(                                # Call the custom scatter_markers method
        vs=ts_df['Vaccinations'], 
        cs=ts_df['Cases'], 
        s=50,
        col=team_colours[ts], 
        ms=label_values,
        l=ts_names[ts]
    )
    
# Delivery strategy legend, made using the labels added when calling ax.scatter()
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels,handles))
leg1 = ax.legend(by_label.values(), by_label.keys(), title="Delivery strategy (colour):")
plt.gca().add_artist(leg1)

# Team allocation strategy legend, made manually by adding scatter points and scatter markers to scatters,sc_letters
handles, labels = scatters, sc_letters
by_label = OrderedDict(zip(labels,handles))
leg2 =ax.legend(by_label.values(), by_label.keys(), loc ='lower right', title="Team strategy (letter):")

fig.savefig("results/strategies/team_vaccine_scatter_targeted_mono.pdf",bbox_inches='tight')
plt.show()