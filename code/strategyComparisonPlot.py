import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
import pandas as pd
from collections import OrderedDict

STANDARDIZE_COLUMNS = False
SAVE_FIGURE = True
SHOW_ONCE = False
vacc_type = 'targeted'

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
df_mono = pd.read_csv(f'results/strategies/strategy_comparisons_{vacc_type}_monocentric.csv')
df_poly = pd.read_csv(f'results/strategies/strategy_comparisons_{vacc_type}_polycentric.csv')
df_city = pd.read_csv(f'results/strategies/strategy_comparisons_{vacc_type}_city.csv')
df_rural = pd.read_csv(f'results/strategies/strategy_comparisons_{vacc_type}_rural.csv')

# Remove uncapped delivery strategy from each df
df_mono = df_mono[df_mono['Delivery strategy'] != 'uncapped']
df_poly = df_poly[df_poly['Delivery strategy'] != 'uncapped']
df_city = df_city[df_city['Delivery strategy'] != 'uncapped']
df_rural = df_rural[df_rural['Delivery strategy'] != 'uncapped']

if STANDARDIZE_COLUMNS:
    # Change the columns to percentage change from average
    for df in [df_mono, df_poly, df_city, df_rural]:
        # Standardize each column
        for col in ['Cases','Deaths','Vaccinations','Drone Deliveries']:
            col_mean = df[col].mean()
            df[col] = (df[col] - col_mean) / col_mean * 100
    print("Standardized columns.")
else:
    print("Using unstandardized columns.")

# Combine the three into a single network
df_capped = pd.concat([df_mono, df_poly, df_city, df_rural])
print(df_capped.head())

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
def plot_network(network_type):
    # Prepare the figure for scatters
    global ax, scatters, sc_letters
    fig, ax = plt.subplots()
    scatters = []
    sc_letters = []
    ax.set_ylabel("Total number of cases")
    ax.set_xlabel("Total number of vaccines given")
    df_filtered = df_capped[df_capped['Network type'] == network_type]

    for ts in ['I','S','N','EPE','I/N','spread']:
        ts_df = df_filtered[df_filtered['Team strategy']==ts]           # list of team alc strat shortnames
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

    ax.set_xlim(ax.get_xlim()[0],ax.get_xlim()[1]*1.23)
    ax.set_title(f"Comparison of allocation strategy pairs - {network_type}, {vacc_type}")
        
    # Delivery strategy legend, made using the labels added when calling ax.scatter()
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels,handles))
    leg1 = ax.legend(by_label.values(), by_label.keys(), loc ='upper right', title="Team strategy (colour):")
    plt.gca().add_artist(leg1)

    # Team allocation strategy legend, made manually by adding scatter points and scatter markers to scatters,sc_letters
    handles, labels = scatters, sc_letters
    by_label = OrderedDict(zip(labels,handles))
    leg2 = ax.legend(by_label.values(), by_label.keys(), loc ='lower right', title="Delivery strategy (letter):")

    if SAVE_FIGURE:
        fig.savefig(f"results/strategies/team_vaccine_scatter_{vacc_type}_{network_type}.pdf",bbox_inches='tight')
    
    plt.show()
    plt.close()
    
    if SHOW_ONCE:
        quit()

for network_type in ['monocentric','polycentric','city','rural']:
    plot_network(network_type)