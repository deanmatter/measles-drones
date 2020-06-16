import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
import pandas as pd

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
        plt.scatter(v,c,s=50,color=col,marker=m,label=l) 

df = pd.read_csv('results/strategies/strategy_comparisons.csv')
df_capped = df[df['Delivery strategy'] != 'uncapped']
print(df_capped)

#plt.xlim()
#plt.ylim()
plt.title("Comparison of team and vaccine resource allocation strategies")
plt.ylabel("")
plt.xlabel("")

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

