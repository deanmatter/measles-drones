import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
mpl.rcParams['figure.dpi'] = 150

# Read in the file with all drone sensitivity details
filename = "results/sensitivity/team_investigation.csv"
df = pd.read_csv(filename)
print(df.head())

# Set plot characteristics
plt.xlabel('Number of teams')
plt.ylabel('Change in number of cases per epidemic (%)')
plt.title('Impact of number of teams on cases')

# Plot a curve for each population level
for vacc_targeting in df['VaccType'].unique():
    df_filtered = df[df['VaccType']==vacc_targeting]
    base_cases = df_filtered[df_filtered['Teams']==15]['Cases'].values[0]
    
    df_scaled = df_filtered
    df_scaled['Cases'] = df_filtered['Cases'] / base_cases * 100 - 100

    plt.plot(
        df_scaled['Teams'],
        df_scaled['Cases'],
        label=vacc_targeting,
        marker='.'
        )

plt.legend(title="Vaccination type")
plt.savefig(f"results/sensitivity/absN vaccStrat/numTeams_cases.pdf",bbox_inches='tight')
plt.show()
plt.close()