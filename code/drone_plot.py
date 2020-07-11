import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
mpl.rcParams['figure.dpi'] = 150

# Read in the file with all drone sensitivity details
filename = "results/sensitivity/drone_investigation.csv"
df = pd.read_csv(filename)
print(df.head())

# Set plot characteristics
plt.xlabel('Number of drones')
plt.ylabel('Change in number of cases per epidemic (%)')
plt.title('Impact of number of drones on cases, for various population levels')

# Plot a curve for each population level
for pop_multiplier in df['PopMultiplier'].unique():
    df_filtered = df[df['PopMultiplier']==pop_multiplier]
    base_cases = df_filtered[df_filtered['drones']==5]['Cases'].values[0]
    
    df_scaled = df_filtered
    df_scaled['Cases'] = df_filtered['Cases'] / base_cases * 100 - 100

    plt.plot(
        df_scaled['drones'],
        df_scaled['Cases'],
        label=str(int(pop_multiplier*100)) + '%'
        )

plt.legend(title="Population level")
plt.savefig(f"results/sensitivity/absN vaccStrat/numDrones.pdf",bbox_inches='tight')
plt.show()
plt.close()