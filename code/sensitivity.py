import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
untargeted_base = [178297,5679,724251,8782]                  # The base results set
targeted_base = [148102,4718,462797,5798]

param_names = {
    "exposedDays": "Days exposed",
    "infectiousDays": "Days infectious",
    "R0": "$R_{0}$",
    "migrationIntensity": "Amount of migration",
    "deathRate": "Mortality rate",
    "vaccineEffectiveness": "Vaccine efficacy - susceptible",
    "prophylaxis72hrSuccessRate": "Vaccine efficacy - exposed",
    "numTeams": "Number of response teams",
    "workingMinutesPerDay": "Working hours per day",
    "maxVaccsTeamDay": "Daily vaccinations per team",
    "droneVaccineCapacity": "Vaccines per drone delivery",
    "monoDaysPotency": "Days of vaccine potency",
    "interventionCaseRatio": "Epidemic detection delay",
    "interventionLeadTime": "Intervention delay in days",
    "numberOfDrones": "Number of drones available",
    "droneSpeed": "Average drone flight speed",
    "flightLaunchTime": "Launch time per drone flight",
    "vaccinationRate": "Population vaccination rate",
    "maxDistance": "Diameter of network",
    "populationMultiplier": "Total population in network",
    "numTeams targeted": "Number of response teams, targeted",
    "numTeams untargeted": "Number of response teams, untargeted"
}

def plot_param_set(p_set):
    xs = [-50,-30,-20,-10,10,20,30,50]
    xlabel = f"Percentage change in parameter"
    metrics = ["cases","deaths","vaccines","drone deliveries"]

    for i, metric in enumerate(metrics):
        if metric in ["deaths", "drone deliveries"]:
            # These metrics are not plotted for sensitivity by default
            continue

        title = f"Sensitivity of {metric} to changes in parameters"
        ylabel = f"Percentage change in number of {metric} (%)"
        xlabel = f"Percentage change in parameter value (%)"
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Iterate through parameters, plotting this metric for each parameter
        for p in p_set:
            # Split the y calculation based on whether the parameter is targeted vaccination or not
            if p[0] != 'numTeams targeted':
                # ys is an array of y-values: percentage change in parameter value
                ys = [(float(p[i+1+4*j])/untargeted_base[i]-1)*100 for j in range(len(xs))]
            else:
                ys = [(float(p[i+1+4*j])/targeted_base[i]-1)*100 for j in range(len(xs))]

            # Set the y-scale based on the max y value
            curr_low, curr_high = plt.ylim()
            plt.ylim(min(-20,np.min(ys)-5,curr_low), max(20,np.max(ys)+5,curr_high))

            # Plot the values, with the label being the parameter name
            plt.plot(xs, ys, label=param_names[p[0]])

        plt.legend()
        plt.savefig(f"results/sensitivity/absN vaccStrat/{p[0]}_{metric}.pdf",bbox_inches='tight')
        #plt.show()
        plt.close()

# Input data from csv
data = open("results/sensitivity/monocentric.csv","r")

# Skip the header row. For each parameter set, plot the sensitivity chart
lineCount = 0
current_param_set = []
for line in data:
    lineCount += 1
    if lineCount == 1:
        continue
    lineContents = line.strip().split(",")
    if lineContents[0] == "":                           # If line empty, plot and reset parameter set
        plot_param_set(current_param_set)
        current_param_set = []
    else:                                               # Else plot this parameter set
        current_param_set.append(lineContents)
