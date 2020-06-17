import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
base = [178220,5677,724477,8772]                  # The base results set

param_names = {
    "":""
}

def plot_param_set(p_set):
    xs = [-50,-30,-20,-10,10,20,30,50]
    xlabel = f"Percentage change in parameter"
    metrics = ["cases","deaths","vaccines","drone deliveries"]

    print(p_set)

    for i in range(len(metrics)):
        metric = metrics[i]
        title = f"Sensitivity of {metric} to changes in parameters"
        ylabel = f"Percentage change in {metric} (%)"
        xlabel = f"Percentage change in parameter value (%)"
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Iterate through parameters, plotting this metric for each parameter
        for p in p_set:
            # ys is an array of y-values: percentage change in parameter value
            ys = [(float(p[i+1+4*j])/base[i]-1)*100 for j in range(len(xs))]
            print(ys)

            # Plot the values, with the label being the parameter name
            plt.plot(xs, ys, label=p[0])

        #plt.savefig(f"results/sensitivity/{metric}.pdf",bbox_inches='tight')
        plt.legend()
        plt.show()


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
