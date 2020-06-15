import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150
base = [178220,5677,724477,8772]                  # The base results set

def plot_param_set(p_set):
    xs = [-50,-30,-20,-10,10,20,30,50]
    xt = ["-50%","-30%","-20%","-10%","+10%","+20%","+30%","+50%"]
    xlabel = f"Percentage change in parameter"
    metrics = ["cases","deaths","vaccines","drone deliveries"]

    for i in range(len(metrics)):
        metric = metrics[i]
        title = f"Sensitivity of {metric} to changes in parameters"
        ylabel = f"Percentage change in {metric}"
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(xt)

        # Iterate through parameters, plotting this metric for each parameter
        for p in p_set:
            ys = [p[i+1+4*j]/base[i] for j in range(len(xs))]
            print(ys)
            # Ensure that the correct values are selected
            quit()
            plt.plot(xs, ys)

        #plt.savefig(f"results/sensitivity/{metric}.pdf",bbox_inches='tight')
        plt.show()


# Input data from csv
data = open("results/sensitivity/monocentric.csv","r")

# Skip the header row. For each parameter set, plot the sensitivity chart
lineCount = 0
for line in data:
    lineCount += 1
    if lineCount == 1:
        continue
    current_param_set = []
    if line == "":                          # If line empty, plot previous param set
        plot_param_set(current_param_set)
        current_param_set = []
    else:
        lineContents = line.strip().split(",")
        current_param_set.append(lineContents)
