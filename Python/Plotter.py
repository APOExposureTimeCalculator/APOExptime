import matplotlib.pyplot as plt
import numpy as np

def plotter(Observation):
    if Observation.isImager == 1:
        filterDATA = Observation.SN
        filterNum = len(filterDATA)
        filterSN = [row[0] for row in filterDATA]
        filter_names = [row.replace('prime_SN', "'") for row in [row[1] for row in filterDATA]]
        exptime = str(Observation.exptime)


        ind = np.arange(filterNum)  # the x locations for the groups
        width = 0.35  # the width of the bars
        barcolors = ["cornflowerblue", "green", "red", "firebrick", "maroon"]

        plt.figure(figsize=(10.0, 5.0))
        plt.bar(ind, filterSN, color=barcolors, width=0.35)
        plt.ylabel(r"$\frac{S}{N}$", fontsize=18, rotation=0, labelpad=20)
        plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
        plt.title(r"$\frac{S}{N}$  vs.  $\lambda$  for " + exptime + " seconds", y=1.08, fontsize=20)
        plt.xticks(ind, filter_names, fontsize=18)
        plt.yticks(fontsize=13)
        plt.ylim(0, 1.2 * np.max(filterSN))
        for i, v in enumerate(filterSN):
            plt.text(ind[i], int(v) + 1, str(int(v)), horizontalalignment="center", verticalalignment="bottom", fontsize=15)
        plt.grid(True)

        plt.show()

