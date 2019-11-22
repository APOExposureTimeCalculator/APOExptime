import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import rgb2hex
from matplotlib import cm


def makeplots(Observation):
    """Function to visualize the output of the ``APOinputclasses.Observation`` class.

    This function has different outputs based on whether the instrument used is a spectrograph or an imager. For an
    imager, the function will output a bar plot that shows the S/N per bandpass and a plot that shows the different
    noise sources as a function of wavelength. If it is a spectrograph, then the function will plot the S/N per
    wavelength and it also plots the noise sources as a function of wavelength.

    Parameters
    -----------
    Observation : object
        The ``APOinputclasses.Observation``class.

    """
    if Observation.isImager == 1:

        # Set colors
        cmap = cm.get_cmap("rainbow", filterNum)
        hexvals = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            hexvals.append(rgb2hex(rgb))

        width = 0.35  # width of SN/exptime plot bars
        width2 = 0.3  # width of noise plot bars

        if Observation[0] == 1:
            filterDATA = Observation.SN
            filterSN = [row[0] for row in filterDATA]
            filterNum = len(filterDATA)
            filter_names = [row.replace('prime_SN', "'") for row in [row[1] for row in filterDATA]]
            exptime = str(Observation.exptime)
            ind = np.arange(filterNum)  # the x locations for the bars

            plt.figure(figsize=(10.0, 5.0))
            plt.bar(ind, filterSN, color=hexvals, edgecolor="black", width=width)
            plt.title(r"$\frac{S}{N}$  vs.  $\lambda$  for " + exptime + " seconds", y=1.08, fontsize=20)
            plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            plt.ylabel(r"$\frac{S}{N}$", fontsize=18, rotation=0, labelpad=20)
            plt.xticks(ind, filter_names, fontsize=18)
            plt.yticks(fontsize=13)
            plt.ylim(0, 1.2 * np.max(filterSN))
            plt.grid(True)
            for i, v in enumerate(filterSN):
                plt.text(ind[i], int(v) + 1, str(int(v)), horizontalalignment="center", verticalalignment="bottom",
                         fontsize=15)
            plotname = "SNfromTime" + exptime + "s.png"
            plotname2 = "Noise" + plotname[2:]
            plt.savefig(plotname)

        if Observation[0] == 0:
            ind = np.arange(filterNum)  # the x locations for the groups

            plt.figure(figsize=(10.0, 5.0))
            plt.bar(ind, filterTime, color=hexvals, edgecolor="black", width=width)
            plt.title(r"Exposure Time  vs.  $\lambda$  for $\frac{S}{N}=$" + SN, y=1.08, fontsize=20)
            plt.ylabel(r"$t$ ($s$)", fontsize=18, rotation=0, labelpad=25)
            plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            plt.xticks(ind, filter_names, fontsize=18)
            plt.yticks(fontsize=13)
            plt.ylim(0, 1.2 * np.max(filterTime))
            plt.grid(True)
            for i, v in enumerate(filterSN):
                plt.text(ind[i], int(v) + 1, str(int(v)), horizontalalignment="center", verticalalignment="bottom",
                         fontsize=15)
            plotname = "TimefromSN" + SN + ".png"
            plotname2 = "Noise" + plotname[4:]
            plt.savefig(plotname)

        else:
            plotname2 = "Noise.png"
            print("unknown method error")

        sourcenoise = [10, 25, 30, 20, 5.0]
        bgnoise = [15, 20, 27, 30, 45]
        rdnoise = [3.0, 7.0, 5.0, 6.0, 7.0]

        plt.subplots(figsize=(10.0, 5.0))

        bar1 = plt.bar(ind - width2, sourcenoise, color=hexvals, edgecolor="black", hatch="\\", width=width2,
                       label="Source")
        bar2 = plt.bar(ind, bgnoise, color=hexvals, edgecolor="black", hatch="+", width=width2, label="Background")
        bar3 = plt.bar(ind + width2, rdnoise, color=hexvals, edgecolor="black", hatch="//", width=width2,
                       label="Readout")

        plt.ylabel(r"Arbitrary", fontsize=18, labelpad=20)
        plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
        plt.title("Noise Sources", y=1.08, fontsize=20)
        plt.xticks(ind, ["u'", "g'", "r'", "i'", "z'"], fontsize=18)
        plt.yticks(fontsize=13)
        plt.ylim(0, 1.2 * np.max([np.max(sourcenoise), np.max(bgnoise), np.max(rdnoise)]))
        for i, v in enumerate(sourcenoise):
            plt.text(ind[i] - width2, v + 1, str(v), horizontalalignment="center", verticalalignment="bottom",
                     fontsize=15)
        for i, v in enumerate(bgnoise):
            plt.text(ind[i], v + 1, str(v), horizontalalignment="center", verticalalignment="bottom", fontsize=15)
        for i, v in enumerate(rdnoise):
            plt.text(ind[i] + width2, v + 1, str(v), horizontalalignment="center", verticalalignment="bottom",
                     fontsize=15)
        plt.grid(True)
        plt.legend(loc="upper center", bbox_to_anchor=(1.16, 0.98), fontsize=15, shadow=True)
        plt.savefig(plotname2)

    if Observation.isImager == 0:
        if Observation[0] == 1:
            plt.figure(figsize=(10.0, 5.0))
            plt.plot(x, y1, color="mediumseagreen", label=r"$\frac{S}{N}$")
            plt.title(r"$\frac{S}{N}$  vs.  $\lambda$  for " + exptime + " seconds", y=1.08, fontsize=20)
            plt.ylim(-1, 1.5 * np.median(y1))
            plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            plt.ylabel(r"$\frac{S}{N}$", fontsize=18, rotation=0, labelpad=20)
            plt.xticks(fontsize=13)
            plt.yticks(fontsize=13)
            plt.grid(True)
            plt.show()

        if Observation[0] == 0:
            plt.figure(figsize=(10.0, 5.0))
            plt.plot(x, y2, color="mediumslateblue", label="Exposure Time")
            plt.title(r"Exposure Time vs.  $\lambda$  for  $\frac{S}{N}=$" + SN, y=1.08, fontsize=20)
            plt.ylim(-1, 1.5 * np.median(y1))
            plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            plt.ylabel(r"$t$ ($s$)", fontsize=18, rotation=0, labelpad=25)
            plt.xticks(fontsize=13)
            plt.yticks(fontsize=13)
            plt.grid(True)
            plt.show()

        else:
            plotname2 = "Noise.png"
            print("unknown method error")

        plt.figure(figsize=(10.0, 5.0))
        plt.plot(x, y1, color="green", label=r"Source")
        plt.plot(x, y2, color="blue", linestyle="--", label=r"Background")
        plt.plot(x, y3, color="black", linestyle="-.", label=r"Readout")
        plt.title("Noise Sources", y=1.08, fontsize=20)
        plt.ylim(-1, 1.5 * np.max([np.median(y1), np.median(y2), np.median(y3)]))
        plt.xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
        plt.ylabel(r"Arbitrary", fontsize=18, labelpad=20)
        plt.legend(loc="upper center", bbox_to_anchor=(1.16, 0.98), fontsize=15, shadow=True)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.grid(True)
        plt.show()

    else:
        print("unknown instrument error")
