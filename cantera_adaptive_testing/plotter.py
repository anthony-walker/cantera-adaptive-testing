import matplotlib as mpl
import matplotlib.pyplot as plt
import ruamel.yaml, warnings, os, operator
import numpy as np
import random
# change font
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'


def plot_precon_species_barchart(labels, y, xend, manual_max=None):
    # make plots
    x = np.arange(len(labels))  # the label locations
    width = 0.5  # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(x, y, width, color='#7570b3', label="speed-up")
    # labels and ticks
    if manual_max is None:
        ymax =np.ceil(max(y) + max(y) * 0.05)
    else:
        ymax = manual_max
    yticks = np.linspace(0, ymax, 5)
    if ymax > 1:
        ylbls, yticks = zip(*[(str(int(yt)), int(yt)) for yt in yticks])
    else:
        ylbls, yticks = zip(*[(str(yt), yt) for yt in yticks])
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylbls, fontsize=14)
    trim = list(zip(x, labels))
    newlabels = [(trim[i][0], trim[i][1]) for i in range(0, len(trim)-xend, 2)]
    newlabels += [(trim[i][0], trim[i][1]) for i in range(len(trim)-xend, len(trim), 1)]
    newx, newlabels = zip(*newlabels)
    ax.set_xticks(newx)
    ax.set_xticklabels(newlabels, fontsize=14)
    ax.set_xlabel('Threshold', fontsize=14)
    return fig, ax


def plot_single_comparison(self, x1, y1, x2, y2, xL, yL, test_name, class_name, save_dir=None, replace=False):
    # plot the data
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x1, y1, color='tab:blue')
    ax.plot(x2, y2, color='tab:orange', linestyle=":")
    ax.set_ylabel(yL)
    ax.set_xlabel(xL)
    plt.legend(["Preconditioned", "Standard"])
    joinDir = self.figDir if save_dir is None else save_dir
    save_name = os.path.join(joinDir, "{:s}-vs-{:s}-{:s}-{:s}".format(xL, yL, test_name, class_name))
    if not replace:
        ctr = 1
        while os.path.exists(save_name+"-{:d}.pdf".format(ctr)):
            ctr += 1
        plt.savefig(save_name+"-{:d}.pdf".format(ctr))
    else:
        plt.savefig(save_name+".pdf")
    plt.close()


def plot_multi_comparison(self, x1, y1, x2, y2, xL, yL, test_name, class_name, leg_labs=[]):
    # plot the data
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    rows, cols = np.shape(y1)
    for i in range(cols):
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]
        ax.plot(x1, y1[:, i], color=color[0])
        ax.plot(x2, y2[:, i], color=color[0], linestyle=":")
    # plt.show()
    ax.set_ylabel(yL)
    ax.set_xlabel(xL)
    if leg_labs:
        plt.legend(leg_labs)
    save_name = "./figures/{:s}-vs-{:s}-{:s}-{:s}".format(xL, yL, test_name, class_name)
    ctr = 1
    while os.path.exists(save_name+"-{:d}.pdf".format(ctr)):
        ctr += 1
    plt.savefig(save_name+"-{:d}.pdf".format(ctr))
    plt.close()


def plot_species(self, reactor, x1, y1, x2, y2, xL, yL, test_name, class_name):
    dname = str.join("-", (class_name, test_name, "species-plots"))
    dname = os.path.join(self.figDir, dname)
    rows, cols = np.shape(y1)
    if not os.path.isdir(dname):
        os.mkdir(dname)
    for i in range(cols):
        if (np.mean(y1[:, i]) > self.net1.rtol):
            name = reactor.component_name(i+1)
            plot_single_comparison(self, x1, y1[:, i], x2,  y2[:, i], xL, "MF-{:s}".format(name), test_name+"-{:s}".format(name), class_name, save_dir=dname, replace=True)

