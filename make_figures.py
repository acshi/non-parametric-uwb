#!/usr/bin/python3
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math
import csv

plt.rcParams.update({'font.size': 14})

sweeps = ["by-measurements", "by-delay-sigma", "by-delay-mean"]
sweep_labels = ["Number of measurements per edge", r"$\sigma_{\mathrm{delay}}$ ($\mathrm{\mu s}$)", r"$\mu_{\mathrm{delay}}$ ($\mathrm{\mu s}$)"]
delay_means = [0, 0.51, 0.511, 0.512, 0.513, 0.514, 0.515, 0.516, 0.517, 0.518, 0.519, 0.52, 0.521, 0.522]
delay_sigmas = [3.3e-4*128, 3.3e-4*64, 3.3e-4*32, 3.3e-4*16, 3.3e-4*8, 3.3e-4*4, 3.3e-4*2, 3.3e-4, 3.3e-4/2, 3.3e-4/4]
x_variables = [[], delay_sigmas, delay_means]

evaluations = ["evaluate-los", "evaluate-nlos", "evaluate-nlos2", "evaluate-nlos3"]
evaluation_labels = ["Dataset 1", "Dataset 2", "Dataset 3", "Dataset 4"]

for sweep_i in range(len(sweeps)):
    sweep = sweeps[sweep_i]
    plt.clf()
    for evaluation_i in range(len(evaluations)):
        evaluation = evaluations[evaluation_i]
        data = np.genfromtxt("evaluations/{}_{}.csv".format(evaluation, sweep))
        V = data[:, 0]
        rmse = data[:, 1]

        uniq_vs = np.unique(V)
        rmse_by_v = [rmse[V == v] for v in uniq_vs]

        plt.ylabel("RMSE (m)")
        plt.xlabel(sweep_labels[sweep_i])

        if sweep == "by-measurements":
            plt.xscale("log")
            plt.violinplot(rmse_by_v, uniq_vs, widths=[n/3 for n in uniq_vs], showmeans=True)
            plt.ylim(ymax = 0.5)
        else:
            x_vars = [x_variables[sweep_i][int(v)] for v in uniq_vs]
            if sweep == "by-delay-sigma":
                plt.xscale("log")
                plt.violinplot(rmse_by_v, x_vars, widths=[n/3 for n in x_vars], showmeans=True)
            else:
                plt.violinplot(rmse_by_v, x_vars, widths=(np.max(x_vars) - np.min(x_vars)) / len(x_vars), showmeans=True)

    colors = ['blue', 'orange', 'green', 'red'] #, 'purple', 'brown', 'violet']
    fake_handles = [mlines.Line2D([], [], color=c) for c in colors]
    if sweep == "by-delay-sigma":
        plt.axvline(3.3e-4, color="gray", linestyle = "dashed")
        plt.legend(fake_handles, evaluation_labels, loc="upper left")
    else:
        plt.legend(fake_handles, evaluation_labels)
    if sweep == "by-delay-mean":
        plt.axvline(0.516, color="gray")
        plt.axvline(0.516 + 3.3e-4, linestyle = "dashed", color="gray")
        plt.axvline(0.516 - 3.3e-4, linestyle = "dashed", color="gray")
    plt.tight_layout()
    # plt.show()
    plt.savefig("figures/{}.pdf".format(sweep), bbox_inches = "tight", pad_inches = 0)

plt.clf()
scatter_xs = []
scatter_ys = []
for evaluation_i in range(len(evaluations)):
    evaluation = evaluations[evaluation_i]
    data = np.genfromtxt("evaluations/{}_by-delay-mean.csv".format(evaluation))
    V = data[:, 0]
    rmse = data[:, 1]
    delays = data[:, 2:]
    for j in range(len(V)):
        scatter_xs += [delay_means[int(V[j])]] * len(delays[j])
        scatter_ys = np.append(scatter_ys, delays[j])
uniq_vs = np.unique(scatter_xs)
rmse_by_v = [scatter_ys[scatter_xs == v] for v in uniq_vs]
plt.violinplot(rmse_by_v, uniq_vs, widths=(0.523-0.509) / len(uniq_vs), showmeans=True)
plt.plot([0.51, 0.522], [0.51, 0.522], 'r-', linestyle = "dotted")
plt.plot([0.51, 0.522], [0.51 + 3.3e-4, 0.522 + 3.3e-4], 'r-', linestyle = "dashed")
plt.plot([0.51, 0.522], [0.51 - 3.3e-4, 0.522 - 3.3e-4], 'r-', linestyle = "dashed")
plt.xlabel("Prior $\mu_{\mathrm{delay}}$ ($\mathrm{\mu s}$)")
plt.ylabel("Posterior delays ($\mathrm{\mu s}$)")
plt.xlim(0.509, 0.523)
plt.ylim(0.507, 0.525)
plt.axhline(0.516, color="gray")
plt.tight_layout()
# plt.show()
plt.savefig("figures/delay_mean_recovery.pdf", bbox_inches = "tight", pad_inches = 0)

plt.clf()
methods = ["", "--first-peak", "--max-peak", "--triangulation", "--fixed-delays"]
method_labels = ["Our method", "First peak", "Max peak", "Triangulation", "Fixed delays"]
for method_i in range(len(methods)):
    method = methods[method_i]
    data = []
    for evaluation_i in range(len(evaluations)):
        evaluation = evaluations[evaluation_i]
        data += [np.genfromtxt("evaluations/{}_{}.csv".format(evaluation, method[2:]))]
    data = np.array(data)
    plt.violinplot(data.T, np.add([1, 2, 3, 4], (method_i-2)*0.13), widths=0.18, showmeans=True)
plt.ylim(ymax = 0.45, ymin = 0)
plt.ylabel("RMSE (m)")
plt.xlabel("Data set")
colors = ['blue', 'orange', 'green', 'red', 'purple']
fake_handles = [mlines.Line2D([], [], color=c) for c in colors]
plt.legend(fake_handles, method_labels, loc="upper left")
plt.yticks(np.arange(0, 0.41, 0.1))
plt.xticks(range(1, 5))
plt.tight_layout()
plt.savefig("figures/method_comparison.pdf", bbox_inches = "tight", pad_inches = 0)
# plt.show()

# Here we make figures of difficult measurement PDFs
pdf_files = ["./datasets/dataset4.txt", "./datasets/dataset4.txt", "./datasets/dataset5.txt"]
pdf_as = [81, 81, 54]
pdf_bs = [14, 12, 8]
pdf_is = [0, 1, 2]
plt.clf()
plt.rcParams.update({'font.size': 10})
for pdf_i_i in range(len(pdf_is)):
    pdf_i = pdf_is[pdf_i_i]
    pdf_tprops = []
    pdf_a = pdf_as[pdf_i]
    pdf_b = pdf_bs[pdf_i]
    with open(pdf_files[pdf_i]) as f:
        for line in f:
            cols = line.split()
            tprop = float(cols[0])
            a = int(cols[1])
            b = int(cols[2])

            if a == pdf_a and b == pdf_b:
                pdf_tprops += [tprop]
    print("Found {} measurements for {} {} in {}".format(len(pdf_tprops), pdf_a, pdf_b, pdf_files[pdf_i]))
    pdf_tprops = np.array(pdf_tprops)

    done = False
    while not done:
        std = pdf_tprops[:].std()
        median = np.median(pdf_tprops)
        outliers = np.arange(0, len(pdf_tprops))[np.logical_or(pdf_tprops < 0.5, np.abs(pdf_tprops - median) / std > 4)]
        pdf_tprops = np.delete(pdf_tprops, outliers)
        print("Removed {} outliers".format(len(outliers)))
        done = len(outliers) == 0

    sigma = 1.3e-4
    inv_sigma = 1.0 / sigma

    min_tprop = min(pdf_tprops)
    max_tprop = max(pdf_tprops)
    min_plot = min_tprop - 4 * sigma
    max_plot = max_tprop + 4 * sigma
    plt.subplot(2, len(pdf_is), pdf_i_i+1)
    plt.xlim(min_plot, max_plot)
    plt.hist(pdf_tprops, bins = 100)
    if pdf_i_i == 0:
        plt.ylabel("Counts")
    if True or pdf_i_i == 2:
        ticks = 3 if pdf_i_i != 1 else 2
        tick_spacing = int((max_plot - min_plot) / (ticks - 0.5) / 0.001) * 0.001
        if tick_spacing == 0:
            tick_spacing = 0.001
        first_tick = int(((max_tprop + min_tprop) / 2 - tick_spacing) / tick_spacing) * tick_spacing
        while first_tick < min_plot:
            first_tick += 0.001
        plt.xticks(first_tick + np.arange(ticks) * tick_spacing)

    # convolve with gaussian to get our estimated pdf
    min_x = pdf_tprops.min() - 6 * sigma
    max_x = pdf_tprops.max() + 6 * sigma
    spacing = 2.5e-6
    n = int((max_x - min_x) / spacing + 0.5) + 1
    pdf_xs = np.zeros(n)
    pdf_ys = np.zeros(n)
    y_sum = 0
    print("Evaluating convolution at {} points".format(n))
    for j in range(n):
        x = min_x + spacing * j
        a = (x - pdf_tprops) * inv_sigma
        y = np.exp(-0.5 * a * a).sum()
        y_sum += y
        pdf_xs[j] = x
        pdf_ys[j] = y
    plt.subplot(2, len(pdf_is), len(pdf_is)+pdf_i_i+1)
    plt.xlim(min_plot, max_plot)
    plt.xlabel("$t_{prop}$ ($\mathrm{\mu s}$)")
    if pdf_i_i == 0:
        plt.ylabel("PDF")
    if True or pdf_i_i == 2:
        ticks = 3 if pdf_i_i != 1 else 2
        tick_spacing = int((max_plot - min_plot) / (ticks - 0.5) / 0.001) * 0.001
        if tick_spacing == 0:
            tick_spacing = 0.001
        first_tick = int(((max_tprop + min_tprop) / 2 - tick_spacing) / tick_spacing) * tick_spacing
        while first_tick < min_plot:
            first_tick += 0.001
        plt.xticks(first_tick + np.arange(ticks) * tick_spacing)
    plt.plot(pdf_xs, pdf_ys)
plt.tight_layout()
# plt.show()
plt.savefig("figures/pdfs.pdf", bbox_inches = "tight", pad_inches = 0)
