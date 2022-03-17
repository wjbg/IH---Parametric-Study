# Scripts to generate figures from matlab data.
import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np

def new_ticks(A):
    loc = np.arange(len(A))
    labels = ["%.2f" % a for a in A]
    return loc, labels

# Baseline values
mat = loadmat("base_values.mat")
baseline = {"I": float(mat["I0"][0][0]),
            "delta": mat["delta0"][0][0],
            "freq": float(mat["f0"][0][0]),
            "mu_r": float(mat["mu_r0"][0][0]),
            "sigma": float(mat["sigma0"][0][0]),
            "P": mat["power0"][0][0],
            "flux": mat["power0"][0][0]}
mu_0 = mat["mu0_const"][0][0]
t_plate = 5E-3

# Frequency
mat = loadmat("frequency.mat")
freq = mat["freq"][0]
P = mat["induced_power"][0]
flux = mat["flux_in_coil"][0]
delta_over_t = mat["delta"][0]/t_plate

with plt.style.context("report_style.mplstyle"):
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twiny()
    fig.subplots_adjust(bottom=0.3)

    scatter = ax1.loglog(freq, P/baseline["P"], "ko", markersize=4)
    ax1.grid(visible=True, which="major")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("Normalized induced power [-]")

    f_left = np.polyfit(np.log10(freq[:3]),
                        np.log10(P[:3]/baseline["P"]), 1)
    f_right = np.polyfit(np.log10(freq[-3:]),
                         np.log10(P[-3:]/baseline["P"]), 1)
    ax1.plot(np.array([8, 350]),
             np.array([8, 350])**f_left[0] * 10**f_left[1],
             "k", linewidth=0.5)
    ax1.text(10, 0.0008, r"$\angle \approx 2$")
    ax1.plot(np.array([4E5, 1.2E7]),
             np.array([4E5, 1.2E7])**f_right[0] * 10**f_right[1],
             "k", linewidth=0.5)
    ax1.text(1E6, 1.05, r"$\angle \approx \frac{1}{2}$ ")
    ylim = ax1.get_ylim()

    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.25))
    ax2.set_frame_on(True)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xscale("log")
    ax2.plot(delta_over_t, P/baseline["P"], linestyle="None")
    ax2.invert_xaxis()
    ax2.set_xlabel(r"$\delta / t$")

    plt.savefig("../img/frequency.png", dpi=600)


# Conductivity
mat = loadmat("conductivity.mat")
sigma = mat["sigma"][0]
P = mat["induced_power"][0]
flux = mat["flux_in_coil"][0]
delta_over_t = mat["delta"][0]/t_plate

with plt.style.context("report_style.mplstyle"):
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twiny()
    fig.subplots_adjust(bottom=0.3)

    scatter = ax1.loglog(sigma, P/baseline["P"], "ko", markersize=4)
    ax1.grid(visible=True, which="major")
    ax1.set_xlabel("Conductivity [S/m]")
    ax1.set_ylabel("Normalized induced power [-]")
    ax1.set_ylim(ylim)

    f_left = np.polyfit(np.log10(sigma[:5]),
                        np.log10(P[:5]/baseline["P"]), 1)
    f_right = np.polyfit(np.log10(sigma[-3:]),
                         np.log10(P[-3:]/baseline["P"]), 1)
    ax1.plot(np.array([900, 15E3]),
             np.array([900, 15E3])**f_left[0] * 10**f_left[1],
             "k", linewidth=0.5)
    ax1.text(4500, 0.01, r"$\angle \approx 1$")
    ax1.plot(np.array([5E7, 1.2E9]),
             np.array([5E7, 1.2E9])**f_right[0] * 10**f_right[1],
             "k", linewidth=0.5)
    ax1.text(1E8, 0.1, r"$\angle \approx -\frac{1}{2}$ ")

    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.25))
    ax2.set_frame_on(True)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xscale("log")
    ax2.plot(delta_over_t, P/baseline["P"], linestyle="None")
    ax2.invert_xaxis()
    ax2.set_xlabel(r"$\delta / t$")

    plt.savefig("../img/conductivity.png", dpi=600)

# Current
mat = loadmat("current.mat")
I = mat["I"][0]
P = mat["induced_power"][0]
flux = mat["flux_in_coil"][0]
delta_over_t = mat["delta"][0]/t_plate

with plt.style.context("report_style.mplstyle"):
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    fig.subplots_adjust(bottom=0.15)

    scatter = ax1.loglog(I, P/baseline["P"], "ko", markersize=4)
    ax1.grid(visible=True, which="major")
    ax1.set_xlabel("Current [A]")
    ax1.set_ylabel("Normalized induced power [-]")

    f = np.polyfit(np.log10(I),
                   np.log10(P/baseline["P"]), 1)
    ax1.plot(np.array([0.05, 2000]),
             np.array([0.05, 2000])**f[0] * 10**f[1],
             "k", linewidth=0.5)
    ax1.text(3, 2, r"$\angle = 2$")

    plt.savefig("../img/current.png", dpi=600)
