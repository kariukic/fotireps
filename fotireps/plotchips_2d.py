#!/usr/bin/env python

import numpy as np
from numpy import ma
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import paramsfile


# Python3 code to take the output binary files from CHIPS and plot them with Bryna's plotting code. The code computes and applies the normalisations to take the units from Jy^2Hz^2sr^2 to mK^2 Mpc^3 and transposes and folds the data across the kpar = 0 line. The code will plot the output from two runs (or two pols)

# *********** IMPORTANT **********
# This is code appropriate for crosspower spectra, where the cross, flag and resid power are all differenced. Here, the flagpower and residpower do not represent noise power, but instead represent noise uncertainty. The weights can be used to compute the noise power in 2D, and also then used as the optimal estimator to get the noise uncertainty on the final 1D plot.


def compute_power(args):
    # grab parameters from paramsfile
    lowerfreq = float(paramsfile.lowerfreq)
    # wedge = float(paramsfile.wedge)
    wedge = args.wedge
    mx = float(paramsfile.mx)
    mn = float(paramsfile.mn)
    chan_width = float(paramsfile.chan_width)
    deltat = float(paramsfile.deltat)
    umax = float(paramsfile.umax)

    # print(mx, mn, lowerfreq)

    # set file sizes

    Nchan = int(args.N_chan)
    Nkperp = int(args.N_kperp)
    Neta = int(Nchan / 2)

    lperp = np.zeros(Nkperp)
    for i in range(0, Nkperp - 1):
        lperp[i] = float(i) * umax * 1.1 / float(Nkperp) * 2.0 * np.pi

    lperp[0] = lperp[1] / 2.0

    # read-in data

    filename = args.basedir + "crosspower_" + args.filestub + ".dat"
    train_xf = open(filename, "rb")
    crosspower = np.fromfile(train_xf, dtype=np.float32)
    crosspower = np.reshape(crosspower, (Nkperp, Nchan))
    # crosspower = np.transpose(crosspower, (1,0))

    filename = args.basedir + "outputweights_" + args.filestub + ".dat"
    train_w = open(filename, "rb")
    weights = np.fromfile(train_w, dtype=np.float32)
    weights = np.reshape(weights, (Nkperp, Nchan))
    # weights = np.transpose(weights, (1,0))

    # normalise by the weights

    crosspower = crosspower / (weights + 0.00001)
    crosspower = crosspower[0:Nkperp, Neta:Nchan]
    weights = weights[0:Nkperp, Neta:Nchan]

    # print(crosspower[0:Nkperp-1,0])
    # print(crosspower[0:Nkperp-1,20])

    # print(np.shape(weights))

    # *****************************

    # Define parameters

    # define cosmological units and conversions

    freq = np.zeros(Nchan)
    for i in range(0, Nchan - 1):
        freq[i] = lowerfreq + float(i) * chan_width

    hubble = 0.704 * 100.0  # km/s/Mpc
    c = 2.995e5  # km/s
    omega_matter = 0.272
    omega_baryon = 0.046
    omega_lambda = 0.7
    z = 6.8

    f21 = c * 1.0e3 / 0.21

    z = f21 / freq[Neta] - 1.0

    n = 0.5  # scale-invariant spectrum
    deltah = 1.94e-5 * omega_matter ** (-0.785 - 0.05 * np.log10(omega_matter)) * 0.93
    bw = float(Nchan) * chan_width  # Hz
    Ez = np.sqrt(omega_matter * (1.0 + z) ** 3 + omega_lambda)

    DM = 4424.7 + 224.7 * z  # transverse comoving distance z

    boltz = 1.38e-23
    D = 4.6  # m

    # lambd = 2.997e8/freq   #m
    # beam area calculated from int (beam**2) for each coarse channel

    beam_area_per_chan1 = 80.0e3 * 0.07597

    mpc_conversion = (
        DM ** 2 * 2.995e5 * (1.0 + z) ** 2 / 100.0 / f21 / Ez
    )  # Mpc**3/sr.Hz

    M = Nchan / 2

    normm = 1.0

    lamb = c * 1.0e3 / lowerfreq
    beam_area_steradian1 = beam_area_per_chan1 / 80.0e3
    Ae = 21.0

    obs_volume1 = (beam_area_steradian1) * chan_width * Nchan  # sr. Hz

    jy2_Nchan__to__jy2hz2 = chan_width ** 2 * Nchan
    jy2__to__K2sr2 = (lamb ** 2 / (2.0 * boltz * 1.0e26)) ** 2

    normalisation1 = (
        jy2_Nchan__to__jy2hz2 * jy2__to__K2sr2 / obs_volume1 * 1.0e6 * mpc_conversion
    )

    # For the weights, this mapping allows the weights to be equated with the noise **power**. Noise uncertainty in 2D is obtained by dividing the sigma**2 by sqrt(Nchanall). this then matches the amplitude of the flagpower and residpower noise **uncertainties**.

    normalisation_w1 = normalisation1

    Tsys = 220.0  # K

    expected_noise = (
        2.0 * boltz * Tsys * 1.0e26 / D ** 2 / np.sqrt(bw / Nchan * deltat)
    )  # Jy
    expected_noise = (expected_noise) ** 2  # square ->Jy**2

    factor = c * (1.0 + z) ** 2 / (2.0 * np.pi * 100.0 * f21 * Ez)

    # parametrize in h units - (/Mpc -> h /Mpc)
    hfactor = hubble / 100.0  # Hubble constant used in calcs.

    decoherence_factor = 2.0

    weight_data = weights / (normalisation_w1) ** 2 * 36.0 * np.sqrt(Nchan) / 2.0
    crosspower = crosspower * decoherence_factor * normalisation1

    eta = np.zeros(int(Neta))
    for i in range(0, int(Neta) - 1):
        eta[i] = float(i) - 0.5

    eta[0] = eta[1] / 2.0

    kpa = eta / bw / factor
    # print(eta[0],eta[1],eta[3],bw,factor)
    kper = lperp / DM

    conv_factor = 1.0 / (1.0 / DM * 2.0 * np.pi)
    ns_conv = eta / bw * 1.0e9

    # 1.06 here gives the horizon line for the maximum pointing
    wedgeedge = 1.06 * DM / 2.0 / np.pi / factor / freq[0]

    Nchancut = Nchan

    # delay_params = [ns_conv[2]-ns_conv[1] , ns_conv[Nchancut/2-1]]
    kpar_bin = kpa[2] - kpa[1]

    return crosspower, kper, kpa, Nkperp, Neta


## Plotting

# print(mn,mx)


def plot_power(crosspower, kper, kpa, Nkperp, Neta, outputmode, outdir):
    mx = float(paramsfile.mx)
    mn = float(paramsfile.mn)

    fig = plt.figure(figsize=(8, 6))
    left, bottom, width, height = 0.15, 0.15, 0.75, 0.75
    ax = fig.add_axes([left, bottom, width, height])

    # start, stop, n_values = np.log10(mn), np.log10(mx), 100

    x_vals = kper[2:Nkperp]
    y_vals = kpa[0 : Neta - 2]
    X, Y = np.meshgrid(x_vals, y_vals)

    # print(np.shape(crosspower[2:Nkperp,0:Neta-2]))
    # print(np.shape(X))

    # print(kper[2],kper[Nkperp-2])
    # print(kpa[0],kpa[1],kpa[Neta-2])

    crosspower = ma.masked_where(crosspower <= 0, crosspower)

    z = np.swapaxes((crosspower[2:Nkperp, 0 : Neta - 2]), 1, 0)

    setlen = int(np.log10(mx / mn) + 1)

    lev_exp = np.arange(np.floor(np.log10(z.min()) - 1), np.ceil(np.log10(z.max()) + 1))
    levs = np.power(10, lev_exp)

    cp = plt.contourf(
        X,
        Y,
        z,
        levels=np.logspace(np.log10(mn), np.log10(mx), 40, base=10),
        locator=ticker.LogLocator(),
        vmin=mn,
        vmax=mx,
        cmap=cm.Spectral_r,
    )
    cbar = plt.colorbar(cp, label=r"P(k) mK$^2$ $h^{-3}$ Mpc$^3$", format="%.0e")
    cbar.set_ticks(np.logspace(np.log10(mn), np.log10(mx), setlen))
    cbar.update_ticks()

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    ax.set_title("Crosspower (mK$^2$ $h^{-3}$ Mpc$^3$)", size=20)
    ax.set_xlabel(r"k$_\bot$ ($h$Mpc$^{-1}$)", fontsize=18)
    ax.set_ylabel(r"k$_\parallel$ ($h$Mpc$^{-1}$)", fontsize=18)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlim=[0.008, 0.4])
    ax.set(ylim=[0.005, kpa[Neta - 2]])

    if outputmode == "png":
        print("Printing to file")
        plt.savefig(outdir + "/crosspower.png")

    if outputmode == "screen":
        print("Printing to screen")
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filestub")
    parser.add_argument(
        "--basedir", default="/astro/mwaeor/MWA/output/", required=False
    )
    parser.add_argument("--N_kperp", type=int, default=80, required=False)
    parser.add_argument("--N_chan", type=int, default=384, required=False)
    parser.add_argument(
        "--output", default="png", required=False, help="mode ('screen','png')"
    )
    parser.add_argument("--outputdir", default=".", required=False)
    parser.add_argument("--wedge", default=0, required=False)
    args = parser.parse_args()
    outputmode = args.output
    outdir = args.outputdir

    # easy way to define args
    # args = Namespace(
    #     basedir="/astro/mwaeor/MWA/output/", #this path is the default
    #     filestub=yy_0.iter.FILETAG1_FILETAG2" #keep the "yy_0.iter" part as it is. maybe only change yy to xx
    #     N_kperp=int(80), #default
    #     N_chan=int(384), #default
    #     wedge=wedge, #either 0 or 3.5 for all modes and window modes only respectively. #invalid for 2d plot anyway.
    #     output = "plot type"
    #     outputdir = "path_to_place_the_plot"
    # )

    # print(args.echo)
    crosspower, kper, kpa, Nkperp, Neta = compute_power(args)
    plot_power(crosspower, kper, kpa, Nkperp, Neta, outputmode, outdir)
