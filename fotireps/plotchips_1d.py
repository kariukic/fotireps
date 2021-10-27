#!/usr/bin/env python

import numpy as np
from numpy import ma
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
    omega_matter = float(paramsfile.omega_matter)
    omega_baryon = float(paramsfile.omega_baryon)
    omega_lambda = float(paramsfile.omega_lambda)

    mapval = cm.viridis
    print(mapval)

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

    filename = args.basedir + "residpower_" + args.filestub + ".dat"
    train_xf = open(filename, "rb")
    residpower = np.fromfile(train_xf, dtype=np.float32)
    residpower = np.reshape(residpower, (Nkperp, Nchan))

    filename = args.basedir + "outputweights_" + args.filestub + ".dat"
    train_w = open(filename, "rb")
    weights = np.fromfile(train_w, dtype=np.float32)
    weights = np.reshape(weights, (Nkperp, Nchan))

    # normalise by the weights

    crosspower = crosspower / (weights + 1.0e-8)
    residpower = residpower / (weights + 1.0e-8)

    # cut negative frequencies

    crosspower = crosspower[0:Nkperp, Neta:Nchan]
    residpower = residpower[0:Nkperp, Neta:Nchan]
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

    f21 = c * 1.0e3 / 0.21

    z = f21 / freq[Neta] - 1.0

    n = 0.5  # scale-invariant spectrum
    deltah = 1.94e-5 * omega_matter ** (-0.785 - 0.05 * np.log10(omega_matter)) * 0.93
    bw = Nchan * chan_width  # Hz
    Ez = np.sqrt(omega_matter * (1.0 + z) ** 3 + omega_lambda)

    DM = 4424.7 + 224.7 * z  # transverse comoving distance to redshift z

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

    # multiply powers by normalisations to get the units right

    weight_data = weights / (normalisation_w1) ** 2 * 36.0 * np.sqrt(Nchan) / 2.0
    crosspower = crosspower * decoherence_factor * normalisation1
    residpower = residpower * decoherence_factor * normalisation1

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

    ## Average data to 1D including a wedge cut

    Netaa = 300.0

    kmax = np.sqrt(kper[Nkperp - 1] ** 2 + kpa[Neta - 1] ** 2)

    ptot = np.zeros(int(Netaa))
    psig = np.zeros(int(Netaa))
    pmeas = np.zeros(int(Netaa))
    num = np.zeros(int(Netaa))
    numedge = np.zeros(int(Netaa))

    ktot_bins = np.zeros(int(Netaa + 1))
    for i in range(0, int(Netaa)):
        ktot_bins[i] = float(i) * 1.7163922 / float(97.0 / 4.0)

    noise_obs = np.sqrt(weight_data)

    for i in range(1, Nkperp - 1):
        for j in range(0, Neta - 1):

            for k in range(0, int(Netaa) - 1):

                if (
                    np.sqrt(kper[i] ** 2 + kpa[j] ** 2) > ktot_bins[k]
                    and np.sqrt(kper[i] ** 2 + kpa[j] ** 2) <= ktot_bins[k + 1]
                ):
                    if (
                        kpa[j] > kper[i] * wedge
                        and kper[i] <= 10.8
                        and kper[i] > 0.00
                        and kpa[j] > 0.0
                    ):

                        ptot[k] = ptot[k] + (noise_obs[i, j]) ** 2
                        pmeas[k] = (
                            pmeas[k] + (crosspower[i, j]) * (noise_obs[i, j]) ** 2
                        )
                        numedge[k] = numedge[k] + 1.0

    pmeas = pmeas / ptot

    pmeas = pmeas * ktot_bins[0 : int(Netaa)] ** 3 / 2.0 / np.pi ** 2

    kernel = 2.0

    ptot = (ktot_bins[0 : int(Netaa)] ** 3 / np.sqrt(ptot) / 2.0 / np.pi ** 2) * kernel

    # print(pmeas)
    return pmeas, ptot, ktot_bins, Netaa, wedge


##############
## Plotting ##
##############

# print(mn,mx)


def plot_power(pmeas, ptot, ktot_bins, Netaa, outputmode, outdir, wedge=None):
    fig = plt.figure(figsize=(8, 6))
    left, bottom, width, height = 0.15, 0.15, 0.75, 0.75
    ax = fig.add_axes([left, bottom, width, height])

    # start, stop, n_values = np.log10(mn), np.log10(mx), 100

    # print(np.shape(crosspower[2:Nkperp,0:Neta-2]))
    # print(np.shape(X))

    # print(kper[2],kper[Nkperp-2])
    # print(kpa[0],kpa[1],kpa[Neta-2])

    plt.plot(
        ktot_bins[0 : int(Netaa)], pmeas, "r--", ktot_bins[0 : int(Netaa)], ptot, "bs"
    )

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    ax.set_title("Power (mK$^2$)", size=20)
    ax.set_xlabel(r"k ($h$Mpc$^{-1}$)", fontsize=18)
    ax.set_ylabel(r"$\Delta^2$ ($h$mK$^{2}$)", fontsize=18)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlim=[0.008, 2.0])
    ax.set(ylim=[1.0e2, 1.0e9])

    plotname = "/1dpower_window_only.png" if wedge == 3.5 else "/1dpower_allmodes.png"

    if outputmode == "png":
        print("Printing to file")
        plt.savefig(outdir + plotname)

    if outputmode == "screen":
        print("Printing to screen")
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("basedir")
    parser.add_argument("filestub")
    parser.add_argument("N_kperp", type=int)
    parser.add_argument("N_chan", type=int)
    parser.add_argument("output", help="mode ('screen','png')")
    parser.add_argument("outputdir")
    parser.add_argument("wedge")
    args = parser.parse_args()
    outputmode = args.output
    outdir = args.outputdir

    # easy way to define args
    # args = Namespace(
    #     basedir="/astro/mwaeor/MWA/output/", #this path is the default
    #     filestub=yy_0.iter.FILETAG1_FILETAG2" #keep the "yy_0.iter" part as it is. maybe only change yy to xx
    #     N_kperp=int(80), #default
    #     N_chan=int(384), #default
    #     wedge=wedge, #either 0 or 3.5 for all modes and window modes only respectively
    #     output = "plot type"
    #     outputdir = "path_to_place_the_plot"
    # )

    # print(args.echo)
    pmeas, ptot, ktot_bins, Netaa, wedge = compute_power(args)
    plot_power(pmeas, ptot, ktot_bins, Netaa, outputmode, outdir, wedge)
