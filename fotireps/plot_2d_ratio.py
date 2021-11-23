import numpy as np
import matplotlib.pyplot as plt
from argparse import Namespace

import paramsfile
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import cmocean

from plotchips_2d import compute_power as power2d

mx = float(paramsfile.mx)
mn = float(paramsfile.mn)
Nchan = int(384)
Nkperp = int(80)
Neta = int(Nchan / 2)


def colorbar(mappable):
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.07)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


def get_2dpower(FILETAG1, FILETAG2, pol="yy"):
    args = Namespace(
        basedir="/astro/mwaeor/MWA/output/",
        filestub=f"{pol}_0.iter.{FILETAG1}_{FILETAG2}",
        N_kperp=int(80),
        N_chan=int(384),
        wedge=0,
    )
    crosspawa, kper, kpa, _Nkperp, _Neta = power2d(args)

    return crosspawa, kper, kpa


def plot_ratio(power_list2d, kper, kpa, titles, plotname=None):
    fig, axlist = plt.subplots(1, 3, sharey=False, figsize=(18, 6), dpi=150)
    for i, (ax, crosspower) in enumerate(zip(axlist, power_list2d)):
        x_vals = kper[2:Nkperp]
        y_vals = kpa[0 : Neta - 2]
        X, Y = np.meshgrid(x_vals, y_vals)

        z = np.swapaxes((crosspower[2:Nkperp, 0 : Neta - 2]), 1, 0)

        if i == 2:
            mnn = 0.00000001
            mxx = 2

            nmm = colors.Normalize(vmin=mnn, vmax=mxx)
            cp = ax.pcolor(X, Y, z, cmap=cmocean.cm.balance, norm=nmm)
        else:
            mx = 1.0e14  # 1.0e14  # maximum power value to show in crosspower plot
            mn = 1.0e7
            nm = colors.LogNorm(vmin=mn, vmax=mx)
            cp = ax.pcolor(X, Y, z, cmap=cmocean.cm.balance, norm=nm)

        ideal_cax = colorbar(cp)
        ideal_cax.ax.tick_params(axis="both", which="major", labelsize=15)
        if i != 2:
            ideal_cax.ax.set_ylabel(r"P(k) mK$^2$ $h^{-3}$", size=16)
        # cbar = plt.colorbar(
        #     cp, label=r"P(k) mK$^2$ $h^{-3}$ Mpc$^3$", format="%.0e", cax=ax
        # )
        # cbar.set_ticks(np.logspace(np.log10(mn), np.log10(mx), setlen))
        # cbar.update_ticks()

        # ax.xticks(fontsize=16)
        # ax.yticks(fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=16)
        # ax.tick_params(axis='both', which='minor', labelsize=8)

        ax.set_title(titles[i], size=20)
        ax.set_xlabel(r"k$_\bot$ ($h$Mpc$^{-1}$)", fontsize=18)
        ax.set_ylabel(r"k$_\parallel$ ($h$Mpc$^{-1}$)", fontsize=18)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set(xlim=[0.008, 0.4])
        ax.set(ylim=[0.005, kpa[Neta - 2]])
    fig.tight_layout()
    plt.savefig(plotname, bbox_inches="tight")
    plt.show()
    return


def plot_xx_yy_ratio(filetag1, filetag2="multichips", prefix="alldata"):
    power2d_1, kper, kpa = get_2dpower(filetag1, filetag2, pol="xx")
    power2d_2, kper, kpa = get_2dpower(filetag1, filetag2, pol="yy")

    ratio2d = power2d_1 / power2d_2

    power_list2d = [power2d_1, power2d_2, ratio2d]

    label1 = f"{filetag1[:2]}_xx"
    label2 = f"{filetag1[:2]}_yy"
    titles = [
        label1,
        label2,
        f"Ratio ({label1}/{label2})",
    ]
    plot_ratio(
        power_list2d,
        kper,
        kpa,
        titles,
        plotname=f"{plotdir}/{label1}_vs_{label2}_{prefix}_ratio_power_2d.pdf",
    )

    return


if __name__ == "__main__":
    import itertools

    plotdir = "/astro/mwaeor/kchege/final_paper2_analysis/final_paper_plots"
    polas = ["yy", "xx"]

    ffff = [
        # "A1_patch_1000_1iono_1peel",
        # "A2_patch_4000_1iono_1peel",
        "B1_patch_1000_1iono_peel_1000",
        "C1_patch_1000_iono_1000_peel_1000",
        # "D1_SF_patch_1000_iono_1000_peel_1000",
        "B2_patch_4000_1iono_peel_4000",
        "C2_patch_4000_iono_4000_peel_4000",
        # "D2_SF_patch_4000_iono_4000_peel_4000",
        "B3_patch_4000_1iono_2sigma_peel",
        "C3_patch_4000_iono_2sigma_peel_4000",
        # "D3_SF_patch_4000_iono_2sigma_peel_4000",
    ]
    for ftag in ffff:
        plot_xx_yy_ratio(ftag)
        for t, typ in enumerate([1, 2, 3, 4]):
            iono_ftag = f"{ftag[:2]}_revised_type{typ}"
            plot_xx_yy_ratio(iono_ftag, prefix=f"itype{typ}")

    import sys

    sys.exit(0)

    filetags = list(itertools.combinations(ffff, 2))

    for couple in filetags:
        for pola in polas:
            FILETAG1_1, FILETAG1_2 = couple[0], "multichips"
            FILETAG2_1, FILETAG2_2 = couple[1], "multichips"

            label1 = f"{couple[0][:2]}_{pola}"
            label2 = f"{couple[1][:2]}_{pola}"

            power2d_1, kper, kpa = get_2dpower(FILETAG1_1, FILETAG1_2, pol=pola)
            power2d_2, kper, kpa = get_2dpower(FILETAG2_1, FILETAG2_2, pol=pola)

            ratio2d = power2d_1 / power2d_2

            power_list2d = [power2d_1, power2d_2, ratio2d]

            titles = [
                label1,
                label2,
                f"Ratio ({label1}/{label2})",
            ]

            plot_ratio(
                power_list2d,
                kper,
                kpa,
                titles,
                plotname=f"{plotdir}/{label1}_vs_{label2}_alldata_ratio_power_2d.pdf",
            )

    # typ = 3
    # for typ in [1, 2, 4]:
    #     itag = "D1_SF_patch_1000_iono_1000_peel_1000"
    #     couple = [f"D1_revised_type{typ}", f"D1_revised_type{typ}"]

    # couple = [f"D1_revised_type{typ}", f"D1_revised_type{typ}"]

    # for t, typ in enumerate([1, 2, 3, 4]):
    #     iono_filetags = [
    #         (f"{itag[0][:2]}_revised_type{typ}", "multichips") for itag in ffff
    #     ]
    #     # labels = [itag[0][:2] for itag in iono_filetags]

    #     iono_filetags = list(itertools.combinations(ffff, 2))

    #     for couple in iono_filetags:

    # couple = [
    #     # "D1_SF_patch_1000_iono_1000_peel_1000",
    #     # "D1_SF_patch_1000_iono_1000_peel_1000",
    #     "A1_patch_1000_1iono_1peel",
    #     "A2_patch_4000_1iono_1peel",
    # ]

    # fig, axlist = plt.subplots(1, 3, sharey=False, figsize=(18, 6), dpi=150)
    # for i, (ax, crosspower) in enumerate(zip(axlist, power_list2d)):
    #     x_vals = kper[2:Nkperp]
    #     y_vals = kpa[0 : Neta - 2]
    #     X, Y = np.meshgrid(x_vals, y_vals)

    #     # crosspower = ma.masked_where(crosspower <= 0, crosspower)

    #     z = np.swapaxes((crosspower[2:Nkperp, 0 : Neta - 2]), 1, 0)

    #     if i == 2:
    #         mnn = 0.00000001
    #         mxx = 2
    #         # mnn = -10 ** 2
    #         # mxx = 10 ** 2
    #         # nmm = colors.SymLogNorm(linthresh=0.03, vmin=mnn, vmax=mxx)
    #         nmm = colors.Normalize(vmin=mnn, vmax=mxx)
    #         cp = ax.pcolor(X, Y, z, cmap=cmocean.cm.balance, norm=nmm)
    #         #
    #         # cp = ax.contourf(
    #         #     X,
    #         #     Y,
    #         #     z,
    #         #     levels=np.logspace(np.log10(mnn), np.log10(mxx), 40, base=10),
    #         #     vmin=mnn,
    #         #     vmax=mxx,
    #         #     cmap=cmocean.cm.tarn,
    #         # )
    #     else:
    #         mx = 1.0e14  # 1.0e14  # maximum power value to show in crosspower plot
    #         mn = 1.0e7
    #         nm = colors.LogNorm(vmin=mn, vmax=mx)
    #         cp = ax.pcolor(X, Y, z, cmap=cmocean.cm.balance, norm=nm)

    #     # cp = ax.contourf(
    #     #     X,
    #     #     Y,
    #     #     z,
    #     #     levels=np.logspace(np.log10(mn), np.log10(mx), 40, base=10),
    #     #     locator=ticker.LogLocator(),
    #     #     vmin=mn,
    #     #     vmax=mx,
    #     #     cmap=cm.Spectral_r,
    #     # )

    #     ideal_cax = colorbar(cp)
    #     ideal_cax.ax.tick_params(axis="both", which="major", labelsize=15)
    #     if i != 2:
    #         ideal_cax.ax.set_ylabel(r"P(k) mK$^2$ $h^{-3}$", size=16)
    #     # cbar = plt.colorbar(
    #     #     cp, label=r"P(k) mK$^2$ $h^{-3}$ Mpc$^3$", format="%.0e", cax=ax
    #     # )
    #     # cbar.set_ticks(np.logspace(np.log10(mn), np.log10(mx), setlen))
    #     # cbar.update_ticks()

    #     # ax.xticks(fontsize=16)
    #     # ax.yticks(fontsize=16)
    #     ax.tick_params(axis="both", which="major", labelsize=16)
    #     # ax.tick_params(axis='both', which='minor', labelsize=8)

    #     ax.set_title(titles[i], size=20)
    #     ax.set_xlabel(r"k$_\bot$ ($h$Mpc$^{-1}$)", fontsize=18)
    #     ax.set_ylabel(r"k$_\parallel$ ($h$Mpc$^{-1}$)", fontsize=18)
    #     ax.set_xscale("log")
    #     ax.set_yscale("log")
    #     ax.set(xlim=[0.008, 0.4])
    #     ax.set(ylim=[0.005, kpa[Neta - 2]])
    # fig.tight_layout()
    # plt.savefig(
    #     f"/astro/mwaeor/kchege/final_paper2_analysis/final_paper_plots/{label1}_vs_{label2}_alldata_ratio_power_2d.pdf",
    #     bbox_inches="tight",
    # )
    # plt.show()
