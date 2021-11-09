import os

import cmasher as cmr
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from cthulhu.plot_tools import plot_tec, plot_vector_arrows, setup_subplot
from cthulhu.reconstruct import Obsid
from mpl_toolkits.axes_grid1 import make_axes_locatable

sequential_colors = cmr.get_sub_cmap("cmr.lilac", 0.2, 0.8)

sns.set()
sns.set_theme(context="paper", palette="bright", style="whitegrid", font_scale=1.5)
# sns.set_style("whitegrid")

# Take 5 colors from rainforest in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    "cmr.rainforest", 4, cmap_range=(0.15, 0.85), return_fmt="hex"
)


def colorbar(mappable):
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


def plot_predictions(
    ra1,
    dec1,
    ra_predecvals,
    dec_predecvals,
    ra2,
    dec2,
    xgrid,
    ra_ygrid,
    dec_ygrid,
    rshift2,
    dshift2,
    dfluxes2,
    name=None,
):
    ra_ygrid *= 60
    dec_ygrid *= 60
    rshift2 *= 60
    dshift2 *= 60
    ra_predecvals = np.array(ra_predecvals) * 60
    dec_predecvals = np.array(dec_predecvals) * 60

    xobs2 = np.array([ra2, dec2]).T
    fig, ax = plt.subplots(2, 2, figsize=(14, 10), dpi=250)
    ax1, ax2, ax3, ax4 = (
        ax[0, 0],
        ax[0, 1],
        ax[1, 0],
        ax[1, 1],
        # ax[2, 0],
        # ax[2, 1],
    )

    ax1.pcolormesh(*xgrid, ra_ygrid, shading="gouraud")
    rap = ax1.scatter(*xobs2.T, c=rshift2, s=dfluxes2 ** 2, ec="g", facecolors="none")

    _rapred = ax1.scatter(
        ra1,
        dec1,
        s=7.0,  # df.flux ** 1.5,
        # ec="k",
        marker="*",
        linewidth=0,
    )
    cbar = plt.colorbar(rap, ax=ax1)
    cbar.set_label("Offsets [arcmins]")

    ax2.pcolormesh(*xgrid, dec_ygrid, shading="gouraud")  # shading="gouraud"
    decp = ax2.scatter(*xobs2.T, c=dshift2, s=dfluxes2 ** 2, ec="g", facecolors="none")
    _decpred = ax2.scatter(
        ra1,
        dec1,
        s=7.0,
        # ec="k",
        marker="*",
        linewidth=0,
    )
    cbar = plt.colorbar(decp, ax=ax2)
    cbar.set_label("Offsets [arcmins]")

    ax1.set_xlabel("Right Ascension [deg]")
    ax1.set_ylabel("Declination [deg]")
    ax2.set_xlabel("Right Ascension [deg]")
    ax2.set_ylabel("Declination [deg]")
    # plt.setp(ax2.get_yticklabels(), visible=False)

    ax3.plot(
        rshift2,  # df.true_ra_shift,
        rshift2,  # df.true_ra_shift,
        "r",
        # markersize="3",
        lw=1,
        linestyle="--",
        label="Fiducial RA offsets",
        zorder=3,
    )

    _bb = sns.scatterplot(
        x=rshift2,
        y=ra_predecvals,
        hue=dfluxes2,
        size=dfluxes2,
        sizes=(10, 200),
        # hue_norm=(0, 7),
        palette=sequential_colors,
        # facecolors="none",
        legend=False,
        label="Inferred RA offsets",
        marker="p",
        zorder=2.5,
        linewidth=0,
        ax=ax3,
    )

    ax4.plot(
        dshift2,
        dshift2,
        "r",
        # markersize="3",
        linestyle="--",
        lw=1,
        label="Fiducial Dec offsets",
        zorder=3,
    )

    _dd = sns.scatterplot(
        x=dshift2,  # df.true_dec_shift,
        y=dec_predecvals,
        hue=dfluxes2,  # df.flux,
        size=dfluxes2,  # df.flux,
        sizes=(10, 200),
        linewidth=0,
        # hue_norm=(0, 7),
        palette=sequential_colors,
        # facecolors="none",
        marker="p",
        legend=False,
        label="Inferred Dec offsets",
        zorder=2.5,
        ax=ax4,
    )

    ax3.set_xlabel("Fiducial RA offsets [arcmins]")
    ax3.set_ylabel("Inferred RA offsets [arcmins]")
    ax4.set_xlabel("Fiducial Dec offsets [arcmins]")
    ax4.set_ylabel("Inferred Dec offsets [arcmins]")

    ax3.axhline(y=0, xmin=0, xmax=1, color="k", linestyle="--", zorder=10)
    ax3.plot(
        rshift2,
        ra_predecvals - rshift2,
        marker="X",
        markersize="4",
        linestyle="None",
        color=colors[2],
        alpha=0.6,
        zorder=1,
        label="Residuals",
    )

    ax3.legend(loc="best")

    ax4.axhline(y=0, xmin=0, xmax=1, color="k", linestyle="--", zorder=10)
    ax4.plot(
        dshift2,
        dec_predecvals - dshift2,
        marker="X",
        markersize="4",
        linestyle="None",
        color=colors[2],
        alpha=0.6,
        zorder=1,
        label="Residuals",
    )
    ax4.legend(loc="best")

    fig.tight_layout()
    if name:
        plt.savefig(name)
    else:
        plt.show()
    return


def plot_cthulhu_tec(
    ra,
    dec,
    ra_shifts,
    dec_shifts,
    radius=20,
    frequency=154.235,
    name="g2s_tec_screen.png",
):
    _fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    o = Obsid((ra, dec, ra_shifts, dec_shifts), radius=radius, frequency=frequency)
    o.reconstruct_tec(resolution=400)
    o.obsid_metric()
    setup_subplot(
        axis=ax,
        title=f"""TEC  Med: {round(o.metrics[0][0], 4)}, pca: {round(o.metrics[1][0], 4)} Metric: {round(o.metric, 4)}""",
    )

    plot_vector_arrows(axis=ax, obsid=o, norm=True)
    tec = plot_tec(axis=ax, obsid=o, colourbar=False)  # vlim=(0, 2)

    label = "TEC [TECU]"
    cb = colorbar(tec)
    cb.set_label(label, fontsize=14, labelpad=15)
    plt.savefig(name, dpi=300)
    return


if __name__ == "__main__":
    pass
