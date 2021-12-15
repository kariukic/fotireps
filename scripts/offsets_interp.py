#!/usr/bin/env python3
import argparse
import json
import logging
import os
import sqlite3
import sys
import time

import numpy as np
import pandas as pd
import yaml
from fotireps.misc_utils import (
    get_sky_temp,
    iqr_bounds,
    nrmse,
    precise_dist,
    rmse,
    thermal_noise,
)
from fotireps.slurm_bash_scripts_writer import write_run_shift_command
from fotireps.plotting import plot_cthulhu_tec, plot_predictions
from scipy import spatial
from scipy.interpolate import CloughTocher2DInterpolator, RBFInterpolator, griddata
from yaml import SafeLoader as SafeLoader

__author__ = "Kariuki Chege"
__version__ = "0.1.0"
__date__ = "2021-09-27"


# TODO add second yaml

logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(module)s:%(levelname)s %(message)s",
)
log = logging.getLogger("OI")
logging_level = logging.INFO  # logging.DEBUG if args.debug else logging.INFO
log.setLevel(logging_level)
log.info("Karibu :-). \n This is OI from Fotireps %s-%s", __version__, __date__)


def loadfile(data_file):
    if data_file.split(".")[-1] == "csv":
        df = pd.read_csv(data_file)
    else:
        with open(data_file, "r") as f:
            unpacked = yaml.load(f, Loader=SafeLoader)

        (
            flux,
            ra,
            dec,
            name,
            dshift,
            rshift,
            ampscales,
            mudshift,
            murshift,
            dshift_dt,
            rshift_dt,
        ) = ([] for i in range(11))

        sourcelist = unpacked["sources"]

        for source in sourcelist:
            name.append(source)
            dec.append(unpacked["sources"][source]["dec"])
            ra.append(unpacked["sources"][source]["ra"])
            flux.append(unpacked["sources"][source]["flux_density"])
            # medians
            dshift.append(np.nanmedian(unpacked["sources"][source]["dec_shifts"]))
            rshift.append(np.nanmedian(unpacked["sources"][source]["ra_shifts"]))
            # means
            mudshift.append(np.nanmean(unpacked["sources"][source]["dec_shifts"]))
            murshift.append(np.nanmean(unpacked["sources"][source]["ra_shifts"]))
            # all timestamps offsets
            dshift_dt.append(unpacked["sources"][source]["dec_shifts"])
            rshift_dt.append(unpacked["sources"][source]["ra_shifts"])

            ampscales.append(np.nanmedian(unpacked["sources"][source]["amp_scales"]))

        df = pd.DataFrame(
            list(
                zip(
                    name,
                    ra,
                    dec,
                    flux,
                    dshift,
                    rshift,
                    ampscales,
                    mudshift,
                    murshift,
                    dshift_dt,
                    rshift_dt,
                )
            ),
            columns=[
                "name",
                "ra",
                "dec",
                "flux",
                "dshift",
                "rshift",
                "ampscale",
                "mudshift",
                "murshift",
                "dshift_dt",
                "rshift_dt",
            ],
        )

    return df


def get_cleanest_bright_sources(dif):
    z = len(dif)
    if "offsets" not in dif.keys():
        dif["offsets"] = precise_dist(
            np.zeros_like(dif.rshift.to_numpy()),
            np.zeros_like(dif.dshift.to_numpy()),
            dif.rshift.to_numpy(),
            dif.dshift.to_numpy(),
        )

    _, maxfence = iqr_bounds(dif.offsets.to_numpy())
    dif = dif[dif["offsets"] <= maxfence]

    x = len(dif)

    amps_minfence, amps_maxfence = iqr_bounds(dif.ampscale.to_numpy())

    dif = dif[(dif["ampscale"] >= amps_minfence) & (dif["ampscale"] <= amps_maxfence)]

    y = len(dif)

    log.info("total sources being outlier-checked: %s", z)
    log.info("offsets outliers: %s", z - x)
    log.info("ampscales outliers: %s", x - y)
    log.info("total outliers: %s", z - y)
    log.info("cleanest bright sources left: %s", y)

    return dif.name.to_list()


def label_dataframe(df, flux_threshold=None, boundaries=None):
    min_ra, max_ra, min_dec, max_dec = boundaries
    df["fov_bounds"] = np.where(
        (df["ra"] >= min_ra)
        & (df["ra"] <= max_ra)
        & (df["dec"] >= min_dec)
        & (df["dec"] <= max_dec),
        "in",
        "out",
    )

    df["flux_group"] = np.where(df["flux"] >= flux_threshold, "bright", "faint")

    log.info(
        "Total sources brighter than 2 sigma: %s", len(df[df["flux_group"] == "bright"])
    )
    log.info(
        "Total sources fainter than 2 sigma: %s", len(df[df["flux_group"] == "faint"])
    )

    bright_clean_sources = get_cleanest_bright_sources(
        df[(df["fov_bounds"] == "in") & (df["flux_group"] == "bright")]
    )

    offsets_quality = [
        "good" if source in bright_clean_sources else "bad"
        for source in df.name.to_list()
    ]

    df["offsets_quality"] = offsets_quality

    return df


def interpolate_offsets(
    ra2,
    dec2,
    rshift2,
    dshift2,
    boundaries,
    resolution=200,
    method=None,
):
    """Interpolate give RA and Dec offsets onto a regular grid.

    Parameters
    ----------
    ra2 : [type]
        [description]
    dec2 : [type]
        [description]
    rshift2 : [type]
        [description]
    dshift2 : [type]
        [description]
    boundaries : [type]
        [description]
    resolution : int, optional
        [description], by default 200
    method : [type], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """
    min_ra, max_ra, min_dec, max_dec = boundaries
    xgrid = np.mgrid[
        min_ra : max_ra : resolution * 1j, min_dec : max_dec : resolution * 1j
    ]
    xflat = xgrid.reshape(2, -1).T
    xobs = np.array([ra2, dec2]).T
    if method == "rbf":
        raflat = RBFInterpolator(
            xobs,
            rshift2,
            # kernel="cubic",
            # neighbors=1,
            # epsilon=1,
            # smoothing=1,
            # degree=1,
        )(xflat)
        ra_ygrid = raflat.reshape(resolution, resolution)

        decflat = RBFInterpolator(
            xobs,
            dshift2,
            # kernel="cubic",
            # neighbors=5,
            # epsilon=1,
            # smoothing=1,
            # degree=1,
        )(xflat)
        dec_ygrid = decflat.reshape(resolution, resolution)

    elif method == "clough":
        X = np.linspace(min_ra, max_ra, resolution)
        Y = np.linspace(min_dec, max_dec, resolution)
        X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation

        ra_interp = CloughTocher2DInterpolator(
            list(zip(ra2, dec2)), rshift2, fill_value=0
        )
        dec_interp = CloughTocher2DInterpolator(
            list(zip(ra2, dec2)), dshift2, fill_value=0
        )

        ra_ygrid = ra_interp(X, Y)
        dec_ygrid = dec_interp(X, Y)
    else:
        log.info("Interpolation method required.")
        exit()

    return ra_ygrid, dec_ygrid, xgrid, xflat


def predictions(ra, dec, ra_ygrid, dec_ygrid, xflat, resolution=200):
    """infer RA and Dec offsets from interpolated RA and Dec grids. \n Uses `spatial.KDTree.query`, a quick nearest neighbour lookup method obtained from `"https://stackoverflow.com/questions/53257607/get-closest-coordinate-in-2d-array"`

    Parameters
    ----------
    ra : [type]
        [description]
    dec : [type]
        [description]
    ra_ygrid : [type]
        [description]
    dec_ygrid : [type]
        [description]
    xflat : [type]
        [description]
    resolution : int, optional
        [description], by default 200

    Returns
    -------
    [type]
        [description]
    """

    ra_predecvals, dec_predecvals = [], []
    # ra_distances, dec_distances = [], []
    for xy in zip(ra, dec):
        _ra_distance, ra_index = spatial.KDTree(xflat).query(list(xy))
        _dec_distance, dec_index = spatial.KDTree(xflat).query(list(xy))
        # ra_distances.append(ra_distance)
        # dec_distances.append(dec_distance)
        ra_predecvals.append(ra_ygrid.reshape(resolution * resolution, 1)[ra_index])
        dec_predecvals.append(dec_ygrid.reshape(resolution * resolution, 1)[dec_index])

    return np.array(ra_predecvals)[:, 0], np.array(dec_predecvals)[:, 0]


def write_to_json(sources, ra_shifts, dec_shifts, name):
    """Writes RA and Dec offsets into a json file to be used in shifting RTS-style sourcelists

    Parameters
    ----------
    sources : [type]
        [description]
    ra_shifts : [type]
        [description]
    dec_shifts : [type]
        [description]
    name : [type]
        [description]
    """
    log.info("Total sources to write into json: %s", len(sources))
    shifts_dict = {}
    for source, rshift, dshift in zip(sources, ra_shifts, dec_shifts):
        shifts_dict[str(source)] = {"ra": rshift, "dec": dshift}

    log.info("Total sources written into json file: %s", len(shifts_dict))
    with open(name, "w", encoding="utf-8") as f:
        json.dump(shifts_dict, f, ensure_ascii=False, indent=4)
    return


def add_missing_sources(file1, file2):
    """Add missing sources in one rts peel sourcelist to the other. Needed because not all peel sources are includedd in cthulhu yaml file.

    Parameters
    ----------
    file1 : [type]
        [description]
    file2 : [type]
        [description]
    """

    with open(file1) as file1:
        alllines1 = file1.readlines()

    with open(file2) as file1:
        alllines2 = file1.readlines()

    lines1 = [l[:21] for l in alllines1 if l.startswith("SOURCE")]
    lines2 = [l[:21] for l in alllines2 if l.startswith("SOURCE")]

    missing = []
    for li in lines1:
        if li not in lines2:
            missing.append(li)

    log.info(f"{len(missing)} missing in sourcelist, included")

    with open(file2, "a") as the_file:
        for src in missing:
            for ll in alllines1:
                if ll.startswith(src):
                    ind = alllines1.index(ll)
                    while not alllines1[ind].startswith("ENDSOURCE"):
                        the_file.write(f"{alllines1[ind]}")
                        ind += 1
                    the_file.write("ENDSOURCE\n")
    return


def run_shifting_bash_script(
    srclist,
    metafits,
    shifts_json_filename,
    exclude_missing_sources,
):
    """Run a source shifting bash script `run_shift_command.sh`

    Parameters
    ----------
    srclist : [type]
        [description]
    metafits : [type]
        [description]
    shifts_json_filename : [type]
        [description]
    """
    peelsrclist = srclist.split("/")[-1].replace(
        "peel", "shifted_peel"
    )  # output name of the shifted peel sourcelist
    patchsrclist = srclist.split("/")[-1].replace(
        "peel", "shifted_patch"
    )  # output name of the shifted patch sourcelist
    write_run_shift_command()
    shcomm = f"sh ./run_shift_command.sh {srclist} {peelsrclist} {patchsrclist} {shifts_json_filename} {metafits} >> shifting_srclists.out"
    os.system(shcomm)
    time.sleep(60)
    if not exclude_missing_sources:
        add_missing_sources(srclist, peelsrclist)

    return


def shift_sources_without_interpolation(
    obsid,
    yamlfile=None,
    metafits=None,
    srclist=None,
    exclude_missing_sources=None,
):
    """Shifts sources with offsets directly from the yaml file. Only outlier offsets are changed to zeros

    Parameters
    ----------
    obsid : [type]
        [description]
    yamlfile : [type], optional
        [description], by default None
    metafits : [type], optional
        [description], by default None
    srclist : [type], optional
        [description], by default None
    """
    dif = loadfile(yamlfile)
    all_sources = dif.name.to_list()

    if "offsets" not in dif.keys():
        dif["offsets"] = precise_dist(
            np.zeros_like(dif.rshift.to_numpy()),
            np.zeros_like(dif.dshift.to_numpy()),
            dif.rshift.to_numpy(),
            dif.dshift.to_numpy(),
        )

    _, maxfence = iqr_bounds(dif.offsets.to_numpy())
    non_outlier_offsets_dif = dif[dif["offsets"] <= maxfence]

    amps_minfence, amps_maxfence = iqr_bounds(dif.ampscale.to_numpy())

    non_outlier_ampscales_dif = dif[
        (dif["ampscale"] >= amps_minfence) & (dif["ampscale"] <= amps_maxfence)
    ]

    all_good_sources = list(
        set(
            non_outlier_offsets_dif.name.to_list()
            + non_outlier_ampscales_dif.name.to_list()
        )
    )

    log.info("%s outliers found", len(all_sources) - len(all_good_sources))

    dif.set_index("name", inplace=True)

    final_shifts = [
        (
            dif.loc[source, "rshift"],
            dif.loc[source, "dshift"],
        )
        if source in all_good_sources
        else (0.0, 0.0)  # TODO Decide if to use the average offset
        for source in all_sources
    ]

    write_to_json(
        np.array(all_sources),
        [x[0] for x in final_shifts],
        [x[1] for x in final_shifts],
        f"{obsid}_offsets_no_interp.json",
    )

    if srclist:
        if metafits:
            run_shifting_bash_script(
                srclist,
                metafits,
                f"{obsid}_offsets_no_interp.json",
                exclude_missing_sources=exclude_missing_sources,
            )

    return


def save_to_db(df, ra_ygrid, dec_ygrid, boundaries, name):
    """Write all sources info, and the interpolation info onto an sqlite database.

    Parameters
    ----------
    df : [type]
        [description]
    ra_ygrid : [type]
        [description]
    dec_ygrid : [type]
        [description]
    boundaries : [type]
        [description]
    name : [type]
        [description]
    """
    index = [i for i in range(1, ra_ygrid.shape[0] + 1)]
    ragrid_df = pd.DataFrame(ra_ygrid, index=index, columns=index)
    decgrid_df = pd.DataFrame(dec_ygrid, index=index, columns=index)
    bounds_df = pd.Series(
        boundaries.T,
        index=["min_ra", "max_ra", "min_dec", "max_dec"],
    )

    conn = sqlite3.connect(name)

    names = ["source_data", "ra_grid", "dec_grid", "bounds"]
    for nm, dief in enumerate([df, ragrid_df, decgrid_df, bounds_df]):
        dief.to_sql(names[nm], conn, if_exists="replace", index=False)
    conn.close()
    return


def perfomance_stats(
    true_ra_shifts, true_dec_shifts, inferred_ra_shifts, inferred_dec_shifts
):
    """Quantify the control accuracy of the interpolation using RMSE and minmax normalized RMSE. Both values are in percentage.

    Parameters
    ----------
    true_ra_shifts : [type]
        [description]
    true_dec_shifts : [type]
        [description]
    inferred_ra_shifts : [type]
        [description]
    inferred_dec_shifts : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    rastats, decstats = [
        (rmse(*x) * 100, nrmse(*x) * 100)
        for x in [
            (true_ra_shifts, inferred_ra_shifts),
            (true_dec_shifts, inferred_dec_shifts),
        ]
    ]
    return (rastats[0], decstats[0]), (rastats[1], decstats[1])


def interp_and_shift_sources(
    obsid,
    method="rbf",
    n_max=None,
    radius=None,
    resolution=None,
    yamlfile=None,
    metafits=None,
    srclist=None,
    plot=True,
    pointing_center=None,
    shift_sourcelist=False,
    write_db=False,
    exclude_missing_sources=None,
    integration_time=None,
):
    """Interpolate offsets and output results.

    Parameters
    ----------
    obsid : [type]
        [description]
    method : str, optional
        [description], by default "rbf"
    n_max : [type], optional
        [description], by default None
    radius : [type], optional
        [description], by default None
    resolution : [type], optional
        [description], by default None
    yamlfile : [type], optional
        [description], by default None
    metafits : [type], optional
        [description], by default None
    srclist : [type], optional
        [description], by default None
    plot : bool, optional
        [description], by default True
    pointing_center : tuple, optional
        [description], by default (0.0, -27.0)
    shift_sourcelist : bool, optional
        [description], by default False
    write_db : bool, optional
        [description], by default False
    """
    df = loadfile(yamlfile)

    all_sources = df.name.to_list()
    log.info(f"Total sources in input dataset: {len(all_sources)}")
    log.info(
        f"All sources flux density range, (min, max), Jy: {(np.min(df.flux), np.max(df.flux))}"
    )

    boundaries = np.array(
        [
            ((-1 * radius) + pointing_center[0]),
            radius + pointing_center[0],
            ((-1 * radius) + pointing_center[1]),
            (radius + pointing_center[1]),
        ]
    )

    log.info(
        f"Interpolation grid size: {radius}deg by {radius}deg \n Resolution: {resolution}, Grid centre: {pointing_center} \n Grid bounds (min_ra, max_ra, min_dec, max_dec): {boundaries}"
    )

    tisis = get_sky_temp(obsid) + 50  # sky temp + 50k MWA receiver temp
    log.info("Obsid Tsys: %s", tisis)
    two_sigma_thermal_noise_level = (
        thermal_noise(Tsys=tisis, integration_time=integration_time) * 2
    )

    df1 = label_dataframe(
        df, flux_threshold=two_sigma_thermal_noise_level, boundaries=boundaries
    )

    within_bounds_df = df1.loc[lambda dd: dd["fov_bounds"] == "in", :]

    log.info("Sources within interpolation boundaries: %s", len(within_bounds_df))
    log.info(
        "Sources outside interpolation boundaries: %s",
        len(df1.loc[lambda dd: dd["fov_bounds"] == "out", :]),
    )

    df2 = within_bounds_df.loc[
        (within_bounds_df["flux_group"] == "bright")
        & (within_bounds_df["offsets_quality"] == "good")
    ]

    log.info("Total interpolants to be used: %s", len(df2))
    log.info(
        "Interpolants flux density range, (min, max), Jy: (%s, %s)",
        np.min(df2.flux),
        np.max(df2.flux),
    )

    ra_ygrid, dec_ygrid, xgrid, xflat, = interpolate_offsets(
        df2.ra.to_numpy(),
        df2.dec.to_numpy(),
        df2.rshift.to_numpy(),
        df2.dshift.to_numpy(),
        boundaries,
        resolution=resolution,
        method=method,
    )

    inferred_ra_offsets, inferred_dec_offsets = predictions(
        within_bounds_df.ra,
        within_bounds_df.dec,
        ra_ygrid,
        dec_ygrid,
        xflat,
        resolution=resolution,
    )

    offset_predictions = [
        (
            inferred_ra_offsets[within_bounds_df["name"].to_list().index(sourcename)],
            inferred_dec_offsets[within_bounds_df["name"].to_list().index(sourcename)],
        )
        if sourcename in within_bounds_df["name"].to_list()
        else (np.nan, np.nan)
        for sourcename in df1["name"].to_list()
    ]

    df1["inferred_ra_offsets"] = [x[0] for x in offset_predictions]
    df1["inferred_dec_offsets"] = [x[1] for x in offset_predictions]

    log.info("Writing final combination of offsets based on source params ...")
    # reindexing for easier final shifts reading
    # inplace=False makes sure we do not modify df1 itself
    reind_df1 = df1.set_index("name", inplace=False)

    final_shifts = [
        (
            reind_df1.loc[source, "rshift"],
            reind_df1.loc[source, "dshift"],
        )
        if source in df2.name.to_list()
        else (
            reind_df1.loc[source, "rshift"],
            reind_df1.loc[source, "dshift"],
        )
        if reind_df1.loc[source, "fov_bounds"] == "out"
        else (
            reind_df1.loc[source, "inferred_ra_offsets"],
            reind_df1.loc[source, "inferred_dec_offsets"],
        )
        for source in all_sources
    ]

    df1["final_ra_offsets"] = [x[0] for x in final_shifts]
    df1["final_dec_offsets"] = [x[1] for x in final_shifts]

    log.info("Interpolation Done")

    control_df = df1[df1.name.isin(df2.name.to_list())]
    radec_rmse, radec_nrmse = perfomance_stats(
        control_df["rshift"],
        control_df["dshift"],
        control_df["inferred_ra_offsets"],
        control_df["inferred_dec_offsets"],
    )
    log.info(
        "Control sources control interpolation RMSE (ra, dec): %s",
        radec_rmse,
    )
    log.info("Control sources interpolation normalized RMSE (ra, dec): %s", radec_nrmse)

    if radec_rmse[0] < 1000:
        log.info("Interpolation converged.")
    else:
        log.info("Interpolation failed to converge.")

    name = f"{obsid}_2sigmainterp_{len(df1)}total_{resolution}res"

    log.info("Outputs basename %s", name)

    if plot:
        log.info("Plotting...")
        plot_predictions(
            within_bounds_df.ra,
            within_bounds_df.dec,
            control_df["inferred_ra_offsets"],
            control_df["inferred_dec_offsets"],
            df2.ra.to_numpy(),
            df2.dec.to_numpy(),
            xgrid,
            ra_ygrid,
            dec_ygrid,
            df2.rshift.to_numpy(),
            df2.dshift.to_numpy(),
            df2.flux.to_numpy(),
            f"{name}.pdf",
        )

        plot_cthulhu_tec(
            df1.dropna().ra,
            df1.dropna().dec,
            df1.dropna()["inferred_ra_offsets"],
            df1.dropna()["inferred_dec_offsets"],
            radius=radius,
            name=f"{name}_g2stec.pdf",
        )

    if write_db:
        log.info("Writing outputs to sqlite database...")
        try:
            save_to_db(df1, ra_ygrid, dec_ygrid, boundaries, f"{name}.sqlite")
        except:
            print("something bad happened couldn't save interpolation info to database")

    if shift_sourcelist:
        log.info("Writing json file for rts sourcelist shifting ...")
        write_to_json(
            all_sources,
            df1["final_ra_offsets"],
            df1["final_dec_offsets"],
            f"{name}.json",
        )

        log.info("Shifting rts sourcelist ...")
        run_shifting_bash_script(
            srclist, metafits, f"{name}.json", exclude_missing_sources
        )
    log.info("Kwaheri :-)")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ideal visibilities Simulation Set Up")
    group1 = parser.add_argument_group("general")

    group1.add_argument(
        "--obsid",
        type=str,
        default="",
        help="obsid",
    )
    group1.add_argument(
        "--method",
        choices=["rbf", "others"],
        default="rbf",
        type=str,
        required=False,
        help="interpolation method",
    )
    group1.add_argument("--metafits", required=False, help="obsid metafits")
    group1.add_argument("--srclist", required=False, help="Peel sourcelist")
    group1.add_argument(
        "--yamlfile",
        required=True,
        help="yamlfile with offsets for n1 and n2 number of sources. This is the one used to do the interpolation.",
    )
    group1.add_argument(
        "--yamlfile2",
        required=False,
        help="This yamlfile should just provides the ra and decs for the n3/n4 source offsets to be predicted",
    )
    group1.add_argument(
        "--output_directory",
        required=False,
        help="directory to save stuff in",
    )

    group1.add_argument(
        "--resolution",
        "-r",
        required=False,
        type=int,
        default=250,
        help="Resolution",
    )

    group1.add_argument(
        "--n_max",
        action="store",
        default=1000,
        required=False,
        type=int,
        help="nsources to interpolate offsets of",
    )
    group1.add_argument(
        "--radius",
        action="store",
        default=20,
        required=False,
        type=int,
        help="Radius/length of interpolation area",
    )

    group1.add_argument(
        "--integration_time",
        action="store",
        default=8,
        required=False,
        type=int,
        help="Integration time for computing thermal noise and sources SNR",
    )

    group1.add_argument(
        "--interp",
        "-i",
        action="store_true",
        help="Perform interpolation",
    )

    group1.add_argument(
        "--plot",
        "-p",
        default=True,
        required=False,
        action="store_true",
        help="Make performance plots",
    )

    group1.add_argument(
        "--write_db",
        "-l",
        default=False,
        required=False,
        action="store_true",
        help="Write all interpolation outputs and source information into a database",
    )

    group1.add_argument(
        "--shift_sourcelist",
        "-m",
        default=False,
        required=False,
        action="store_true",
        help="Output shifted RTS peel and patch sourcelists",
    )
    group1.add_argument(
        "--exclude_missing_sources",
        dest="exclude_missing_sources",
        action="store_true",
        help="The provided sourcelist is ready and does not need srclist_by_beam to be run on it. [default: false]",
    )
    group1.add_argument(
        "--pointing_center",
        default=(0.0, -27.0),
        required=False,
        action="store_true",
        help="Output shifted RTS peel and patch sourcelists",
    )

    args = parser.parse_args()

    assert args.obsid is not None

    if args.shift_sourcelist and not args.srclist:
        logging.error("sourcelist needed for 'shift_sourcelist' option ")
        sys.exit()

    if args.srclist:
        try:
            assert os.path.isfile(args.srclist)
        except AssertionError() as error:
            raise AssertionError(
                f"{args.srclist} seems not to be a valid sourcelist path."
            ) from error

    output_directory = (
        args.output_directory if args.output_directory else os.path.abspath(".")
    )

    # Get into the specified output/working directory
    if args.output_directory and not os.path.exists(output_directory):
        os.system(f"mkdir -p {output_directory}")

    # make sure we actually have output_directory
    try:
        assert os.path.isdir(output_directory)
        os.chdir(output_directory)
    except Exception as error:
        raise AssertionError(
            f"{output_directory} could not be found or made for some unknown reason."
        ) from error

    if args.interp:
        assert args.method is not None

        if not args.metafits:
            try:
                log.info("Trying to download a metafits file for the obsid")
                import urllib.request

                metafits = f"{args.obsid}.metafits"
                urllib.request.urlretrieve(
                    f"http://ws.mwatelescope.org/metadata/fits?obs_id={args.obsid}",
                    metafits,
                )
            except Exception:
                log.error("Metafits file is needed and download attempt failed.")
                sys.exit()
        else:
            metafits = args.metafits

        interp_and_shift_sources(
            args.obsid,
            method=args.method,
            n_max=args.n_max,
            radius=args.radius,
            resolution=args.resolution,
            pointing_center=args.pointing_center,
            yamlfile=args.yamlfile,
            metafits=metafits,
            srclist=args.srclist,
            plot=args.plot,
            shift_sourcelist=args.shift_sourcelist,
            write_db=args.write_db,
            exclude_missing_sources=args.exclude_missing_sources,
            integration_time=args.integration_time,
        )

    else:
        log.info("No interpolation being done!")
        shift_sources_without_interpolation(
            args.obsid,
            yamlfile=args.yamlfile,
            metafits=args.metafits,
            srclist=args.srclist,
            exclude_missing_sources=args.exclude_missing_sources,
        )
