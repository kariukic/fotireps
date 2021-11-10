#!/usr/bin/env python3

"""Friends of the Ionosphere, RTS EoR Processing Suite"""
import sys
import os

import logging
import numpy as np
from argparse import ArgumentParser

from fotireps.slurm_bash_scripts_writer import CookScripts, write_tiles_to_flag_file

__author__ = "Kariuki Chege"
__version__ = "0.1"
__date__ = "2021-10-21"

# TODO find a better way to do this
VIRTUAL_ENV = "source /astro/mwaeor/kchege/virtenvs/karis/bin/activate"


def main():
    """40REPS"""
    parser = ArgumentParser(
        "40REPS",
        description="""an RTS EoR Processing Suite, for Friends Of The Ionosphere, 
        meant to save us from near 40 repetitions of trying to run all the steps.""",
    )
    group1 = parser.add_argument_group("Configuration options")
    group1.add_argument(
        "--obsid",
        required=True,
        help="obsid to run jobs on",
    )
    group1.add_argument("--metafits", type=str, required=False, help="obsid metafits")
    # group1.add_argument(
    #     "--job",
    #     nargs="+",
    #     required=True,
    #     # choices=[
    #     #     "download",
    #     #     "cotter",
    #     #     "rts_setup",
    #     #     "rts_run",
    #     #     "cthulhu",
    #     #     "chips",
    #     #     "all",
    #     # ],
    #     help="job to run",
    # )
    group1.add_argument(
        "--rts",
        dest="rts",
        action="store_true",
        help="rts mode. [default: false]",
    )
    group1.add_argument(
        "--download",
        dest="download",
        action="store_true",
        help="download the obsid data from mwa asvo. [default: false]",
    )
    group1.add_argument(
        "--cotter",
        dest="cotter",
        action="store_true",
        help="Run cotter on gpu box files. [default: false]",
    )

    group1.add_argument(
        "--cthulhu",
        dest="cthulhu",
        action="store_true",
        help="Run cthulhu as well mode. [default: false]",
    )
    group1.add_argument(
        "--offset_offsets",
        dest="offset_offsets",
        action="store_true",
        help="Run offsets analysis. [default: false]",
    )
    group1.add_argument(
        "--chips",
        dest="chips",
        action="store_true",
        help="Run chips. For individual obsid or multiple obsids [default: false]",
    )

    group1.add_argument(
        "--fire",
        action="store_true",
        help="Run all the jobs in the right order on garrawarla",
    )

    group1.add_argument(
        "--data_format",
        choices=["uvfits", "gpu_boxes"],
        default="gpu_boxes",
        type=str,
        required=False,
        help="input data format",
    )
    group1.add_argument(
        "--boxes_path",
        type=str,
        required=False,
        help="Path to Directory containing GPU box files for the obsid, must be provided for rts jobs",
    )

    group2 = parser.add_argument_group("RTS calibration options")

    group2.add_argument(
        "--no_srclist_by_beam",
        dest="no_srclist_by_beam",
        action="store_true",
        help="The provided sourcelist is ready and does not need srclist_by_beam to be run on it. [default: false]",
    )

    group2.add_argument(
        "--patch",
        default=0,  # 1000
        type=int,
        help="Run DI with this number of sky model sources (patch)",
    )
    group2.add_argument(
        "--patch_time_config",
        type=int,
        nargs=2,
        default=None,  # (32, 2),
        help="""Change the default time calibration parameters for the patch only. The needed params are (integration_time, CorrDumpTime) e.g 32 2. If you just want to run defaults ignore this.
        The order is important.
        """,
    )

    group2.add_argument(
        "--peel",
        type=int,
        nargs=3,
        default=None,  # (0, 1000, 1000),
        help="""Run DD calibration with this number of sky model sources (num_full_dd_cals, numiono,  numpeel)
        The order is important.
        """,
    )

    group2.add_argument(
        "--fov_cutoff",
        type=float,
        default=30.0,
        required=False,
        help="Radius for srclist cutoff",
    )

    group2.add_argument(
        "--corrdumptime",
        action="store",
        default=2,
        required=False,
        type=int,
        help="Correlator time resolution",
    )
    group2.add_argument(
        "--integration_time",
        action="store",
        default=8,
        required=False,
        type=int,
        help="Integration time",
    )
    group2.add_argument(
        "--subbands",
        nargs="+",
        required=False,
        default=[int(sb) for sb in np.arange(1, 25)],
        help="subbands, can be a list of values in 1 to 24",
    )
    group2.add_argument(
        "--flag_tiles",
        nargs="+",
        required=False,
        default=None,
        help="IDs of tiles to flag",
    )
    group2.add_argument(
        "--srclist",
        required=False,
        default="/astro/mwaeor/kchege/srclists/srclist_pumav3_EoR0aegean_fixedEoR1pietro+ForA_phase1+2.txt",
        help="model/sourcelist path",
    )

    group2.add_argument(
        "--di_gains_path",
        required=False,
        default=None,
        help="path to DI (patch) gains solution files with .dat extensions. Needed if you are not running patch in the same run as peel",
    )
    group2.add_argument(
        "--dd_logs_path",
        required=False,
        default=None,
        help="path to DD log files with .log extensions. Needed if you are not running peel in the same run as cthulhu",
    )

    group3 = parser.add_argument_group("Cthulhu options")
    group3.add_argument(
        "--radius",
        type=float,
        default=25.0,
        help="Radius for cthulhu iono_radius",
    )

    group6 = parser.add_argument_group("CHIPS options")
    group6.add_argument(
        "--uvfits_path",
        type=str,
        required=False,
        help="Path to Directory containing visibilities uvfits files of which to calculate the power spectrum.",
    )
    group6.add_argument(
        "--band",
        choices=["low", "high"],
        default="low",
        type=str,
        required=False,
        help="low or high frequency band",
    )
    group6.add_argument(
        "--eorfield",
        choices=[0, 1],
        default=0,
        type=int,
        required=False,
        help="low or high frequency band",
    )

    group2.add_argument(
        "--chipstags",
        type=str,
        nargs=2,
        default=None,
        help="Two UNIQUE tags to name the output chips files. I repeat, UNIQUE. Otherwise, your output files might overwrite someone else's chips output of the same obsid with the same name. Sorry.",
    )

    group4 = parser.add_argument_group("Output options")
    group4.add_argument(
        "--output_directory",
        required=False,
        help="directory to save stuff in",
    )

    group5 = parser.add_argument_group("Extra options")
    group5.add_argument(
        "--email",
        required=False,
        action="store_true",
        help="Notify you via email if your jobs fail",
    )

    args = parser.parse_args()

    # logging configuration
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format="%(module)s:%(levelname)s %(message)s",
    )
    log = logging.getLogger("")
    # logging_level = logging.DEBUG if args.debug else logging.INFO
    logging_level = logging.INFO
    log.setLevel(logging_level)
    log.info("Karibu. This is FotIREPS %s-(%s)", __version__, __date__)

    output_directory = (
        args.output_directory if args.output_directory else os.path.abspath(".")
    )
    # set some paths straight
    if args.di_gains_path:
        di_gains_path = os.path.abspath(args.di_gains_path)
    if args.dd_logs_path:
        dd_logs_path = os.path.abspath(args.dd_logs_path)

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

    # Initiate the bash scripts writer object
    scripts_writer = CookScripts(
        obsid=args.obsid,
        sourcelist=args.srclist,
        boxes_path=args.boxes_path,
        virtual_env=VIRTUAL_ENV,
    )

    # Remove jobs from this list as you implement them
    jobs_not_implemented_yet = ["download", "cotter", "offset_offsets"]

    if args.download or args.cotter or args.offset_offsets:
        logging.warning(
            "You have asked for one or more of these unimplemented jobs; %s, they won't run.",
            jobs_not_implemented_yet,
        )
    # Add this up as you implement all the other jobs
    implemented = [args.rts, args.cthulhu, args.chips]
    if not any(implemented):
        logging.error(
            "No jobs left to run. Maybe all jobs you asked for are unimplemented."
        )
        sys.exit()

    rts_jobs_available = False
    if args.rts:
        rtsJobs = [args.patch, args.peel]
        if any(rtsJobs):  # if there is at least one rts jobs, proceed to preocess it.
            rts_jobs_available = True
            # make sure we actually have the supposed gpu boxes directory
            try:
                assert os.path.isdir(args.boxes_path)
            except Exception as error:
                raise AssertionError(
                    f"{args.boxes_path} is not a valid directory. Exiting!"
                ) from error
            # Next check if we have a metafits file, try looking for it in some possible locations. if not found we can't run the rts.
            if not args.metafits:
                metafitsfile = f"{args.obsid}.metafits"
                possiblepath1 = f"{args.boxes_path}/{metafitsfile}"
                possiblepath2 = f"{os.path.abspath('.')}/{metafitsfile}"
                possiblepath3 = f"{args.output_directory}/{metafitsfile}"
                for ppath in [possiblepath1, possiblepath2, possiblepath3]:
                    if os.path.isfile(ppath):
                        scripts_writer.metafits = os.path.abspath(ppath)
                        logging.info(
                            "Found metafits %s. Make sure its the correct one or provide the right one!",
                            ppath,
                        )
                        break

                # IF we do not find a metafits file, ask for it and exit
                if not scripts_writer.metafits:
                    logging.error("No metafits file provided or found! Exiting")
                    sys.exit()

            if args.patch:
                scripts_writer.run_patch = True
                DD_calibrators_numbers = args.peel
            if args.peel:
                scripts_writer.run_peel = True
                if not args.patch:
                    # ask for the directory with patch outputs to use
                    # if not in output directory, copy them to here.
                    # run the peel
                    try:
                        assert os.path.isdir(di_gains_path)
                    except Exception as error:
                        raise AssertionError(
                            "If you are not running the patch, a valid path to a directory containing DI (patch) solutions is needed."
                        ) from error
                    # TODO check if we can do it without copying to save the space
                    os.system(f"cp {di_gains_path}/*.dat {os.path.abspath('.')}")

                if args.no_srclist_by_beam:
                    with open(args.srclist) as sourcelist:
                        all_lines = sourcelist.readlines()
                        total_sources = len(
                            [l for l in all_lines if l.startswith("ENDSOURCE")]
                        )
                    if total_sources < np.max(DD_calibrators_numbers):
                        log.info(
                            "The sourcelist has less source (%s) than one of the provided peel params sources (%s). Setting %s as the new peel number.",
                            total_sources,
                            *args.peel,
                            total_sources,
                        )

                        DD_calibrators_numbers = list(
                            map(
                                lambda x: total_sources if x > total_sources else x,
                                DD_calibrators_numbers,
                            )
                        )

            scripts_writer.rts_setup(
                DI_calibrators_numbers=args.patch,
                DD_calibrators_numbers=DD_calibrators_numbers,
                integration_time=args.integration_time,
                CorrDumpTime=args.corrdumptime,
                cutoff=args.fov_cutoff,
                subbands=args.subbands,
                patch_time_config=args.patch_time_config,
                no_srclist_by_beam=args.no_srclist_by_beam,
            )
            # If needed, write the 'flagged_tiles.txt' file needed by the rts to flag
            if args.flag_tiles:
                logging.info("Writing tiles %s to flagged_tiles.txt", args.flag_tiles)
                write_tiles_to_flag_file(args.flag_tiles)

            scripts_writer.rts_run()
        else:
            logging.info(
                "You requested for rts mode without specifying any peel or patch (or both) task to do. Will assume no rts work to be done"
            )
            if (
                len(implemented) == 1
            ):  # means that --rts was the only specified job argument
                logging.info("No jobs left to run.")
                sys.exit()
            else:
                rts_jobs_available = False

    if args.cthulhu:
        if not args.peel:
            # ask for the directory with peel log files for cthulhu to use
            try:
                assert os.path.isdir(dd_logs_path)
            except Exception as error:
                raise AssertionError(
                    "If you are not running peel in the same run as cthulhu, a valid path to a directory containing peel log files is needed."
                ) from error
            # TODO check if we can do it without copying to save the space
            # os.system(f"cp {dd_logs_path}/*.dat {os.path.abspath('.')}")
        else:
            dd_logs_path = os.path.abspath(".")

        scripts_writer.cthulhu(logs_dir=dd_logs_path)

    if args.chips:
        # We need a path to uvfits file if we are not running DI or DD in the same run
        if not args.peel and not args.patch:
            if not args.uvfits_path:
                logging.error("Uvfits files needed for chips to work.")
                sys.exit()
            else:
                # make sure we actually have the supposed uvfits directory
                try:
                    assert os.path.isdir(args.uvfits_path)
                except Exception as error:
                    raise AssertionError(
                        f"{args.uvfits_path} seems not to be a valid directory with the needed chips uvfits files."
                    ) from error
            uvfits_path = os.path.abspath("args.uvfits_path")
        else:
            # if we are running the patch or peel in this same run, then we should have uvfits file in this directory
            uvfits_path = os.path.abspath(".")

        # here we just want to add a prefix to the expected uvfits file becaues patch uvfits are renamed by default.
        uvfits_prefix = f"{args.obsid}_patch" if args.patch and not args.peel else None

        # make sure chips tags are provided
        if not args.chipstags:
            logging.error(
                "For chips, it's crucial that you provide 2 unique tags to label the output chips power spectrum files. eg., '--chipstags unique_tag1 unique_tag2'"
            )
            sys.exit()

        scripts_writer.chips(
            uvfits_path=uvfits_path,
            tags=args.chipstags,
            band=args.band,
            field=args.eorfield,
            uvfits_prefix=uvfits_prefix,
        )

    # Finally write a scheduled bash script to run all the requested jobs in sequence and according to their interdependency

    jobs = [
        job
        for job, requested in zip(["rts", "cthulhu", "chips"], implemented)
        if requested
    ]
    if rts_jobs_available:
        jobs = ["rts_setup", "rts_run", *jobs]
        jobs.remove("rts")
    else:
        if "rts" in jobs:
            jobs.remove("rts")

    scripts_writer.run_jobs(jobs)

    if args.fire:
        os.system("sbatch run_jobs.sh")

    logging.info("All done. FOTIREPS scripts writer signing out, Kwaheri(-:")


if __name__ == "__main__":
    main()
