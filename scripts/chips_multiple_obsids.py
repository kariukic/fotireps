"""Run chips on multiple obsids"""
import os
import sys
import logging
import glob
from argparse import ArgumentParser
from fotireps.chips_runners import (
    write_run_clean_intermed,
    write_run_all_astro_script,
    write_make_obs_astro_script,
)
from fotireps.slurm_bash_scripts_writer import add_permissions

##steps
# get the directory with the needed uvfits for each obsid
# add the "tag1/tag2" identifier for each of the obsids into a chips input text file
# tag1 should be the same for all, say the iteration run like A1_1000_patch_only
# tag2 should be the unique obsid identifier

# for each obsid, make a subdirectory in /astro/mwaeor/MWA/data/{tag1} named after tag2.
# run ln sf linking each obsid uvfits to its respective subdirectory


def chips_config(obsids, uvfits_dirs, runtag, patchrun=False):

    input_file = f"{runtag}_chips_inputfile.txt"

    with open(input_file, "w") as text_file:
        for _i, (obsid, obsid_uvfits_directory) in enumerate(zip(obsids, uvfits_dirs)):

            unique_obsid_directory = f"{runtag}/{obsid}"
            # for each obsid, make a directory in /astro/mwaeor/MWA/data/
            comm1 = f"mkdir -p /astro/mwaeor/MWA/data/{unique_obsid_directory}"
            os.system(comm1)

            # If uvfits are from a peel run they have the uvdump_*.uvfits name expected by chips
            # However if they are from a patch, chances are they were rename to have an 'obsid_patch' prefix.
            # we have to remove it in the symlinks read by chips
            if not patchrun:
                comm2 = f"ln -sf {obsid_uvfits_directory}/uvdump_*.uvfits /astro/mwaeor/MWA/data/{unique_obsid_directory}"
                os.system(comm2)

            else:
                all_patch_uvfits = glob.glob1(
                    obsid_uvfits_directory, f"{obsid}_patch_uvdump*.uvfits"
                )
                for fits in all_patch_uvfits:
                    fits_symlink_name = fits[len(f"{obsid}_patch_") :]

                    comm2 = f"ln -sf {obsid_uvfits_directory}/{fits} /astro/mwaeor/MWA/data/{unique_obsid_directory}/{fits_symlink_name}"
                    os.system(comm2)

            text_file.write(unique_obsid_directory + "\n")


def write_chips_runners(chips_input_file, runtag, band=None, field=None):
    frequency_band = 0 if band == "low" else 1
    identifier = f"{runtag}_multichips"

    write_run_clean_intermed("multichips", identifier)
    write_run_all_astro_script()
    write_make_obs_astro_script()

    with open(f"{identifier}.sh", "w") as script:
        script.write(
            f"""\
#!/bin/bash -l
#SBATCH --job-name={identifier}
#SBATCH --output={identifier}_%A.out
#SBATCH --clusters=garrawarla
#SBATCH --account=mwaeor
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jameskariuki31@gmail.com
module use /pawsey/mwa/software/python3/modulefiles
module load python-singularity
source /astro/mwaeor/kchege/virtenvs/karis/bin/activate
set -eux

sh make_obs_astro.sh {chips_input_file} {runtag} {frequency_band} multichips 0 1 {field}
sh run_all_astro.sh {chips_input_file} {runtag}_multichips 0 1

                    """
        )
    add_permissions(f"{identifier}.sh")
    return


if __name__ == "__main__":
    parser = ArgumentParser(
        "chips_multiple_obsids.py",
        description="running chips on a group of obsids",
    )

    group1 = parser.add_argument_group("general")
    group1.add_argument(
        "obsidsfile",
        # nargs="+",
        help="text file listing the obsids one per line",
    )

    group1.add_argument(
        "--uvfits_directory",
        default=None,
        type=str,
        required=True,
        help="A path to a directory with sub_directories labelled according to each obsid. In these sub_directories, let there be the uvfits to be read into chips",
    )

    group1.add_argument(
        "--runtag",
        default=None,
        type=str,
        required=True,
        help="The name (tag1) to label the chips run. Needed for chips to run. Tag2 is 'multichips' by default.",
    )

    group1.add_argument(
        "--band",
        choices=["low", "high"],
        default="low",
        type=str,
        required=False,
        help="low or high frequency band",
    )
    group1.add_argument(
        "--eorfield",
        choices=[0, 1],
        default=0,
        type=int,
        required=False,
        help="EoR0 or EoR1",
    )
    group1.add_argument(
        "--patchrun",
        dest="patchrun",
        action="store_true",
        help="are they uvfits form a patch? [default: false]",
    )

    args = parser.parse_args()

    logger = logging.getLogger(__file__)

    with open(args.obsidsfile) as file:
        lines = file.readlines()
        all_obsids = [line.rstrip() for line in lines]

    if len(all_obsids) == 0:
        logger.error("Obsids needed to process.")
        sys.exit(1)

    logger.info("%s obsids  in total", len(all_obsids))

    if args.patchrun:
        logger.info(
            "Chips does visibilities differencing requiring at least 2 timesteps. Is that true for patch patch uvfits? Otherwise it won't run."
        )

    directories = [
        f"{os.path.abspath(args.uvfits_directory)}/{obsid}" for obsid in all_obsids
    ]

    chips_config(all_obsids, directories, args.runtag, patchrun=args.patchrun)

    write_chips_runners(
        f"{args.runtag}_chips_inputfile.txt",
        args.runtag,
        band=args.band,
        field=args.eorfield,
    )
