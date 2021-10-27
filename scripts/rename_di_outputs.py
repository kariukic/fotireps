"""Renaming DI outpout files before they get over written by the DD outfits files"""
import os
import sys
from argparse import ArgumentParser
import logging


def rename(filename, prefix="patch"):
    """Renaming DI outpout files before they get over written by the DD outfits files

    Parameters
    ----------
    filesnames : [type]
        [description]
    prefix : [type]
        [description]
    """
    filename = os.path.abspath(filename)

    try:
        assert os.path.isfile(filename)
    except Exception as error:
        raise AssertionError(
            f"{filename} is not a valid file, cannot rename. Exiting!"
        ) from error

    newname = filename.replace(
        filename.split("/")[-1], f"{prefix}_{filename.split('/')[-1]}"
    )

    os.system(f"mv {filename} {newname}")

    return


if __name__ == "__main__":
    parser = ArgumentParser(
        "rename_di_outputs.py",
        description="Renaming DI outpout files before they get overwritten by the DD outfits files",
    )

    group1 = parser.add_argument_group("general")
    group1.add_argument(
        "filenames",
        nargs="+",
        help="DI ouptput files to rename",
    )

    group1.add_argument(
        "--prefix",
        "-p",
        default="patch",
        type=str,
        required=False,
        help="The prefix to add to the filenames",
    )
    args = parser.parse_args()

    logger = logging.getLogger(__file__)

    if len(args.filenames) == 0:
        logger.error("Aborting: Need files to process.")
        sys.exit(1)

    for filename in args.filenames:
        rename(filename, args.prefix)
