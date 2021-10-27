import os
import sys

sys.path.extend(["/astro/mwaeor/kchege/fotireps/scripts"])
from scripts_writer import main


obsid = 1098111544
# sample_data_path = "/astro/mwaeor/kchege/fotireps/fotireps/data"
sample_data_path = "/astro/mwaeor/kchege/test_fotireps"
srclist = "/astro/mwaeor/kchege/mset_data/srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt"
boxes_path = f"/astro/mwaeor/kchege/pipeline2/boxesnflags/{obsid}"
obsid = 1098111544

DD_calibrators_numbers = (
    0,
    1000,
    1000,
)
num_full_dd_cals, num_iono, num_peel = DD_calibrators_numbers
numpatch = 1000
tint = 28
fovcut = 30
jobs = ["rts_setup", "rts_run", "cthulhu"]


def test_main():
    # os.system(f"rm -r {sample_data_path}/rts*.sh {sample_data_path}/cthulhu*.sh")
    # try:
    os.system(
        f"scripts_writer.py --obsid {obsid} --srclist {srclist} --job rts_setup rts_run cthulhu --output_directory {sample_data_path} --boxes_path {boxes_path} --numpatch {numpatch} --num_full_dd_cals {num_full_dd_cals} --numpeel {num_peel} --numiono {num_iono} --fov_cutoff {fovcut} --integration_time {tint}"
    )
    # except Exception as e:
    #     raise AssertionError("Could not cook scripts")
    return None
