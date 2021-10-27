import os
import sys

sys.path.extend(["/astro/mwaeor/kchege/fotireps/scripts"])
from rename_di_outputs import rename

sample_data_path = "/astro/mwaeor/kchege/fotireps/fotireps/data"
obsid = 1098111544
prefix = "patch"


def dont_test_rename():
    try:
        dummy_patch_uvfits_file = f"{sample_data_path}/uv_dump_test.uvfits"
        if not os.path.isfile(dummy_patch_uvfits_file):
            os.system(f"touch {dummy_patch_uvfits_file }")
        assert os.path.isfile(dummy_patch_uvfits_file)
        rename(dummy_patch_uvfits_file, prefix)

        assert os.path.isfile(f"{sample_data_path}/{prefix}_uv_dump_test.uvfits")

    except Exception:
        raise AssertionError("Failing to rename patch files")
    return