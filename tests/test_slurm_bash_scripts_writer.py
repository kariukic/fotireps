import os
import numpy as np
from fotireps.slurm_bash_scripts_writer import CookScripts, add_permissions

VIRTUAL_ENV = "source /astro/mwaeor/kchege/virtenvs/karis/bin/activate"

sample_data_path = "/astro/mwaeor/kchege/fotireps/fotireps/data"
os.chdir(sample_data_path)

obsid = 1098111544
boxes_path = f"/astro/mwaeor/kchege/pipeline2/boxesnflags/{obsid}"
metafitsfile = f"{boxes_path}/{obsid}.metafits"
# srclist = "/pawsey/mwa/software/python3/srclists/master/srclist_pumav3_EoR0aegean_fixedEoR1pietro+ForA_phase1+2.txt"
srclist = "/astro/mwaeor/kchege/mset_data/srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt"

DI_calibrators_numbers = 1000
DD_calibrators_numbers = (
    # args.numprepeel,
    0,
    1000,
    1000,
)
integration_time = 28
CorrDumpTime = 2
cutoff = 30
subbands = [int(sb) for sb in np.arange(1, 25)]
jobs = [
    "rts_setup",
    "rts_run",
    "cthulhu",
]


def dont_test_cook_scripts():
    scripts_writer = CookScripts(
        obsid=obsid,
        sourcelist=srclist,
        metafits=metafitsfile,
        boxes_path=boxes_path,
        virtual_env=VIRTUAL_ENV,
    )
    scripts_writer.rts_setup(
        DI_calibrators_numbers=DI_calibrators_numbers,
        DD_calibrators_numbers=DD_calibrators_numbers,
        integration_time=integration_time,
        CorrDumpTime=CorrDumpTime,
        cutoff=cutoff,
        subbands=subbands,
    )

    scripts_writer.rts_run()
    scripts_writer.cthulhu()
    scripts_writer.run_jobs(jobs)

    # os.system("sh run_jobs.sh")