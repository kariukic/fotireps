"""This script writes bash scripts to run job on hpc slurm"""
import os
from fotireps.chips_runners import *


def add_permissions(bash_script):
    """Ensure permissions are sensible!

    Parameters
    ----------
    bash_script : [type]
        [description]
    """
    with open(bash_script, mode="a") as script:
        if bash_script == "run_jobs.sh":
            script.writelines(
                "echo 'All jobs fired successfully on garrawarla. Hoorooo!!'\n"
            )
            return
        else:
            script.writelines(
                """
            find . -user $USER -type d -exec chmod g+rwx,o+rx,o-w {} \;\n"""
            )
            script.writelines(
                "find . -user $USER -type f -exec chmod g+rw,o+r,o-w {} \;\n"
            )
            script.writelines("echo 'job finished successfully.'\n")
    return


def write_tiles_to_flag_file(tile_ids_to_flag):
    with open("flagged_tiles.txt", "w") as f:
        for item in tile_ids_to_flag:
            f.write("%s\n" % item)
    return


def write_run_shift_command():
    with open("run_shift_command.sh", "w") as run_shift_command:
        run_shift_command.write(
            """\
#!/bin/bash -l
module load hyperdrive
source /astro/mwaeor/kchege/virtenvs/karis/bin/activate

srclist=$1
peelsrclist=$2
patchsrclist=$3
shiftsjason=$4
metafits=$5

# removed this command because I now use srcshift-by-beam on the shifted peel srclist to make the shifted patch srclist
# srclist shift -v --collapse-into-single-source "${srclist}" "${shiftsjason}" "${patchsrclist}" --metafits "${metafits}"

hyperdrive srclist-shift -v "${srclist}" "${shiftsjason}" "${peelsrclist}" --include-unshifted-sources

echo "done"
                """
        )
    return


class CookScripts:
    def __init__(
        self,
        obsid="",
        run_patch=False,
        run_peel=False,
        sourcelist="",
        metafits=None,
        boxes_path="",
        virtual_env="",
        mail=None,
    ):
        self.obsid = obsid
        self.run_patch = run_patch
        self.run_peel = run_peel
        self.sourcelist = sourcelist
        self.metafits = metafits
        self.boxes_path = boxes_path
        self.virtual_env = virtual_env
        self.mail = mail

    def new_sh_script(self, job):
        """[summary]

        Parameters
        ----------
        job : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        script_name = f"{job}.sh"
        with open(script_name, mode="w") as script:
            script.writelines("#!/bin/bash -l\n")
            script.writelines(f"#SBATCH --job-name={job[:2]}{self.obsid}\n")
            script.writelines(f"#SBATCH --output={job}_%A.out\n")
            script.writelines("#SBATCH --clusters=garrawarla\n")
            script.writelines("#SBATCH --account=mwaeor\n")
            script.writelines("#SBATCH --export=NONE\n")
            if self.mail:
                script.writelines("#SBATCH --mail-type=FAIL \n")
                script.writelines("#SBATCH --mail-user={self.mail}\n")
        return script_name

    def rts_setup(
        self,
        DI_calibrators_numbers=None,
        DD_calibrators_numbers=None,
        integration_time=8,
        CorrDumpTime=2,  # tHIS IS THE time resolution of the input data, IN THIS CASE THE GPU BOX FILES
        cutoff=None,
        subbands=None,  # 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
        patch_time_config=None,
        no_srclist_by_beam=None,
        no_cotter_flags=None,
        use_fee_beam=None,
    ):
        rts_setup_script = self.new_sh_script("rts_setup")
        with open(rts_setup_script, mode="a") as script:
            script.writelines("#SBATCH --time=00:20:00\n")
            script.writelines("module use /pawsey/mwa/software/python3/modulefiles\n")
            script.writelines("module load python-singularity\n")
            script.writelines("module load mongoose/v0.2.5\n")
            script.writelines(f"{self.virtual_env}\n")
            script.writelines("set -eux\n")

            extras = ["--no-cotter-flags", "--use-fee-beam"]
            extra_options = ""
            for i, option in enumerate([no_cotter_flags, use_fee_beam]):
                if option:
                    extra_options += f" {extras[i]}"
            if extra_options:
                print(f"Using these extra options: {extra_options}")

            if self.run_patch:

                # The number of sources to make the patch(DI) sky model
                num_patch = DI_calibrators_numbers

                if not no_srclist_by_beam:
                    script.writelines(
                        f"srclist_by_beam.py -n {num_patch} --srclist {self.sourcelist} --metafits {self.metafits} --cutoff={cutoff} \n"
                    )

                if patch_time_config:
                    (
                        integration_time,
                        CorrDumpTime,
                    ) = patch_time_config  # eg 32 s integeration time, 2 seconds from the gpu boxes
                    corrDumpsPerCadence = (
                        integration_time // CorrDumpTime
                    )  # eg 32s integeration time/2 seconds from the gpu boxes
                    numberOfIterations = (
                        64 // integration_time
                    )  # eg.  64/32s integration time = 2

                    script.writelines(
                        f"rts-in-file-generator patch --base-dir {self.boxes_path} --metafits {self.metafits} --srclist *_patch{num_patch}.txt --num-iterations {numberOfIterations} --corr-dumps-per-cadence {corrDumpsPerCadence} --subband-ids {' '.join(str(id) for id in subbands)} --write-vis-to-uvfits {extra_options} -o rts_patch.in \n"
                    )
                else:
                    script.writelines(
                        f"rts-in-file-generator patch --base-dir {self.boxes_path} --metafits {self.metafits} --srclist *_patch{num_patch}.txt --subband-ids {' '.join(str(id) for id in subbands)} --write-vis-to-uvfits {extra_options} -o rts_patch.in \n"
                    )

                # A hack for now since mongoose does not have this command
                script.writelines('echo "DisableDipoleFlags=1" >> rts_patch.in \n')

            if self.run_peel:
                # This part seems necessary for the moongoose to edit the peel infile to tell the rts correctly the integration time we want.
                corrDumpsPerCadence = integration_time // CorrDumpTime
                # This calculates the totalnumber of timestamps
                numberOfIterations = 112 // integration_time
                # The different source numbers that the RTS uses in its DD calibration (Peel)
                # num_prepeel,
                num_full_dd_cals, num_iono, num_peel = DD_calibrators_numbers

                src_npeel = num_peel + 1000
                if not no_srclist_by_beam:
                    # Add 100 sources to numpeel just to make sure the rts veto passes
                    # src_npeel += 1000

                    script.writelines(
                        f"srclist_by_beam.py -n {src_npeel} --srclist {self.sourcelist} --metafits {self.metafits} --no_patch --cutoff={cutoff}\n"
                    )
                # --num-prepeel {num_prepeel} is not needed because the prepeel will always happen determined by the max in the other DD params (either num-peel or niono)
                # if no_cotter_flags:
                #     script.writelines(
                #         f"rts-in-file-generator peel --num-primary-cals {num_full_dd_cals} --num-cals {num_iono} --num-peel {num_peel} --num-iterations {numberOfIterations} --corr-dumps-per-cadence {corrDumpsPerCadence} --base-dir {self.boxes_path} --metafits {self.metafits} --srclist srclist_*_peel*.txt --subband-ids {' '.join(str(id) for id in subbands)}  --no-cotter-flags -o rts_peel.in\n"
                #     )  # TODO fix the * after srclist
                # else:
                script.writelines(
                    f"rts-in-file-generator peel --num-primary-cals {num_full_dd_cals} --num-cals {num_iono} --num-peel {num_peel} --num-iterations {numberOfIterations} --corr-dumps-per-cadence {corrDumpsPerCadence} --base-dir {self.boxes_path} --metafits {self.metafits} --srclist srclist_*_peel*.txt --subband-ids {' '.join(str(id) for id in subbands)} {extra_options} -o rts_peel.in\n"
                )  # TODO fix the * after srclist

                # A hack for now since mongoose does not have this command
                script.writelines('echo "DisableDipoleFlags=1" >> rts_peel.in \n')
        add_permissions(rts_setup_script)

    def rts_run(self):
        rts_run_script = self.new_sh_script("rts_run")
        with open(rts_run_script, mode="a") as script:
            script.writelines("#SBATCH --partition=gpuq\n")
            script.writelines("#SBATCH --gres=gpu:1\n")
            script.writelines("#SBATCH --nodes=25\n")
            script.writelines("#SBATCH --time=01:20:00\n")
            script.writelines("module use /pawsey/mwa/software/python3/modulefiles\n")
            script.writelines("module load python-singularity\n")
            script.writelines(f"{self.virtual_env}\n")
            script.writelines("module load RTS/sla_to_pal\n")
            script.writelines("set -eux\n")
            script.writelines("command -v rts_gpu\n")
            script.writelines("export UCX_MEMTYPE_CACHE=n\n")

            # if flag_tiles:
            #     script.writelines(
            #         f"printf '%s\\n' {' '.join(str(id) for id in flag_tiles)} >flagged_tiles.txt \n"
            #     )

            if self.run_patch:
                script.writelines("date\n")
                script.writelines("srun -n 25 --export=ALL rts_gpu rts_patch.in\n")
                script.writelines("date\n")
                script.writelines("set +e\n")
                script.writelines(
                    f"srun -n 1 -N 1 --export=ALL plot_BPcal_128T.py --both --outname {self.obsid}_BPcal.png\n"
                )
                script.writelines(
                    f"srun -n 1 -N 1 --export=ALL plot_CalSols.py --metafits {self.metafits}\n"
                )

                script.writelines(
                    f"srun -n 1 -N 1 --export=ALL rename_di_outputs.py uvdump_*.uvfits -p {self.obsid}_patch \n"
                )
            if self.run_peel:
                script.writelines("set +e\n")
                script.writelines("date\n")
                script.writelines("srun -n 25 --export=ALL rts_gpu rts_peel.in\n")
                script.writelines("date\n")
        add_permissions(rts_run_script)

    def cthulhu(self, logs_dir=None, radius=25, frequency=200):
        cthulhu_script = self.new_sh_script("cthulhu")
        with open(cthulhu_script, mode="a") as script:
            script.writelines(f"{self.virtual_env}\n")
            script.writelines("set -eux\n")
            # --frequency {frequency} no scaling the values that go into the logfile. for now
            # TODO fix this in cthulhu
            script.writelines(
                f"/astro/mwaeor/kchege/cthulhu/scripts/cthulhu_rts2format.py {logs_dir}/rts_mwa0*.log --filtering\n"  # 154.235
            )
            script.writelines(
                f"/astro/mwaeor/kchege/cthulhu/scripts/cthulhu_wrapper.py -p --frequency {frequency} -r {radius} {self.obsid}.yaml\n"
            )
        add_permissions(cthulhu_script)

        # script.writelines("mkdir rts_output\n")
        # script.writelines("mv *.dat peeled_* restore* *.log *.png  rts_output\n")
        # script.writelines("mkdir cthulhu_output\n")
        # script.writelines("mv plots *.yaml fits_files raw_and_tec cthulhu_output\n")
        # script.writelines("mkdir srclists_used\n")
        # script.writelines("mv *patch*.txt *peel*.txt srclists_used\n")

    def chips(
        self, uvfits_path=None, tags=None, band=None, field=None, uvfits_prefix=None
    ):
        chips_script = self.new_sh_script("chips")
        working_dir = os.path.abspath(".")
        with open(chips_script, mode="a") as script:
            script.writelines(f"{self.virtual_env}\n")
            script.writelines("set -eux\n")

            # script.writelines("mkdir chips\n")
            # script.writelines("cd chips\n")
            # change this with a command to write the script.
            # script.writelines("mv run_clean_intermed.sh chips\n")

            write_run_clean_intermed(self.obsid, f"{tags[0]}_{tags[1]}")
            write_run_all_astro_script()
            write_make_obs_astro_script()
            # script.writelines(
            #     "cp processing/chips_work/run_all_astro.sh /processing/chips_work/make_obs_astro.sh .\n"
            # )
            # tags[0]={obsid}_{npatch}patch_{npeel}peel
            script.writelines(f"echo '{self.obsid}/{tags[0]}' >chips_input_file.txt\n")
            script.writelines(
                f"mkdir -p /astro/mwaeor/MWA/data/{self.obsid}/{tags[0]}\n"
            )
            script.writelines(f"cd /astro/mwaeor/MWA/data/{self.obsid}/{tags[0]}\n")
            if uvfits_prefix:
                script.writelines(
                    f"ln -sf {uvfits_path}/{uvfits_prefix}_uvdump_*.uvfits .\n"
                )
            else:
                script.writelines(f"ln -sf {uvfits_path}/uvdump_*.uvfits .\n")
            script.writelines(f"cd {working_dir}\n")
            frequency_band = 0 if band == "low" else 1
            script.writelines(
                f"sh make_obs_astro.sh chips_input_file.txt {tags[0]} {frequency_band} {tags[1]} 0 1 {field}\n"
            )
            script.writelines(
                f"sh run_all_astro.sh chips_input_file.txt {tags[0]}_{tags[1]} 0 1\n"
            )
        add_permissions(chips_script)

    def run_jobs(self, jobs):
        if len(jobs) > 1:
            # At the moment this is all the jobs that fotireps can run
            jobs_in_sequential_running_order = [
                "download",
                "cotter",
                "rts_setup",
                "rts_run",
                "cthulhu",
                "chips",
            ]
            # Make a dict that assigns each job an index
            order_dict = {}
            for i, b in enumerate(jobs_in_sequential_running_order):
                order_dict[b] = i
            # Get the indices of the jobs needed to be run this time
            jobs_indices = [order_dict[ind] for ind in jobs]
            # sort the jobs
            ordered_jobs_list = [
                job
                for _idx, job in sorted(
                    zip(jobs_indices, jobs), key=lambda pair: pair[0]
                )
            ]
        else:
            ordered_jobs_list = jobs
        # Write the first job withou a "dependency:afterok"
        with open("run_jobs.sh", "w") as jobs_runner:
            job = ordered_jobs_list[0]
            jid = f"{job}_id"
            jid2 = f"{job}_id[3]"
            jsub = f"{job}_sub"
            jobs_runner.write(
                f"""\
#! /bin/bash
{jsub}="sbatch --nice=2 {job}.sh"
{jid}=($(${{{jsub}}}))
{jid}=${{{jid2}}} \n
        """
            )
        # Write the other jobs in a loop with "dependency:afterok"
        for n, job in enumerate(ordered_jobs_list[1:]):
            previous_job_index = ordered_jobs_list.index(job) - 1
            previous_job_id = f"{ordered_jobs_list[previous_job_index]}_id"
            this_job_id = f"{job}_id"
            this_job_id2 = f"{job}_id[3]"
            this_job_sub = f"{job}_sub"
            with open("run_jobs.sh", "a") as jobs_runner:
                jobs_runner.write(
                    f"""\

{this_job_sub}="sbatch --nice=2 --dependency=afterok:${{{previous_job_id}}} {job}.sh"
{this_job_id}=($(${{{this_job_sub}}}))
{this_job_id}=${{{this_job_id2}}}\n
                    """
                )
        add_permissions("run_jobs.sh")
