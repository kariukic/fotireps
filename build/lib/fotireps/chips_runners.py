def write_run_all_astro_script():
    with open("run_all_astro.sh", "w") as run_all_astro:
        run_all_astro.write(
            """\
#!/bin/bash -l

i=0

cnt=0
while read p; do

    if [ "$cnt" -ge ${3} ] && [ "$cnt" -lt ${4} ] ;
    then
 
obsid=$p

OBSIDCUT=$(echo $obsid| cut -d'/' -f 1)
echo $OBSIDCUT

if [ "$i" == '0' ]; then
echo First job
echo run_job${OBSIDCUT}DIFF.sh
JOB_alld=$(sbatch run_job${OBSIDCUT}DIFF.sh | cut -d " " -f 4)
echo $JOB_alld
fi
echo $i
if [ "$i" != '0' ]; then
JOB_alld=$(sbatch --dependency=afterok:$JOB_alld run_job${OBSIDCUT}DIFF.sh | cut -d " " -f 4)
fi

((i++))

    fi
    cnt=$((cnt+1))
done < ${1}

# run lssa stage

JOB_lssaxx=$(sbatch --dependency=afterok:$JOB_alld run_prepdiff.${2}.xx.sh | cut -d " " -f 4)
JOB_lssayy=$(sbatch --dependency=afterok:$JOB_alld run_prepdiff.${2}.yy.sh | cut -d " " -f 4)

# clean-up intermediate files

#sbatch --dependency=afterok:$JOB_lssaxx:$JOB_lssayy run_clean_intermed.diff.${2}.sh
sbatch --dependency=afterok:$JOB_lssaxx:$JOB_lssayy run_clean_intermed.sh

                    """
        )
    return


# make_obs_astro.sh
def write_make_obs_astro_script():
    with open("make_obs_astro.sh", "w") as make_obs_astro:
        make_obs_astro.write(
            """\
#!/bin/bash -l

#source /group/mwaops/setup_mwaops_galaxy.rc
source /astro/mwaeor/MWA/chips/scripts/env_variables_garrawarla.sh

#while read p; do
#  echo $p
#/home/ctrott/pbs_scripts/make_run_script_array.sh $p ${2} ${2}_${4}
#done < ${1}

cnt=0
while read p
do
    if [ "$cnt" -ge ${5} ] && [ "$cnt" -lt ${6} ] ;
    then
        echo $p
	/astro/mwaeor/MWA/chips/scripts/make_run_script_array_astro.sh $p ${2}_${4} ${3} ${7}
    fi
    cnt=$((cnt+1))
done < ${1}


# make prepare_lssa and run lssa PBS scripts

cat > run_prepdiff.${2}_${4}.xx.sh << Endofmessage
#!/bin/bash -l

#SBATCH --export=NONE
#SBATCH --job-name="lssa_xx"
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --output=myjob.%j.o
#SBATCH --error=myjob.%j.e
#SBATCH --account=mwaeor

source /astro/mwaeor/MWA/chips/scripts/env_variables_garrawarla.sh
export OMP_NUM_THREADS=24

printenv

cd /astro/mwaeor/MWA/chips/bin

./prepare_diff ${2}_${4} 384 0 'xx' ${2}_${4} ${3}
./lssa_fg_thermal ${2}_${4} 384 80 'xx' 300. ${2}_${4} 0 ${3} 0 8.

Endofmessage

cat > run_prepdiff.${2}_${4}.yy.sh << Endofmessage
#!/bin/bash -l

#SBATCH --export=NONE
#SBATCH --job-name="lssa_yy"
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --output=myjob.%j.o
#SBATCH --error=myjob.%j.e
#SBATCH --account=mwaeor

source /astro/mwaeor/MWA/chips/scripts/env_variables_garrawarla.sh
export OMP_NUM_THREADS=24

printenv

cd /astro/mwaeor/MWA/chips/bin

./prepare_diff ${2}_${4} 384 0 'yy' ${2}_${4} ${3}
./lssa_fg_thermal ${2}_${4} 384 80 'yy' 300. ${2}_${4} 0 ${3} 0 8.

Endofmessage


exit

                """
        )
    return


def write_run_clean_intermed(obsid, identifier):
    with open("run_clean_intermed.sh", "w") as run_clean_intermed:
        run_clean_intermed.write(
            f"""\
#!/bin/bash -l
#SBATCH --job-name=cl{obsid}
#SBATCH --output=cl{obsid}_%A.out
#SBATCH --clusters=garrawarla
#SBATCH --account=mwaeor
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --export=NONE
#SBATCH --mail-type=FAIL

set -eux

cd /astro/mwaeor/MWA/output
find . -maxdepth 1 -name "*freq*{identifier}*" -user kchege -print0 | xargs -0 rm -r
find . -maxdepth 1 -name '*syslog*{identifier}*' -user kchege -print0 | xargs -0 rm -r

echo 'cleaning finished successfully.'

                """
        )
    return