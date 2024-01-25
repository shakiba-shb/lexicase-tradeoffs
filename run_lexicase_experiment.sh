#!/usr/bin/bash

#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name lexicase-tradeoffs

rdir="results"
ntrials=30 #change to 100
seeds=$(cat seeds.txt | head -n $ntrials)

mkdir -p $rdir

G=(500)
S=(100 200 300)
Dim=(5 15 25)
Damp=(1)
MU=(0.01)
epsilon=(0)

N=128
num_jobs="\j"  # The prompt escape for number of jobs currently running
count=0
for seed in ${seeds[@]} ; do
    for g in ${G[@]} ; do
        for s in ${S[@]} ; do
            for dim in ${Dim[@]} ; do
                for damp in ${Damp[@]} ; do
                    for mu in ${MU[@]} ; do
                        for eps in ${epsilon[@]} ; do
                            while (( ${num_jobs@P} >= N )); do
                                wait -n
                            done
                            ((++count))
                            {
                                job_name="G-${g}_S-${s}_Dim-${dim}_Damp-${damp}_MU-${mu}_Seed-${seed}_eps-${eps}" 
                                job_file="${rdir}/${job_name}.log" 

                                echo "job $count = ${job_name}..." 

                                python single_lexicase_experiment.py \
                                    -G ${g} \
                                    -S ${s} \
                                    -Dim ${dim} \
                                    -Damp ${damp}\
                                    -MU ${mu} \
                                    -Seed ${seed} \
                                    -epsilon ${eps} \
                                    -rdir ${rdir} \
                                    | tee -i ${job_file} >/dev/null

                                echo "$count completed ${job_file}..." 
                            } &
                        done
                    done
                done
            done
        done
    done
done

exit