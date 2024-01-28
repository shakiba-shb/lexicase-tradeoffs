from glob import glob
import os
import sys
import argparse
import itertools as it
import subprocess


parser = argparse.ArgumentParser(description="Submit jobs.",
									add_help=True)
parser.add_argument('-experiment', action='store', type=str, dest='experiments', default='runid-0a2d09ad-13b7-48a3-95a5-226c5466104a.json',
                    help='Population size')   
parser.add_argument('-rdir', action='store', default='results', type=str,
                    help='Results directory')    
parser.add_argument('-datadir', action='store', type=str, 
					default='/home/shakiba/lexicase-tradeoffs/results/selected_files')
parser.add_argument('-Seed', action='store', type=str, dest='SEEDS',
					default='14724,24284,31658,6933,1318,16695,27690,8233,24481,6832,'
					'13352,4866,12669,12092,15860,19863,6654,10197,29756,14289,'
					'4719,12498,29198,10132,28699,32400,18313,26311,9540,20300')
parser.add_argument('-n_trials', action='store', dest='N_TRIALS', default=20,
					type=int, help='Number of trials to run')
parser.add_argument('-n_jobs', action='store', default=1,
					type=int, help='Number of parallel jobs')
parser.add_argument('-mem', action='store', dest='mem', default=1000, type=int,
					help='memory request and limit (MB)')
parser.add_argument('--slurm', action='store_true',
					default=False, help='Run on an slurm HPC')
parser.add_argument('-time', action='store', dest='time', 
					default='72:10:00', type=str, help='time in HR:MN:SS')
args = parser.parse_args()

fig = 'fig7'
datadir = args.datadir
args.slurm = True

print('running cag for experiments in', datadir) 

all_commands = []
job_info = []
rdir = '/'.join([args.rdir, fig])+'/cag_results/'
os.makedirs(rdir, exist_ok=True)

for filename in os.listdir(datadir):
    if filename.endswith('.json'):
        file_path = os.path.join(datadir, filename)

        all_commands.append(
            f'python3 single_cag_experiment.py -experiment {file_path} -rdir {rdir}'
        )

        job_info.append({
            'experiment': os.path.splitext(os.path.basename(file_path))[0],
            'rdir': rdir
        })

print(len(job_info), 'total jobs created')
if args.slurm:
	# write a jobarray file to read commans from
	jobarrayfile = 'cag_jobfiles/joblist.txt'
	os.makedirs(rdir + 'jobfiles', exist_ok=True)
	for i, run_cmd in enumerate(all_commands):

		job_name = job_info[i]['experiment']
		job_file = f'{rdir}/jobfiles/{job_name}.sb'
		out_file = job_info[i]['rdir'] + job_name + '_%J.out'

		batch_script = (
			f"""#!/usr/bin/bash 
#SBATCH -A ecode
#SBATCH --output={out_file} 
#SBATCH --job-name={job_name} 
#SBATCH --ntasks={1} 
#SBATCH --cpus-per-task={1} 
#SBATCH --time={args.time}
#SBATCH --mem={args.mem}

date
module purge
module restore new_modules
source lex/bin/activate
which python3
{run_cmd}

date
"""
		)

		with open(job_file, 'w') as f:
			f.write(batch_script)

		print(run_cmd)
		# sbatch_response = subprocess.check_output(
		# 	[f'sbatch {job_file}'], shell=True).decode()     # submit jobs
		# print(sbatch_response)
