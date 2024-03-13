from glob import glob
import os
import sys
import argparse
import itertools as it
import subprocess


parser = argparse.ArgumentParser(description="Submit jobs.",
									add_help=True)
parser.add_argument('-S', action='store', type=str, dest='Ss', default='100,200,300,400,500,600,700,800,900,1000',
                    help='Population size')
parser.add_argument('-G', action='store', type=str, dest='Gs', default='500',
                    help='G')                
parser.add_argument('-Dim', action='store', type=str, dest='Dims', default='5,10,15,25,50,100,125,150,175,200',
                    help='Dimension size')
parser.add_argument('-Damp', action='store', type=str, dest='Damps', default='1',
                    help='Dampening factor')
parser.add_argument('-MU', action='store', type=str, dest='MUs', default='0.01',
                    help='mutation rate')               
parser.add_argument('-epsilon', action='store', type=str, dest='eps', default='0,1,2,3.5,5',
                    help='epsilon')
parser.add_argument('-epsilon_type', action='store', type=str, dest='ep_types', default='0,1,2,3,4,5',
                    help='epsilon_type')
parser.add_argument('-max_loops', action='store', type=str, dest='loops', default='10000',
                    help='maximum number of loops')       
parser.add_argument('-Seed', action='store', type=str, dest='SEEDS',
					default='14724,24284,31658,6933,1318,16695,27690,8233,24481,6832,'
					'13352,4866,12669,12092,15860,19863,6654,10197,29756,14289,'
					'4719,12498,29198,10132,28699,32400,18313,26311,9540,20300')
parser.add_argument('-rdir', action='store', default='results', type=str,
                    help='Results directory') 
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

n_trials = len(args.SEEDS)
Seeds = args.SEEDS.split(',')[:n_trials]
Gs = args.Gs.split(',')
Ss = args.Ss.split(',')
Dims = args.Dims.split(',')
Damps = args.Damps.split(',')
MUs = args.MUs.split(',')
eps = args.eps.split(',')
ep_types = args.ep_types.split(',')
loops = args.loops.split(',')
args.slurm = True

# write run commands
all_commands = []
job_info = []
rdir = '/'.join([args.rdir])+'/stochastic results/'
os.makedirs(rdir, exist_ok=True)

for s,g,dim,damp,mu,seed,epsilon,epsilon_type,max_loops in it.product(Ss,Gs,Dims,Damps,MUs,Seeds,eps,ep_types,loops):

	all_commands.append(
		f'python3 single_lexicase_experiment.py -S {int(s)} -G {int(g)} -Dim {int(dim)} -Damp {int(damp)} -MU {float(mu)} -epsilon {float(epsilon)} -epsilon_type {int(epsilon_type)} -max_loops {int(max_loops)} -Seed {seed} -rdir {rdir}'
	)

	job_info.append({
         'S':s,
         'G':g,
         'Dim':dim,
         'Damp':damp,
         'MU': mu,  
         'epsilon': epsilon, 
		 'epsilon_type': epsilon_type,
         'max_loops': max_loops,
		 'Seed': seed,
		 'rdir': rdir
        
	})

print(len(job_info), 'total jobs created')
if args.slurm:
	# write a jobarray file to read commans from
	jobarrayfile = 'jobfiles/joblist.txt'
	os.makedirs(rdir + 'jobfiles', exist_ok=True)
	for i, run_cmd in enumerate(all_commands):

		job_name = '_'.join([x + '-' + f'{job_info[i][x]}' for x in
							['S','G','Dim', 'Damp', 'MU', 'epsilon', 'epsilon_type', 'Seed']])
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
