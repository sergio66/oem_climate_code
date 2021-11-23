## submit 40 jobs, one per latbin (each loops over 365 timesteps, SLOW!!)
/bin/rm slurm*.out
sbatch --exclude=cnode[204,225,231,267] --array=1-5      --begin=now        sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=6-10     --begin=now+15     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=11-15    --begin=now+30     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=16-20    --begin=now+45     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=21-25    --begin=now+60     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=26-30    --begin=now+75     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=31-35    --begin=now+90     sergio_matlab_jobB.sbatch 2
sbatch --exclude=cnode[204,225,231,267] --array=36-40    --begin=now+105    sergio_matlab_jobB.sbatch 2

