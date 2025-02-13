
#!/bin/sh

#  ChromHMM.sh
#
#
#  040624 DDMMYY

################################################
#Set options
################################################
#Set resource requirements
#SBATCH --time=24:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rdua5@jhu.edu
#SBATCH --cpus-per-task=24
#SBATCH --job-name ChipAlignRPKM
#SBATCH --mem=90G
#SBATCH --partition=parallel
#SBATCH --account heaswar1
module restore


module load anaconda3
module load ChromHMM


mlStat=$?
echo mlStat : $mlStat
if [ $mlStat -eq 1 ]
then
  echo module load fail. Exiting
  exit 1
fi

module list

scratchDir=$HOME/scr4_heaswar1/Riya
date
echo "running ${SLURM_JOBID} job now"
hostname

#####################
#------Inputs-------#
#####################

#trimmed_data dir

baseDir=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ
mlStat=$?
echo mlStat : $mlStat
if [ $mlStat -eq 1 ]
then
  echo module load fail. Exiting
  exit 1
fi

module list

scratchDir=$HOME/scr4_heaswar1/Riya
date
echo "running ${SLURM_JOBID} job now"
hostname

#####################
#------Inputs-------#
#####################

#trimmed_data dir

baseDir=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ
