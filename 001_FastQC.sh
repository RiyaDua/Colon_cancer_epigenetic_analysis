#!/bin/sh

#  FastQC.sh
#  Chip seq - quality control
#
#  280524 DDMMYY

################################################
#Set options
################################################
#Set resource requirements
#SBATCH --job-name="Fastqc-Riya"
#SBATCH --time=72:00:00
#SBATCH --mail-type=end 
#SBATCH --mail-user=rdua5@jhu.edu 
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G 
#SBATCH --partition=parallel



# Start script
#module load python/3.8.6
#module load fastqc/0.12.1

module load anaconda3

mlStat=$?
echo mlStat : $mlStat
if [ $mlStat -eq 1 ]
then
  echo module load fail. Exiting
  exit 1
fi

module list


conda activate chip


scratchDir=$HOME/scratch
date
echo "running ${SLURM_JOBID} job now"
hostname

####################
#Change-Inputs-here#
####################

baseDir=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ
rawData=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ
cd $rawData

echo "pwd: $(pwd)"

samples=(*.fastq.gz)

echo "---------------------------"
echo "---------SAMPLES-----------"
echo "---------------------------"

echo "samples: ${samples[*]}"

echo "running Fastqc on ${SLURM_JOBID}"

echo "pwd: $(pwd)"

mkdir ${baseDir}/fastqc_reports

mkdir ${baseDir}/fastqc_reports/untrimmedFastqc

echo "fastqc reports dir: ${baseDir}/fastqc_reports"
echo "untrimmed fastqc reports dir: ${baseDir}/fastqc_reports/untrimmedFastqc"


echo ----------------------------------
echo "---------FASTQC-----------------"
echo ----------------------------------

fastqc -o ${baseDir}/fastqc_reports/untrimmedFastqc  -t 24 *.fastq.gz


echo "command run = fastqc -o ${baseDir}/fastqc_reports/untrimmedFastqc -t 24 *.fastq.gz" 

#################

conda deactivate

conda activate multiqc

echo ----------------------------------
echo "--------multiqc-----------------"
echo ----------------------------------

cd ${baseDir}/fastqc_reports/untrimmedFastqc

echo "pwd: $(pwd)"

multiqc .

echo " command run: multiqc ."

#################

conda deactivate

conda activate chip

echo "-------------------------------"
echo "-------Trim Galore!------------"
echo "-------------------------------"

cd ${rawData}

echo "pwd: $(pwd)"

mkdir ${baseDir}/fastqc_reports/trimmedFastqc

echo "trimmed FastqcDir: ${baseDir}/fastqc_reports/trimmedFastqc"

mkdir ${baseDir}/trimmedData

echo "trimmed data dir: ${baseDir}/trimmedData"

trim_galore  --fastqc_args "--outdir ${baseDir}/fastqc_reports/trimmedFastqc" -o ${baseDir}/trimmedData *.fastq.gz

echo "trim_galore --fastqc_args --outdir ${baseDir}/fastqc_reports/trimmedFastqc -o ${baseDir}/trimmedData *.fastq.gz"

#################

conda deactivate

conda activate multiqc

echo "-------------------------------"
echo "-------multiqc on trimmed------"
echo "-------------------------------"

cd ${baseDir}/fastqc_reports/trimmedFastqc

echo "pwd: $(pwd)"

multiqc .

echo "command run: multiqc ."


echo "-----------------------"
echo "---------DONE----------"
echo "-----------------------"

date

