#!/bin/sh

#  FastQC.sh
#  Chip seq - quality control


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


baseDir=$HOME/path/to/my/data
rawData=$HOME/path/to/my/data
cd $rawData

echo "pwd: $(pwd)"
samples=(*.fastq.gz)

echo "samples: ${samples[*]}"
echo "running Fastqc on ${SLURM_JOBID}"
echo "pwd: $(pwd)"

mkdir ${baseDir}/fastqc_reports
mkdir ${baseDir}/fastqc_reports/untrimmedFastqc

echo "fastqc reports dir: ${baseDir}/fastqc_reports"
echo "untrimmed fastqc reports dir: ${baseDir}/fastqc_reports/untrimmedFastqc"


echo "---------FASTQC-----------------"


fastqc -o ${baseDir}/fastqc_reports/untrimmedFastqc  -t 24 *.fastq.gz

echo "command run = fastqc -o ${baseDir}/fastqc_reports/untrimmedFastqc -t 24 *.fastq.gz" 



conda deactivate
conda activate multiqc


echo "--------multiqc-----------------"


cd ${baseDir}/fastqc_reports/untrimmedFastqc
echo "pwd: $(pwd)"
multiqc .

echo " command run: multiqc ."


conda deactivate
conda activate chip


echo "-------Trim Galore!------------"


cd ${rawData}

echo "pwd: $(pwd)"
mkdir ${baseDir}/fastqc_reports/trimmedFastqc
echo "trimmed FastqcDir: ${baseDir}/fastqc_reports/trimmedFastqc"
mkdir ${baseDir}/trimmedData
echo "trimmed data dir: ${baseDir}/trimmedData"
trim_galore  --fastqc_args "--outdir ${baseDir}/fastqc_reports/trimmedFastqc" -o ${baseDir}/trimmedData *.fastq.gz
echo "trim_galore --fastqc_args --outdir ${baseDir}/fastqc_reports/trimmedFastqc -o ${baseDir}/trimmedData *.fastq.gz"



conda deactivate
conda activate multiqc


echo "-------multiqc on trimmed------"


cd ${baseDir}/fastqc_reports/trimmedFastqc
echo "pwd: $(pwd)"
multiqc .
echo "command run: multiqc ."



