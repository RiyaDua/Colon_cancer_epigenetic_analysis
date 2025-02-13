#!/bin/sh

#  computeMatrixmm10.sh
#
#
#  210823 DDMMYY

################################################
#Set options
################################################
#Set resource requirements
#SBATCH --time=24:00:00
#SBATCH --mail-type=end 
#SBATCH --mail-user=rdua5@jh.edu 
#SBATCH --cpus-per-task=6
#SBATCH --job-name Matrix
#SBATCH --mem=50G 
#SBATCH --partition=parallel
#SBATCH --account heaswar1

scratchDir=$HOME/scr4_heaswar1/Riya
date
echo "running ${SLURM_JOBID} job now"
hostname

#module load intel
#module load python/3.8
#module load deeptools/3.1.2

module restore

module load anaconda3

conda activate deepTools

mlStat=$?
echo mlStat : $mlStat
if [ $mlStat -eq 1 ]
then
  echo module load fail. Exiting
  exit 1
fi

module list

############################
baseDir=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ
#baseDir=$HOME/scr4_heaswar1/Sara/Nivi_ChIP_sjt/

echo "baseDir: $baseDir"

TF=("H3K27Ac" "H3K27me3" "H3K4me3")

echo "all TF: ${TF[*]}"

echo "TF selected currently: ${TF[$SLURM_ARRAY_TASK_ID]}"

outputBigWigDir=${baseDir}/bigwig

echo "outputBigWig: $outputBigWigDir"

allGenesBed=${baseDir}/mm10_allgenes.bed

###########################

mkdir ${baseDir}/deepToolsOutputTest/

outputDeepTools=${baseDir}/deepToolsOutputTest

echo "outputDeepTools: ${outputDeepTools}"

echo ---------------------
echo "---ComputeMatrix----"
echo ---------------------
# Get the bigWig files for the specific transcription factor, excluding those with "log2ratio" in their names
bwFiles=(${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw)
bwFiles=("${bwFiles[@]/%log2ratio.bw/}")  # Remove entries containing "log2ratio"

# Create the computeMatrix command
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles[@]} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed 

# Create the plotHeatmap command
plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu --samplesLabel Cre-1-H3K27Ac Cre1-H3K27me3 Cre1-H3K4me3 Cre3-H3K27Ac Cre3-H3K27me3 Cre3-H3K4me3 Cre5-H3K27Ac Cre5-H3K27me3 Cre5-H3K4me3 Ev1-H3K27Ac Ev1-H3K27me3 Ev1-H3K4me3 Ev3-H3K27Ac Ev3-H3K27me3 Ev3-H3K4me3 Ev5-H3K27Ac Ev5-H3K27me3 Ev5-H3K4me3


#bwFiles=$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | grep -v "log2ratio")

#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed 

#echo "command: computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed --samplesLabel "$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | sed 's/month//g')""

#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed

#echo "command: computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed"


echo ---------------------
echo "---PlotProfile-----"
echo ---------------------

#plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu --samplesLabel Cre1-H3K27Ac Cre1-H3K27me3 Cre1-H3K4me3 Cre3-H3K27Ac Cre3-H3K27me3 Cre3-H3K4me3 Cre5-H3K27Ac Cre5-H3K27me3 Cre5-H3K4me3 Ev1-H3K27Ac Ev1-H3K27me3 Ev1-H3K4me3 Ev3-H3K27Ac Ev3-H3K27me3 Ev3-H3K4me3 Ev5-H3K27Ac Ev5-H3K27me3 Ev5-H3K4me3


#echo "command: plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu  --samplesLabel "$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | sed 's/month//g')""


echo --------------------
echo "-----DONE---------"
echo --------------------

date
