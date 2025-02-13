#!/bin/sh

#  matrix.sh
# Use deepTools to generate TSS heatmaps and profiles from bigWig coverage tracks.


scratchDir=$HOME/path
date
echo "running ${SLURM_JOBID} job now"
hostname


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


baseDir=$HOME/path

TF=("TF1" "TF2" "TF3") # hiddedn original names of TF to abide by the research privacy.
echo "all TF: ${TF[*]}"
echo "TF selected currently: ${TF[$SLURM_ARRAY_TASK_ID]}"

outputBigWigDir=${baseDir}/bigwig
echo "outputBigWig: $outputBigWigDir"


allGenesBed=${baseDir}/allgenes.bed


mkdir ${baseDir}/deepToolsOutputTest/
outputDeepTools=${baseDir}/deepToolsOutputTest
echo "outputDeepTools: ${outputDeepTools}"


echo "---ComputeMatrix----"

# Get the bigWig files for the specific transcription factor, excluding those with "log2ratio" in their names
bwFiles=(${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw)
bwFiles=("${bwFiles[@]/%log2ratio.bw/}")  # Remove entries containing "log2ratio"


# Creating the computeMatrix command
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles[@]} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed 

# Creating the plotHeatmap command
plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu --samplesLabel .....

bwFiles=$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | grep -v "log2ratio")
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed 
echo "command: computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${bwFiles} --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed --samplesLabel "$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | sed 's/month//g')""

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed
echo "command: computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -R ${allGenesBed} -S ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw --skipZeros -o ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -p 24 --outFileSortedRegions ${outputDeepTools}/regions_TSS_${TF[$SLURM_ARRAY_TASK_ID]}.bed"


echo "---PlotProfile-----"


plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu --samplesLabel ....
echo "command: plotHeatmap -m ${outputDeepTools}/matrix_${TF[$SLURM_ARRAY_TASK_ID]}_TSS.gz -out ${outputDeepTools}/TSS_${TF[$SLURM_ARRAY_TASK_ID]}_profile-heatmap.png --colorMap RdBu  --samplesLabel "$(ls ${outputBigWigDir}/*${TF[$SLURM_ARRAY_TASK_ID]}.bw | sed 's/month//g')""

