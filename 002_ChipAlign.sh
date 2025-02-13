#!/bin/sh

#  ChipAlign.sh
# bash script for ChIP-seq data alignment and processing using tools like samtools, MACS2, and deepTools on an HPC cluster with SLURM.

module restore
module load samtools/1.10
module load anaconda3
conda activate Align

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





baseDir=$HOME/HOME/path/to/my/data

trimmedDataDir=${baseDir}/trimmedData
indexDir=$HOME/HOME/path/to/my/data
cd $trimmedDataDir
echo "pwd: $(pwd)"
fastq_R1=(*_trimmed.fq.gz)
filename=$(echo ${fastq_R1[$SLURM_ARRAY_TASK_ID]} | cut -f 4,5 -d '_')


echo "trimmedDatadir: $trimmedDataDir"
echo "indexDir: $indexDir"
echo "fastq_R1: ${fastq_R1[$SLURM_ARRAY_TASK_ID]}"
echo "filename= $filename"



outputUnalignedDir=$baseDir/bowtie2Alignment/unalignedReadsfq/
outputMetricsDir=$baseDir/bowtie2Alignment/bowtie2Metrics/
outputSamDir=$baseDir/bowtie2Alignment/alignedUnsortedSam/


cd $baseDir/bowtie2Alignment

echo "pwd: $(pwd)"


conda deactivate
conda activate samtools


echo "-------samtools--------------"

outputBamDir=$baseDir/bowtie2Alignment/alignedSortedBam/

echo "outputBamDir: $outputBamDir"
echo "filename: $filename"


echo "sam to bam conversion"



samtools view -hb -F 0x4  -x "XS" ${outputSamDir}${filename}_aligned_unsorted.sam | samtools sort > ${outputBamDir}${filename}_aligned_sorted.bam
echo "samtools view -hb -F 0x4  -x "XS" ${outputSamDir}${filename}_aligned_unsorted.sam | samtools sort > ${outputBamDir}${filename}_aligned_sorted.bam"


echo "----samtools index-----"


cd ${outputBamDir}
echo "pwd: $(pwd)"

samtools index -b ${filename}_aligned_sorted.bam  ${filename}_aligned_sorted.bam.bai
echo "samtools index -b -@ 24 ${filename}_aligned_sorted.bam  ${filename}_aligned_sorted.bam.bai"

conda deactivate
conda activate macs2


echo "-------------MACS2---------------"


cd ${outputBamDir}

echo "pwd: $(pwd)"


chipFiles=(*month*.bam)
inputFiles=(*Input*.bam)


echo "chipFiles: ${chipFiles[$SLURM_ARRAY_TASK_ID]}"
echo "inputFiles: ${inputFiles[0]}"


name=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | awk -F'_' '{print $1"_"$2}')
echo "filenameForMacsOutput: ${name}"

mkdir ${baseDir}/macs2Results/
outputMacsDir=${baseDir}/macs2Results/


echo "command run:  macs2 callpeak -t ${chipFiles[$SLURM_ARRAY_TASK_ID]} -c ${inputFiles[0]}  -g 1.87e9 -n ${name} --broad --outdir $outputMacsDir"
macs2 callpeak -t ${chipFiles[$SLURM_ARRAY_TASK_ID]} -c ${inputFiles[0]}  -g 1.87e9 -n ${name} --broad --outdir $outputMacsDir

conda deactivate
conda activate deepTools


echo "-----BamToBigWig------"




cd ${baseDir}/bowtie2Alignment/alignedSortedBam/
mkdir ${baseDir}/bigwig_RPKM/
outputBigWigDir=${baseDir}/bigwig_RPKM/
echo "outputBigWigDir: $outputBigWigDir"



bamCoverage --bam ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads 200
echo "bamCoverage --bam ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652906236 --ignoreForNormalization chrX --extendReads"


echo "----bamCompare-------"


cd ${baseDir}/bowtie2Alignment/alignedSortedBam/
chipFiles=(*month*.bam)
inputFiles=(*Input*.bam)


filenames=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | awk -F'_' '{print $1"_"$2}')
echo "filename: ${name}"


echo "chipFile: ${chipFiles[$SLURM_ARRAY_TASK_ID]}"
echo "inputFile: ${inputFiles[0]}"
echo "name: $name"

bamCompare -b1 ${chipFiles[$SLURM_ARRAY_TASK_ID]} -b2 ${inputFiles[0]} -o ${outputBigWigDir}${name}_log2ratio.bw
echo "bamCompare -b1 ${chipFiles[$SLURM_ARRAY_TASK_ID]} -b2 ${inputFiles[0]} -o ${outputBigWigDir}${name}_log2ratio.bw"
echo "finished $SLURM_ARRAY_TASK_ID"

conda deactivate

