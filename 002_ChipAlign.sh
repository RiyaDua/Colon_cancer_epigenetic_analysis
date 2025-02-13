#!/bin/sh

#  ChipAlign.sh
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

#####################
#------Inputs-------#
#####################

#trimmed_data dir 

baseDir=$HOME/scr4_heaswar1/Riya/ftp.activemotif.com/FASTQ

trimmedDataDir=${baseDir}/trimmedData

#indexes
indexDir=$HOME/scr4_heaswar1/Riya

#enrichment

#this script is for histone ChIP peak calling - even though it says TF below
#I did not wish to change the variable name because of the graphing scripts I like to use after that use the same syntax 
#change the 'TF' variable below to your histone (no spaces etc)

#??????????
#TF=H3K27Ac
#make dir
#allGenesBed=${baseDir}/mm10-allCuratedRefSeq.bed

#####################

cd $trimmedDataDir

echo "pwd: $(pwd)"

fastq_R1=(*_trimmed.fq.gz) #_trimmed.fq.gz

#fastq_R2=(*_fixed_2.fq.gz)


filename=$(echo ${fastq_R1[$SLURM_ARRAY_TASK_ID]} | cut -f 4,5 -d '_')


echo "trimmedDatadir: $trimmedDataDir"

echo "indexDir: $indexDir"

echo "fastq_R1: ${fastq_R1[$SLURM_ARRAY_TASK_ID]}"

#echo "fastq_R2: ${fastq_R2[$SLURM_ARRAY_TASK_ID]}"

echo "filename= $filename"



#mkdir ${baseDir}/bowtie2Alignment

#mkdir ${baseDir}/bowtie2Alignment/unalignedReadsfq

#mkdir ${baseDir}/bowtie2Alignment/bowtie2Metrics

#mkdir ${baseDir}/bowtie2Alignment/alignedUnsortedSam



outputUnalignedDir=$baseDir/bowtie2Alignment/unalignedReadsfq/

outputMetricsDir=$baseDir/bowtie2Alignment/bowtie2Metrics/

outputSamDir=$baseDir/bowtie2Alignment/alignedUnsortedSam/



#echo "outputUnalignedFile: ${outputUnalignedDir}${filename}_unalignedReads.fq.gz"

#echo "outputMetricsFile: ${outputMetricsDir}${filename}_bowtie2Metrics.txt"

echo "outputSamFile: ${outputSamDir}${filename}_aligned_unsorted.sam"

echo ----------------------
echo ----------------------
echo ----------------------


#bowtie2 -p 24 -x ${indexDir}/mm10-index/mm10 -U ${fastq_R1[$SLURM_ARRAY_TASK_ID]}  -S ${outputSamDir}${filename}_aligned_unsorted.sam

echo "bowtie2 -p 24 -x ${indexDir}/mm10-index/mm10 -U ${fastq_R1[$SLURM_ARRAY_TASK_ID]}  -S ${outputSamDir}${filename}_aligned_unsorted.sam"

cd $baseDir/bowtie2Alignment

echo "pwd: $(pwd)"


conda deactivate
conda activate samtools

echo -------------------------------
echo "-------samtools--------------"
echo -------------------------------

#mkdir alignedUnsortedBam

#mkdir alignedSortedBam

outputBamDir=$baseDir/bowtie2Alignment/alignedSortedBam/

echo "outputBamDir: $outputBamDir"
echo "filename: $filename"

echo -----------------------
echo "sam to bam conversion"
echo -----------------------

#getting rid of unmapped sequences and multi-mappers. bowtie2 marks mulitmappers as 'XS'

#samtools view -hb -F 0x4  -x "XS" ${outputSamDir}${filename}_aligned_unsorted.sam | samtools sort > ${outputBamDir}${filename}_aligned_sorted.bam

echo "samtools view -hb -F 0x4  -x "XS" ${outputSamDir}${filename}_aligned_unsorted.sam | samtools sort > ${outputBamDir}${filename}_aligned_sorted.bam"

echo ----------------------
echo " unmapped removed"
echo "multimapped removed"
echo ----------------------

echo
echo

echo -------------------------
echo "----samtools index-----"
echo -------------------------

cd ${outputBamDir}

echo "pwd: $(pwd)"


#samtools index -b ${filename}_aligned_sorted.bam  ${filename}_aligned_sorted.bam.bai

echo "samtools index -b -@ 24 ${filename}_aligned_sorted.bam  ${filename}_aligned_sorted.bam.bai"

echo
echo

conda deactivate

conda activate macs2

echo ----------------------------------
echo "-------------MACS2---------------"
echo ----------------------------------

cd ${outputBamDir}

echo "pwd: $(pwd)"


chipFiles=(*month*.bam)

inputFiles=(*Input*.bam)


echo "chipFiles: ${chipFiles[$SLURM_ARRAY_TASK_ID]}"

echo "inputFiles: ${inputFiles[0]}"



##name=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | cut -f 1 -d '_')
name=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | awk -F'_' '{print $1"_"$2}')
echo "filenameForMacsOutput: ${name}"

mkdir ${baseDir}/macs2Results/

outputMacsDir=${baseDir}/macs2Results/


echo "command run:  macs2 callpeak -t ${chipFiles[$SLURM_ARRAY_TASK_ID]} -c ${inputFiles[0]}  -g 1.87e9 -n ${name} --broad --outdir $outputMacsDir"


#macs2 callpeak -t ${chipFiles[$SLURM_ARRAY_TASK_ID]} -c ${inputFiles[0]}  -g 1.87e9 -n ${name} --broad --outdir $outputMacsDir

conda deactivate
conda activate deepTools
echo ------------------------
echo "-----BamToBigWig------"
echo ------------------------



cd ${baseDir}/bowtie2Alignment/alignedSortedBam/

mkdir ${baseDir}/bigwig_RPKM/

outputBigWigDir=${baseDir}/bigwig_RPKM/

echo "outputBigWigDir: $outputBigWigDir"

#bamCoverage -b ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24 

#echo "bamCoverage -b ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24"

bamCoverage --bam ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads 200
echo "bamCoverage --bam ${name}_aligned_sorted.bam -o ${outputBigWigDir}${name}.bw -p 24 --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652906236 --ignoreForNormalization chrX --extendReads"

echo -----------------------
echo "----bamCompare-------"
echo -----------------------

cd ${baseDir}/bowtie2Alignment/alignedSortedBam/

chipFiles=(*month*.bam)

inputFiles=(*Input*.bam)

#filenames=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | cut -f 1 -d '_')
filenames=$(echo ${chipFiles[$SLURM_ARRAY_TASK_ID]} | awk -F'_' '{print $1"_"$2}')
echo "filename: ${name}"


echo "chipFile: ${chipFiles[$SLURM_ARRAY_TASK_ID]}"
echo "inputFile: ${inputFiles[0]}"
echo "name: $name"

bamCompare -b1 ${chipFiles[$SLURM_ARRAY_TASK_ID]} -b2 ${inputFiles[0]} -o ${outputBigWigDir}${name}_log2ratio.bw

echo "bamCompare -b1 ${chipFiles[$SLURM_ARRAY_TASK_ID]} -b2 ${inputFiles[0]} -o ${outputBigWigDir}${name}_log2ratio.bw"



echo "finished $SLURM_ARRAY_TASK_ID"

conda deactivate

echo "-------------"
echo "---DONE------"
echo "-------------"


date
