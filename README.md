
# 🧬 Charting the Epigenetic Dynamics in Colon Carcinogenesis using ChIP-Seq Analysis Techniques

This repository contains the scripts and workflow for the **ChIP-Seq data analysis project** aimed at **mapping epigenetic changes associated with BRAF V600E mutations in colon cancer**.

## 📄 Project Overview
The goal of this project was to **investigate the chromatin landscape and identify epigenetic changes during colon carcinogenesis** driven by **BRAF V600E mutations**. **ChIP-Seq data** was analyzed to profile histone modifications and understand the epigenetic regulation of gene expression in colon cancer.

## 🧪 Key Steps in Analysis Pipeline

### 1️⃣ **Upstream Analysis (Data Preprocessing and Alignment)**
Performed on **HPC cluster (slurm)** using **Linux-based tools**:
- **Quality Control:** `FastQC`
- **Adapter Trimming:** `Trim Galore!`
- **Alignment:** `Bowtie2`
- **SAM to BAM Conversion, Sorting, Indexing:** `Samtools`
- **Peak Calling:** `MACS2`
- **BigWig Coverage Tracks:** `deepTools`

### 2️⃣ **Downstream Analysis (Differential Binding and Visualization)**
Performed in **R** using:
- **Differential Binding Analysis:** `DiffBind`
- **Visualization:** `ggplot2`

## 🛠️ Tools Used
| Phase               | Tools                                             |
|---------------------|---------------------------------------------------|
| Upstream Analysis   | `FastQC`, `Bowtie2`, `MACS2`, `Samtools`, `deepTools`, `Shell`, `HPC (slurm)` |
| Downstream Analysis | `R`, `DiffBind`, `ggplot2`                       |



## 📊 Key Outputs
- **Differential Binding Sites** (associated with BRAF V600E mutations)
- **ChIP-Seq Peak Profiles** for **histone modifications**
- **Heatmaps & Coverage Plots** of ChIP enrichment


## 🚀 Running the Pipeline
1. **Submit FastQC.sh** – Raw data quality control and trimming
2. **Submit ChipAlign.sh** – Read alignment, SAM→BAM conversion, indexing
3. **Submit computeMatrixmm10.sh** – Coverage matrix and heatmaps (deepTools)
4. **Submit ChromHMM.sh** – Chromatin state segmentation
5. **Run DiffBind in R** – Differential binding analysis

## ⚠️ Data Availability
**Due to privacy and publication restrictions**, raw data is **not provided**. This repository contains **generalizable scripts and analysis workflow**.


