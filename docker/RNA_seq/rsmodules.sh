#!/bin/bash

module purge
module load FastQC/0.11.9
module load MultiQC/1.13-gimkl-2022a-Python-3.10.5
module load cutadapt/4.1-gimkl-2022a-Python-3.10.5
module load HISAT2/2.2.1-gimpi-2022a
module load SAMtools/1.15.1-GCC-11.3.0
module load Subread/2.0.3-GCC-11.3.0



echo "Loaded modules:FastQC,MultiQC,cutadapt,HISAT2,SAMtools,Subread"
