#!/bin/bash
#SBATCH -J cellranger_arc
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=64G 
#SBATCH -t 00-00:00:00 
#SBATCH -o %x-%j.log  


/data/wangxin/software/cellranger-arc-2.0.2/cellranger-arc count --id=PBMC \
                     --reference=/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                     --libraries=/data/wangxin/work/MT_bm/10x/PBMC_arc/libraries.csv \
                     --localcores=12 \
                     --localmem=64
