#Preprocess: Old Donor 1 PBMC mtscATAC-seq dataset
cellranger-atac count --id=YFL \
                      --reference=/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/data/wangxin/reads/scATAC/PBMC/YFL \
                      --sample=Y_1,Y_2,Y_3,Y_4,Y_5,Y_6,Y_7,Y_8\
                      --localcores=16 \
                      --localmem=64 


#Preprocess: Old Donor 2 PBMC mtscATAC-seq dataset
cellranger-atac count --id=THQ \
                      --reference=/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/data/wangxin/reads/scATAC/PBMC/THQ \
                      --sample=T_1,T_2,T_3,T_4,T_5,T_6,T_7,T_7,T_8\
                      --localcores=16 \
                      --localmem=64 


