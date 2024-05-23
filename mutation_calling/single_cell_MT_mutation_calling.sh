#!/bin/bash
#Xin Wang, x.wang@siat.ac.cn
#Mitochondrial DNA mutation calling.

gatk="/data/wangxin/software/miniconda3/envs/tzj/opt/gatk-3.8/GenomeAnalysisTK.jar"
ref="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
mtIndels="/data/wangxin/ref/gatk_hg38/hg19.chrM.phase3_callmom-v0_4.20130502.genotypes.vcf"

#To generate a sam-like dict file for reference genome.
picard CreateSequenceDictionary -R ${ref}


#Extract chrM-specific reads and split reads into individual files.
input_bam="/data/wangxin/work/scATAC/PBMC/YFL/outs/possorted_bam.bam"
barcode_file="/data/wangxin/work/scATAC/PBMC/YFL/outs/filtered_peak_bc_matrix/barcodes.tsv"

samtools view -h ${input_bam} chrM |sed -e '27,195d'|samtools view -bS > possorted_bam_chrM.bam
samtools index possorted_bam_chrM.bam
java -Xmx8g -jar ${gatk} -R ${ref} -T RealignerTargetCreator -known ${mtIndels} -nt 8 -I possorted_bam_chrM.bam -o possorted_bam_chrM.target.intervals
java -Xmx10g -jar ${gatk} -R ${ref} -T IndelRealigner -maxReads 10000000 -known ${mtIndels} -I possorted_bam_chrM.bam -targetIntervals possorted_bam_chrM.target.intervals -o possorted_bam_chrM_realign.bam
samtools sort -@ 8 -t CB possorted_bam_chrM_realign.bam -o possorted_bam_chrM_realign_CB_sorted.bam
python /data/wangxin/src/scATAC/04_SplitBam.py -i possorted_bam_chrM_realign_CB_sorted.bam -b ${barcode_file}

#Call mutations from pseudobulk for identifying germline MT mutations.
mt_bed="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/chrM.bed"
bam="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/possorted_bam_chrM_realign.bam"
mpileup_output="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/YFL_all.mpielup"
snv_output="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/YFL_all.snv"
samtools mpileup -l $mt_bed -f $ref -q 20 -Q 30 $bam > $mpileup_output
varscan pileup2snp $mpileup_output --min-var-freq 0.01  --min-reads2 2 > $snv_output


#Call MT mutations in single cells.
bam_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/splitted_bam"
rmdupbam_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/rmdup_bam"
mt_bed="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/chrM.bed"
ref="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
mpileup_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/mpileup"
snv_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/SNV"
counts_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/counts"
src_dir="/data/wangxin/work/MT_bm/PBMC/tzj/YFL/src"
count=0
while read barcode 
do
    count=$((count + 1))
    flag=$((count % 500)) #500 cells per slurm job
    #echo $count
    if [[ $flag == 1 ]]; then
        echo -e "#!/bin/bash\n#SBATCH -J scATAC_$count\n#SBATCH -p all\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH --mem=10G \n#SBATCH -t 00:00:00 \n#SBATCH -o %x-%j.log"
    fi
    input_bam=$bam_dir/${barcode}.bam
    output_bam=$rmdupbam_dir/${barcode}.rmdup.bam
    echo -e "#Processing $input_bam"
    echo picard MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=$input_bam OUTPUT=$output_bam METRICS_FILE=$rmdupbam_dir/${barcode}.metrics
    echo samtools mpileup -l $mt_bed -q 30 -Q 30 -f $ref -x $output_bam \> $mpileup_dir/${barcode}.rmdup.mpileup
    echo varscan pileup2snp $mpileup_dir/$barcode.rmdup.mpileup --min-var-freq 0.01 -min-reads2 2 \> $snv_dir/${barcode}.snv
    echo perl /data/wangxin/src/scATAC/pileup_inf_rj.pl $mpileup_dir/${barcode}.rmdup.mpileup \> $counts_dir/${barcode}.counts
    echo
done < ${barcode_file} > slurm.src

#Split slurm job scripts and run them in parallel.
split -d -l 3008 slurm.src ${src_dir}/snv_src_  #3008 = 8+500*6
for i in `ls ${src_dir}/snv_src_*`
do
    sbatch $i
done

#Calculate mean depth for each cell after the previous jobs are completed.
for i in `ls $counts_dir`
do
    barcode=${i%.*}
    awk -v barcode=$barcode 'BEGIN {sum=0} {sum=sum+$4} END{print barcode","sum/NR}' $counts_dir/$i >> YFL_mean_depth.csv
done

