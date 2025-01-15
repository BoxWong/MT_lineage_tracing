#!/bin/bash
#Xin Wang, x.wang@siat.ac.cn
#Calling MT mutations from SS2 data. 



gatk="/data/wangxin/software/miniconda3/envs/tzj/opt/gatk-3.8/GenomeAnalysisTK.jar"
hg38_ref="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
mtIndels="/data/wangxin/ref/gatk_hg38/chrM.phase3_callmom-v0_4.20130502.genotypes.vcf"
mt_bed="/data/wangxin/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/chrM.bed"

#Making a pseudobulk bam for identifying germline mutations.
sbatch -J merge -N 1 -n 16 -t 00:00:00 --mem=120G -o %x-%j.log --wrap="samtools merge -@ 16 -h sam.header -b bam.list -o embryo.sorted.rmdup.bam"

samtools view -h embryo.sorted.rmdup.bam chrM |samtools view -bS - > embryo_chrM.sorted.rmdup.bam
samtools reheader sam.header embryo_chrM.sorted.rmdup.bam > embryo_chrM.noRG.sorted.rmdup.bam
samtools addreplacerg -r ID:embryo -r LB:L1 -r SM:S1 -o embryo_chrM.sorted.rmdup.bam embryo_chrM.noRG.sorted.rmdup.bam
samtools index embryo_chrM.sorted.rmdup.bam
java -Xmx8g -jar ${gatk} -R ${hg38_ref} -T SplitNCigarReads -U ALLOW_N_CIGAR_READS -L $mt_bed -I embryo_chrM.sorted.rmdup.bam -o embryo_chrM.splitN.sorted.rmdup.bam
java -Xmx8g -jar ${gatk} -R ${hg38_ref} -T RealignerTargetCreator -known ${mtIndels} -L $mt_bed -nt 4 -I embryo_chrM.splitN.sorted.rmdup.bam -o embryo_possorted_chrM.target.intervals
java -Xmx10g -jar ${gatk} -R ${hg38_ref} -T IndelRealigner -maxReads 10000000 -known ${mtIndels} -L $mt_bed -I embryo_chrM.splitN.sorted.rmdup.bam -targetIntervals embryo_possorted_chrM.target.intervals -o embryo_possorted_chrM_realign.bam
samtools index embryo_possorted_chrM_realign.bam
samtools mpileup -l $mt_bed -q 30 -Q 30 -f $hg38_ref -x embryo_possorted_chrM_realign.bam > embryo_rmdup.mpileup
varscan pileup2snp embryo_rmdup.mpileup --min-var-freq 0.01 -min-reads2 2 > embryo.snv
perl /data/wangxin/src/scATAC/pileup_inf_rj.pl embryo_rmdup.mpileup > embryo.counts



#Calling MT mutations in single cells.
mkdir -p bam
mkdir -p mpileup
mkdir -p snv
mkdir -p counts

while read bam; do
        tmp_name=`basename $bam`
        primary_name=`echo ${tmp_name%%.*}`
        samtools view -h ${bam} chrM |sed -e '27,195d'|samtools view -bS > ./bam/${primary_name}_possorted_bam_chrM_noRG.bam
        #Modify bam files by addding RG tags of LB and SM.
        samtools addreplacerg -r ID:$primary_name -r LB:L1 -r SM:S1 -o ./bam/${primary_name}_possorted_bam_chrM.bam ./bam/${primary_name}_possorted_bam_chrM_noRG.bam && rm ./bam/${primary_name}_possorted_bam_chrM_noRG.bam
        samtools index ./bam/${primary_name}_possorted_bam_chrM.bam
        java -Xmx8g -jar ${gatk} -R ${hg38_ref} -T SplitNCigarReads -U ALLOW_N_CIGAR_READS -L $mt_bed -I ./bam/${primary_name}_possorted_bam_chrM.bam -o ./bam/${primary_name}_possorted_chrM.bam
        java -Xmx8g -jar ${gatk} -R ${hg38_ref} -T RealignerTargetCreator -known ${mtIndels} -L $mt_bed -nt 4 -I ./bam/${primary_name}_possorted_chrM.bam -o ./bam/${primary_name}_possorted_chrM.target.intervals
        java -Xmx10g -jar ${gatk} -R ${hg38_ref} -T IndelRealigner -maxReads 10000000 -known ${mtIndels} -L $mt_bed -I ./bam/${primary_name}_possorted_chrM.bam -targetIntervals ./bam/${primary_name}_possorted_chrM.target.intervals -o ./bam/${primary_name}_possorted_chrM_realign.bam
        samtools index ./bam/${primary_name}_possorted_chrM_realign.bam
        samtools mpileup -l $mt_bed -q 30 -Q 30 -f $hg38_ref -x ./bam/${primary_name}_possorted_chrM_realign.bam > ./mpileup/${primary_name}.rmdup.mpileup
        varscan pileup2snp ./mpileup/${primary_name}.rmdup.mpileup --min-var-freq 0.01 -min-reads2 2 > ./snv/${primary_name}.snv
        perl /data/wangxin/src/scATAC/pileup_inf_rj.pl ./mpileup/${primary_name}.rmdup.mpileup > ./counts/${primary_name}.counts
done < bam.list
