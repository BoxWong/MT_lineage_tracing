reads_repo="/syn1/wangxin/reads/gb_rev/Ludwig/colonies/donor1"

mkdir -p star_results sorted_bam rmdup

while read line
do
    #STAR alignment
    echo STAR --genomeLoad LoadAndKeep --genomeDir /data/wangxin/ref/drop_seq/STAR --readFilesIn $reads_repo/${line}_1.fastq.gz $reads_repo/${line}_2.fastq.gz --readFilesCommand zcat --runThreadN 1 --outSAMtype BAM Unsorted --outFileNamePrefix ./star_results/${line}_

#    #sort bam
    echo samtools sort -@ 1 ./star_results/${line}_Aligned.out.bam -o ./sorted_bam/${line}_sorted.bam
    echo samtools index ./sorted_bam/${line}_sorted.bam
    #remove duplicates
    echo gatk MarkDuplicates --CREATE_INDEX true --REMOVE_DUPLICATES true --ADD_PG_TAG_TO_READS false --TMP_DIR . -I ./sorted_bam/${line}_sorted.bam -M ./rmdup/$line.metrics -O ./rmdup/$line.sorted.rmdup.bam


done < donor1_srr_accession.txt
