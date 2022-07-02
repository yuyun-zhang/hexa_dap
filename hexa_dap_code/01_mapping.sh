test -d 02_fastqc || mkdir 02_fastqc 
test -d 02_fastp || mkdir 02_fastp
test -d 03_bwa || mkdir 03_bwa
test -d 04_macs2 || mkdir 04_macs2
test -d 05_meme || mkdir 05_meme
function dap(){
    i=$1
    genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
    genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
    genomelen=14271578887
    index=~/yuyun/genome_data/wheat6_bwaindex/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
    halo=/public/home/zhangyuyun/yuyun/dap-seq/for_wangmeiyue/Halo-merge.sort.q20.rmdup.bam
    halo_bed=/public/home/zhangyuyun/yuyun/dap-seq/for_wangmeiyue/Halo-merge_PE_peaks.bed
    toppeak=6000

    fastqc -o 02_fastqc 01_rawdata/*${i}*_1.fq.gz 01_rawdata/*${i}*_2.fq.gz &
    fastp -w 10 -l 20 -q 20 -r --cut_right_window_size 4 --cut_right_mean_quality 20 --detect_adapter_for_pe \
        -i 01_rawdata/*${i}*_1.fq.gz -I 01_rawdata/*${i}*_2.fq.gz -o 02_fastp/${i}_1_val_1.fq.gz -O 02_fastp/${i}_2_val_2.fq.gz
    fastqc -o 02_fastp 02_fastp/${i}_1_val_1.fq.gz 02_fastp/${i}_2_val_2.fq.gz &
    
    cd 03_bwa
    file1=../02_fastp/${i}_1_val_1.fq.gz
    file2=../02_fastp/${i}_2_val_2.fq.gz
    bwa mem -t 10 ${index} $file1 $file2 | samtools view -bS -@ 10 - | samtools sort -@ 10 - > ${i}.sort.bam
    samtools stats ${i}.sort.bam >  ${i}.sort.bam.stats &
    samtools view -q 20 -b -@ 10  ${i}.sort.bam | samtools rmdup -  ${i}.sort.q20.rmdup.bam
    bedtools bamtobed -i ${i}.sort.q20.rmdup.bam | bedtools sort -i - > ${i}.sort.q20.rmdup.bed &
    samtools stats ${i}.sort.q20.rmdup.bam > ${i}.sort.q20.rmdup.bam.stats 
    scale=`grep 'reads mapped:'  ${i}.sort.q20.rmdup.bam.stats | cut -f3 | awk '{print 1000000/$1}'`
    bedtools genomecov -ibam  ${i}.sort.q20.rmdup.bam  -bg -split -scale ${scale} > ${i}.sort.q20.rmdup.rpm.bedgraph
    ~/yuyun/software/USCS/bin/wigToBigWig -clip ${i}.sort.q20.rmdup.rpm.bedgraph ${genomesize} ${i}.sort.q20.rmdup.rpm.bw &
    cd ../
    
    cd 04_macs2 
    macs2 callpeak -t ../03_bwa/${i}.sort.q20.rmdup.bam -c ${halo} -f BAMPE -g ${genomelen} --nomodel --nolambda -n ${i}_PE
    awk '$8>=10{print $1"\t"$2"\t"$3"\t"$4}' ${i}_PE_peaks.narrowPeak > ${i}_PE_peaks.p10.bed
    bedtools intersect -a ${i}_PE_peaks.p10.bed -b ${halo_bed} -v -sorted -wa > ${i}_PE_peaks.filter.p10.bed
    cut -f 4 ${i}_PE_peaks.filter.p10.bed | sort | join -1 1 - -2 10 <(grep -v '#' ${i}_PE_peaks.xls  | tail -n +3 | sort -k10,10) |\
        awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | sort -k9,9nr -k8,8nr > ${i}_PE_peaks.filter.p10.xls #sort by qvalue and then by enrichment
    cut -f 4  ${i}_PE_peaks.filter.p10.bed | sort | join -1 1 - -2 4 <(sort -k4,4 ${i}_PE_summits.bed) |\
        awk '{print $2"\t"$3"\t"$4"\t"$1}' | bedtools sort -i - > ${i}_PE_summits.filter.p10.bed 
    cd ../

    cd 05_meme
    ln -s ../04_macs2/${i}_PE_peaks.xls ./
    grep -v '#' ${i}_PE_peaks.xls | tail -n +3 | awk '{print $1"\t"$5"\t"$5+1}' | bedtools sort -i - | \
        bedtools slop -i - -g ${genomesize} -b 300 | tee ${i}.summitupd300.bed | seqtk subseq $genome - > ${i}.summitupd300.fa 
    grep -v '#' ${i}_PE_peaks.xls | tail -n +3 | head -n $toppeak | awk '{print $1"\t"$5"\t"$5+1}' | bedtools sort -i - | \
        bedtools slop -i - -g ${genomesize} -b 300 | tee ${i}.top${toppeak}.summitupd300.bed | seqtk subseq $genome - > ${i}.top${toppeak}.summitupd300.fa 
    cd ../
}

j=0
for tf in `cat tf.list`
do
    j=`expr $j + 1`
{
    echo "start $tf"
    dap ${tf}
    echo "finished $tf"
}&
    if [ $j == 3 ]
    then
        echo 'wait'
        j=0
        wait
    fi
done 
wait 

