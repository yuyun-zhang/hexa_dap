###################################################
############## 1. calculate TFAS ##################
###################################################
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
mkdir tmp/
x=0
for i in `cat ../../HC.tf.list`
do 
x=`expr $x + 1`
{
    echo start_$i
    summit=../../01_pass_peak/${i}_PE_summits.filter.p10.bed
    peak=../../01_pass_peak/${i}_PE_peaks.filter.p10.bed

    bed=../../01_bambed/${i}.sort.q20.rmdup.bed
    stats=../../01_bambed/${i}.sort.q20.rmdup.bam.stats 
    m=`grep 'reads mapped:' $stats | awk '{print $4}'` 

    bedtools intersect -a $summit -b $pro -wa -wb | awk '{print $8,$4,($7-$6+1)/2+$6-$2}' | awk '{sub("-","",$3);print}' > tmp/${i}.dis 
    sort -k2,2 tmp/${i}.dis | join -1 2 - -2 4 <(sort -k4,4 $peak) | \
        awk '{print $4,$5,$6,$1,$2,$3}' | tr ' ' '\t' | bedtools coverage -a - -b $bed -counts |\
        awk '{print $5,(($7*1000000000)/(($3-$2)*"'$m'"))*exp(-($6/2000))}' > tmp/${i}.TFAS.tmp  
    awk '{a[$1]=a[$1]+$2}END{for(i in a)print i,a[i]}' tmp/${i}.TFAS.tmp | sort -k1,1 | \
        join -1 1 - -2 1 <(cut -f 4 $pro | sort ) -a 2 | awk 'BEGIN{print "atf '$i'"}{if(NF==1)print $1,"0";else print $0}' | \
        sort -k1,1 > ${i}.TFAS
}&
    if [ $x == 3 ]
    then 
    x=0 
    wait 
    fi 
done 
wait 

###################################################
############## 2. balance type ####################
###################################################
trip=/public/home/zhangyuyun/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_7species_wheat4AB/wheat6.triplet.info.og
for i in `cat ../../HC.tf.list`
do 
    TFAS=../01_TFAS/${i}.TFAS 
    cut -f 2- -d ' ' ${trip} | sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 $TFAS) |\
        sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 $TFAS) | sort -k3,3 |\
        join -1 3 - -2 1 <(sort -k1,1 $TFAS) | awk '{print $3,$2,$1,$4,$5,$6}' > ${i}.TFAS.1-1-1.txt 
    awk '$4>0.25 || $5>0.25 || $6>0.25' ${i}.TFAS.1-1-1.txt > ${i}.TFAS.1-1-1.filter.txt
    n=`cat ${i}.TFAS.1-1-1.filter.txt | wc -l`
    awk '{print $1"-"$2"-"$3,$4/($4+$5+$6),$5/($4+$5+$6),$6/($4+$5+$6)}' ${i}.TFAS.1-1-1.filter.txt > ${i}.TFAS.1-1-1.filter.norm.txt
    Rscript get_balance.r ${i}.TFAS.1-1-1.filter.norm.txt
    if [ $n -ge 200 ]
    then 
        sort -k1,1 ${i}.TFAS.1-1-1.filter.norm.txt.type | join -1 1 - -2 2 <(sed 's/ /-/2;s/ /-/2' $trip | sort -k2,2) |\
            sed 's/-/ /;s/-/ /'  | awk '{print $5,$1,$4"\n"$5,$2,$4"\n"$5,$3,$4}' > ${i}.triplet.type
        echo $i >> HC.tf.trip-gt200.list  
    fi 
done 

###################################################
############## 3. overlap with TE #################
###################################################
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
TE=~/yuyun/genome_data/cs_TE_split.bed
mkdir b-ub_bed/ 
for i in `cat ../02_balance_stats/HC.tf.trip-gt200.list`
do 
    ##### blanced/unbalanced triplet - TFBS 
    summit=../../01_pass_peak/${i}_PE_summits.filter.p10.bed
    grep balanced ../02_balance_stats/${i}.triplet.type | sort -k2,2 | \
        join -1 2 - -2 4 <(sort -k4,4 $pro) | awk '{print $4,$5,$6,$1,$2,$3}' | tr ' ' '\t' |\
        bedtools intersect -a $summit -b - -wa -wb | sort -u | bedtools sort -i - |\
        awk '{print $1,$2,$3,$4,$8,$9,$10}' | tr ' ' '\t' > b-ub_bed/${i}_balance.summit.bed 
    grep -v balanced ../02_balance_stats/${i}.triplet.type | sort -k2,2 | \
        join -1 2 - -2 4 <(sort -k4,4 $pro) | awk '{print $4,$5,$6,$1,$2,$3}' | tr ' ' '\t' |\
        bedtools intersect -a $summit -b - -wa -wb | sort -u | bedtools sort -i - |\
        awk '{print $1,$2,$3,$4,$8,$9,$10}' | tr ' ' '\t' > b-ub_bed/${i}_unbalance.summit.bed 
    
    peak=../../01_pass_peak/${i}_PE_peaks.filter.p10.bed
    bedtools intersect -a $peak -b b-ub_bed/${i}_balance.summit.bed -wa -wb | awk '$4==$8' |\
        cut -f 1-4,9,10,11 > b-ub_bed/${i}_balance.peak.bed 
    bedtools intersect -a $peak -b b-ub_bed/${i}_unbalance.summit.bed -wa -wb | awk '$4==$8' |\
        cut -f 1-4,9,10,11 > b-ub_bed/${i}_unbalance.peak.bed 
    
    ###### TE overlap (with at least 1 TE in triplet promoter)
    bedtools intersect -a b-ub_bed/${i}_balance.summit.bed -b $TE -wa -wb | cut -f 1-7,11 > b-ub_bed/${i}_balance.summit.TE.bed
    bedtools intersect -a b-ub_bed/${i}_unbalance.summit.bed -b $TE -wa -wb | cut -f 1-7,11 > b-ub_bed/${i}_unbalance.summit.TE.bed
    bedtools intersect -a b-ub_bed/${i}_balance.peak.bed -b $TE -wo | cut -f 1-7,11,14 > b-ub_bed/${i}_balance.peak.TE.bed
    bedtools intersect -a b-ub_bed/${i}_unbalance.peak.bed -b $TE -wo | cut -f 1-7,11,14 > b-ub_bed/${i}_unbalance.peak.TE.bed

    ###### stats: TEnumber 
    n=`cut -f 6 b-ub_bed/${i}_balance.summit.TE.bed | sort -u | wc -l`
    cut -f 1,6 b-ub_bed/${i}_balance.summit.TE.bed | sort -u | cut -f 2 | sort | uniq -c | \
        awk '{print "TE"$1}' | sort | uniq -c | awk '{print "'$i' balanced TE",$2,$1,$1/"'$n'"}' >> TEnumber.b+ub.TE.txt 
    n=`cut -f 6 b-ub_bed/${i}_unbalance.summit.TE.bed | sort -u | wc -l`
    cut -f 1,6 b-ub_bed/${i}_unbalance.summit.TE.bed | sort -u | cut -f 2 | sort | uniq -c | \
        awk '{print "TE"$1}' | sort | uniq -c | awk '{print "'$i' unbalanced TE",$2,$1,$1/"'$n'"}' >> TEnumber.b+ub.TE.txt
done 

grep TE stats.TE.relic.ratio | grep balanced -w | awk '$5>20' > temp1 
grep TE stats.TE.relic.ratio | grep unbalanced -w | awk '$5>20' > temp2 
cut -f 1 temp1  -d ' ' | cat - <(cut -f 1 temp2 -d ' ') | sort | uniq -d | fgrep -f -  TEnumber.b+ub.TE.txt -w \
    > TEnumber.b+ub.TE.filter20.txt
rm temp1 temp2

#######################################################
## 4 regulation divergence vs. expression divergence ##
#######################################################
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
trip=~/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_7species_wheat4AB/wheat6.triplet.info

for i in `cat ../../HC.tf.list`
do 
    echo start_$i
    summit=../../01_pass_peak/${i}_PE_summits.filter.p10.bed
    peak=../../01_pass_peak/${i}_PE_peaks.filter.p10.bed
    bed=../../01_bambed/${i}.sort.q20.rmdup.bed
    TFAS=../../04_triplet_balance_relic/02_balance_stats/${i}.TFAS.1-1-1.txt

    ##### target 
    bedtools coverage -a $pro -b $bed -counts | cut -f 4,6 | sort -k1,1 > ${i}.tssupd5k.count 
    join -1 1 <(sort -k1,1 $trip) -2 1 ${i}.tssupd5k.count -a 1 | awk '{if(NF==3)print $0,"0";else print $0}' |\
        sort -k2,2 | join -1 2 - -2 1 ${i}.tssupd5k.count -a 1 | awk '{if(NF==4)print $0,"0";else print $0}' |\
        sort -k3,3 | join -1 3 - -2 1 ${i}.tssupd5k.count -a 1 | awk '{if(NF==5)print $0,"0";else print $0}' |\
        awk '{print $3"-"$2,$4,$5 > "'${i}'.A-B.tssupd5k.count"}{print $3"-"$1,$4,$6 > "'${i}'.A-D.tssupd5k.count"}{print $2"-"$1,$5,$6 > "'${i}'.B-D.tssupd5k.count"}'
    awk '{print $1"-"$2,$4,$5 > "'${i}'.A-B.TFAS"}{print $1"-"$3,$4,$6 > "'${i}'.A-D.TFAS"}{print $2"-"$3,$5,$6 > "'${i}'.B-D.TFAS"}' $TFAS 
    
    for j in A-B A-D B-D 
    do 
        sort -k1,1 ${i}.${j}.tssupd5k.count | join -1 1 - -2 1 <(sort -k1,1 ${i}.${j}.TFAS) |\
            awk 'BEGIN{print "trip trip1 trip2"}$4>0.25||$5>0.25{print $1,$2,$3}' > ${i}.${j}.tssupd5k.filter.count
        Rscript work_edger.r ${i}.${j}.tssupd5k.filter.count 1 1
    done 
done 

#### target matrix 
for j in A-B A-D B-D 
do 
    cat *.${j}.tssupd5k.filter.count | grep -v trip | cut -f 1 -d ' ' | sort -u > ${j}.trip.list 
    for i in `cat ../../HC.tf.list`
    do 
        sort -k1,1 edgeR.${i}.${j}.tssupd5k.filter.count.result | grep -v logFC | awk '{print $1,$2}' |\
            join -1 1 <(sort -k1,1 ${j}.trip.list) -2 1 - -a 1 | awk '{if(NF==1)print $0,"0";else print $0}' |\
            sort -k1,1 | awk 'BEGIN{print "'$i'"}{print $2}' > ${i}.${j}.log2fc.temp 
    done 
    Rscript work_get_abs_sum.r ${j}.log2fc.matrix
done 
rm *temp 

#### expression 
count=/public/home/zhangyuyun/yuyun/dap-seq/01_RNA-seq_CS/wheat6.seedling.count
fpkm=/public/home/zhangyuyun/yuyun/dap-seq/01_RNA-seq_CS/wheat6.seedling.fpkm 
join -1 1 <(sort -k1,1 $trip) -2 1 <(sort -k1,1 $count) -a 1 | awk '{if(NF==3)print $0,"0 0";else print $0}' |\
    sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 $count) -a 1 | awk '{if(NF==5)print $0,"0 0";else print $0}' |\
    sort -k3,3 | join -1 3 - -2 1 <(sort -k1,1 $count) -a 1 | awk '{if(NF==7)print $0,"0";else print $0}' |\
    awk '{print $3"-"$2,$4,$5,$6,$7 > "A-B.expr.count"}{print $3"-"$1,$4,$5,$8,$9 > "A-D.expr.count"}{print $2"-"$1,$6,$7,$8,$9 > "B-D.expr.count"}'
join -1 1 <(sort -k1,1 $trip) -2 1 <(sort -k1,1 $fpkm) -a 1 | awk '{if(NF==3)print $0,"0 0";else print $0}' |\
    sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 $fpkm) -a 1 | awk '{if(NF==5)print $0,"0 0";else print $0}' |\
    sort -k3,3 | join -1 3 - -2 1 <(sort -k1,1 $fpkm) -a 1 | awk '{if(NF==7)print $0,"0";else print $0}' |\
    awk '{print $3"-"$2,$4,$5,$6,$7 > "A-B.expr.fpkm"}{print $3"-"$1,$4,$5,$8,$9 > "A-D.expr.fpkm"}{print $2"-"$1,$6,$7,$8,$9 > "B-D.expr.fpkm"}'

for j in A-B A-D B-D 
do 
    sort -k1,1 ${j}.expr.count | join -1 1 - -2 1 <(sort -k1,1 ${j}.expr.fpkm) |\
        awk 'BEGIN{print "trip trip1-1 trip1-2 trip2-1 trip2-2"}($6+$7)/2>1||($8+$9)/2>1{print $1,$2,$3,$4,$5}' > ${j}.expr.filter.count
    Rscript work_edger.r ${j}.expr.filter.count 2 2 
done 

#### plot: quantile box 
rm plot.regu-log2fc-level.txt 
rm plot.regu-log10p-level.txt 
for j in A-B A-D B-D 
do 
    Rscript work_get_quantile.r ${j}.log2fc.matrix.abssum 10 up 

    for i in {1..10}
    do 
        n1=`awk 'NR=="'$i'"{print $2*1}' ${j}.log2fc.matrix.abssum.quantile`
        n2=`awk 'NR=="'$i'"{print $3*1}' ${j}.log2fc.matrix.abssum.quantile`

        sort -k1,1 ${j}.log2fc.matrix.abssum | join -1 1 - -2 1 <(sort -k1,1 ${j}.expr.log2fc) | awk '{if($3<0)print $1,$2,$3*(-1);else print $0}' | \
            awk -v n1=$n1 -v n2=$n2 '$2>=n1&&$2<=n2{print $0" level'$i'","'$j'"}' >> plot.regu-abslog2fc-level.txt 
    done 
done 


###################################################
############## 5. definition of dTE ###############
###################################################
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
trip=~/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_7species_wheat4AB/wheat6.triplet.txt
TE=~/yuyun/genome_data/cs_TE_split.bed
genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta

mkdir TE_bed/
mkdir blast_fastav2/
mkdir tmp/
x=0
for i in `cat ../../04_triplet_balance_relic/02_balance_stats/HC.tf.trip-gt200.list`
do 
    test -d blast_fastav2/${i}/ || mkdir blast_fastav2/${i}/
    bal=../../04_triplet_balance_relic/03_TE_relic/b-ub_bed/${i}_balance.summit.bed
    unbal=../../04_triplet_balance_relic/03_TE_relic/b-ub_bed/${i}_unbalance.summit.bed

    bedtools intersect -a $TE -b $bal -wa -wb | cut -f 1-4,11,12,13 | sort -u | bedtools sort -i - > TE_bed/${i}.balance.TE.bed 
    bedtools intersect -a $TE -b $unbal -wa -wb | cut -f 1-4,11,12,13 | sort -u | bedtools sort -i - > TE_bed/${i}.unbalance.TE.bed 

    for j in balance unbalance 
    do 
        cut -f 6 TE_bed/${i}.${j}.TE.bed | sort -u | while read line 
        do 
            x=`expr $x + 1`
        {
            grep $line TE_bed/${i}.${j}.TE.bed | bedtools getfasta -fi $genome -bed - -name > blast_fastav2/${i}/${j}.${line}.TE.fa 
            grep $line TE_bed/${i}.${j}.TE.bed | cut -f 4 > tmp/${i}.${j}.${line}.te 
            grep $line TE_bed/${i}.${j}.TE.bed | cut -f 5 | fgrep -f - <(grep $line $trip) -v | cut -f 2 -d ' ' |\
                fgrep -f - $pro | bedtools subtract -a - -b <(grep tmp/${i}.${j}.${line}.te $TE -w) |\
                bedtools getfasta -fi $genome -bed - -name > blast_fastav2/${i}/${j}.${line}.nonTE.fa 
            blastn -query blast_fastav2/${i}/${j}.${line}.TE.fa -subject blast_fastav2/${i}/${j}.${line}.nonTE.fa -outfmt 6 -out blast_fastav2/${i}/${j}.${line}.out 
        }&
        if [ $x == 20 ]
        then 
            x=0
            wait 
        fi 
        done 
        wait 
        echo finished_blast_${i}_${j}
        cat blast_fastav2/${i}/${j}.*.out | awk '$4>50' | cut -f 1,2,7,8,9,10 | sed 's/::/\t/g;s/:/\t/g;s/-/\t/g' | \
            awk '{if($11<$12)print $2,$3+$9,$3+$10,$1,$6,$7+$11,$7+$12,$5;else print $2,$3+$9,$3+$10,$1,$6,$7+$12,$7+$11,$5}' | \
            tr ' ' '\t' > blast_fastav2/${i}.${j}.TE-dTE.bedpe 
        awk '{print $5,$6,$7,$8"::"$4}' blast_fastav2/${i}.${j}.TE-dTE.bedpe | tr ' ' '\t' | bedtools sort -i - > blast_fastav2/${i}.${j}.dTE.bed
        summit=../02_TE_ratio/${i}_${j}.summit.bed
        peak=../02_TE_ratio/${i}_${j}.peak.bed
        bedtools intersect -a $summit -b blast_fastav2/${i}.${j}.dTE.bed -wa -wb | sed 's/::/ /' | awk '$5==$11' | sed 's/ /::/' | \
            cut -f 1-4,5,6,7,11 | sed 's/TraesCS[1-7][ABD]02G[0-9]*:://' > blast_fastav2/${i}.${j}.summit.dTE.bed
        bedtools intersect -a $peak -b blast_fastav2/${i}.${j}.dTE.bed -wo | sed 's/::/ /' | awk '$5==$11' | sed 's/ /::/' |\
            cut -f 1-4,5,6,7,11,12 | sed 's/TraesCS[1-7][ABD]02G[0-9]*:://' > blast_fastav2/${i}.${j}.peak.dTE.bed
    done 
done 
rm -rf tmp/${i}.${j}.${line}.te

##### stats: TE number 
cut -f 1 ../02_TE_ratio/TEnumber.b+ub.TE.filter20.txt -d ' ' | sort -u > HC.TE-triplet-gt20.list 
rm stats.TEnumber.TE+dTE.txt
for i in `cat HC.TE-triplet-gt20.list`
do 
    for j in balance unbalance 
    do 
        TE=../02_TE_ratio/${i}_${j}.summit.TE.bed
        dTE=blast_fastav2/${i}.${j}.summit.dTE.bed #### dTE-TFBS 

        n=`cat $TE $dTE | cut -f 6 | sort -u | wc -l`
        cat $TE $dTE | cut -f 1,6 | sort -u | cut -f 2 | sort | uniq -c | awk '{print "TE+dTE"$1}' | sort | uniq -c |\
            awk '{print "'$i' '$j'",$2,$1,$1/"'$n'"}' >> stats.TEnumber.TE+dTE.txt
    done 
done 

###################################################
############## 6. dTE conservation ################
###################################################
cat ../../04_triplet_balance_triplet_blast/03_TE_blast_promoter/blast_fasta/*.balance.summit.dTE.bed |\
    bedtools sort -i - | cut -f 1-3 | bedtools merge -i - -d 300 > balance.dTE.bed 
#### control1: nonTE promoter TFBS 
trip=~/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_7species_wheat4AB/wheat6.triplet.info
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
TE=~/yuyun/genome_data/cs_TE_split.bed
TFBS=../../01_pass_peak/mergeTFBS_PE_summits.filter.p10.bed

bedtools intersect -a $pro -b $TE -wa -v | cut -f 4 | sort -u > nonTE.gene.list 
join -1 1 nonTE.gene.list -2 1 <(sort -k1,1 $trip) | sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 nonTE.gene.list) |\
    sort -k3,3 | join -1 3 - -2 1 <(sort -k1,1 nonTE.gene.list) | tr ' ' '\n' |\
    fgrep -f - $pro | bedtools intersect -a $TFBS -b - -wa | sort -u | bedtools sort -i - | cut -f 1-3 > control.nonTE.bed 

######## conservation score 
genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
for file in `ls *TE.bed`
do 
    name=`basename $file .bed | sed 's/\./-/'`
    bedtools slop -g $genomesize -i $file -b 2000 | bedtools makewindows -b - -w 50 -i winnum | \
        bedtools sort -i - > ${file/.bed/}.centerupd2k.w50.bed
    dir=~/yuyun/dap-seq/04_epigenetic_data/02_CS/conservation_score/bed/
    x=0
    for chr in `cat chrlist`
    do 
        x=`expr $x + 1`
    {
        bedtools map -a <(awk '$1~"'$chr'"' ${file/.bed/}.centerupd2k.w50.bed) -b ${dir}/${chr}.conserved_score.bed -o mean |\
            awk '{print "'$name'",$4,$5}' > $name.${chr}.conservation.txt 
        echo finished_${name}_${chr} 
    }&
    if [ $x == 21 ]
    then 
        x=0
        wait 
    fi 
    done 
    wait 
    cat $name.chr*.conservation.txt | tee $name.conservation.txt | sed 's/\.$/0/' | \
        awk '{a[$1" "$2]=a[$1" "$2]+$3;b[$1" "$2]=b[$1" "$2]+1}END{for(i in a)print i,a[i]/b[i]}' >> plot.conservation.txt
    rm $name.chr*.conservation.txt
done 

###################################################
############## 7. TE/dTE enrichment ###############
###################################################
dir1=../02_TE_ratio/
dir2=../03_TE_blast_promoter/blast_fasta/
TEinfo=/public/home/zhangyuyun/yuyun/genome_data/TE_annotation/clariTE/wheat6_TE.ingenome.txt
rm enrichment*txt 
for i in `cat ../01_balance_stats/HC.tf.trip-gt200.list`
do 
    ##### TE:length 
    for n in balance unbalance 
    do 
        a=`awk '{s=s+$NF}END{print s}' ${dir1}/${i}_${n}.peak.TE.bed` 
        awk '{a[$8]=a[$8]+$9}END{for(i in a)print i,a[i]/"'$a'"}' ${dir1}/${i}_${n}.peak.TE.bed |\
            sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 $TEinfo) | awk '{print $1,$2,$3,$2/$3}' > ${i}_${n}.TE.temp 
        b=`cut -f 4 ${dir1}/${i}_${n}.summit.TE.bed | sort -u | wc -l` 
        awk '{a[$8]=a[$8]+1}END{for(i in a)print i,a[i],"'$b'"}' ${dir1}/${i}_${n}.summit.TE.bed |\
            sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 ${i}_${n}.TE.temp) |\
            awk '{print "'$i' '$n' "$0}' >> enrichment.TE.txt 
    done 
    ##### dTE:length 
    for n in balance unbalance 
    do 
        a=`awk '{s=s+$NF}END{print s}' ${dir2}/${i}.${n}.peak.dTE.bed` 
        awk '{a[$8]=a[$8]+$9}END{for(i in a)print i,a[i]/"'$a'"}' ${dir2}/${i}.${n}.peak.dTE.bed |\
            sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 $TEinfo) | awk '{print $1,$2,$3,$2/$3}' > ${i}_${n}.dTE.temp 
        b=`cut -f 4 ${dir2}/${i}.${n}.summit.dTE.bed | sort -u | wc -l` 
        awk '{a[$8]=a[$8]+1}END{for(i in a)print i,a[i],"'$b'"}' ${dir2}/${i}.${n}.summit.dTE.bed |\
            sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 ${i}_${n}.dTE.temp) |\
            awk '{print "'$i' '$n' "$0}' >> enrichment.dTE.txt 
    done 
done 
rm *temp 

#######################################################
###### 8. DH divergence vs. expression divergence #####
#######################################################
pro=~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genetssupd5k.strand.longestcds.bed
trip=~/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_7species_wheat4AB/wheat6.triplet.info
expr=~/yuyun/dap-seq/01_RNA-seq_CS/01_DEgenes/wheat6.CK.fpkm

##### expression divergence
sort -k1,1 $trip | join -1 1 - -2 1 <(sort -k1,1 $expr) | sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 $expr) | \
    sort -k3,3 | join -1 3 - -2 1 <(sort -k1,1 $expr) | awk '{print $3"-"$2"-"$1,$4,$5,$6}' > expr.1-1-1.fpkm.txt
awk '$2>1||$3>1||$4>1' expr.1-1-1.fpkm.txt | \
    awk '{print $1,$2/($2+$3+$4),$3/($2+$3+$4),$4/($2+$3+$4)}' > expr.1-1-1.fpkm.filter.norm.txt
Rscript get_distance.r expr.1-1-1.fpkm.filter.norm.txt

##### TE-embedded DH divergence
mkdir tmp/
x=0 
i=DH
cat ../02_dTE_blast_score/blast_fasta/DHS.*balance.summit.dTE.bed ../01_overlap_TE_score/DHS_*balance.summit.TE.bed |\
    cut -f 1-4 | sort -u | bedtools sort -i - > DHS.summit.TE+dTE.bed 
summit=DHS.summit.TE+dTE.bed 
peak=cs_${i}.peak.bed
bed=${i}.bed 
stats=${i}.bam.stats  
m=`grep 'reads mapped:' $stats | awk '{print $4}'` 

bedtools intersect -a $summit -b $pro -wa -wb | awk '{print $8,$4,($7-$6+1)/2+$6-$2}' | awk '{sub("-","",$3);print}' > tmp/${i}.dis 
sort -k2,2 tmp/${i}.dis | join -1 2 - -2 4 <(sort -k4,4 $peak) | \
    awk '{print $4,$5,$6,$1,$2,$3}' | tr ' ' '\t' | bedtools coverage -a - -b $bed -counts | tee tmp/${i}.count |\
    awk '{print $5,(($7*1000000000)/(($3-$2)*"'$m'"))*exp(-($6/5000))}' > tmp/${i}.tmp  

awk '{a[$1]=a[$1]+$2}END{for(i in a)print i,a[i]}' tmp/${i}.tmp | sort -k1,1 | \
    join -1 1 - -2 1 <(cut -f 4 $pro | sort ) -a 2 | awk 'BEGIN{print "atf '$i'"}{if(NF==1)print $1,"0";else print $0}' | \
    sort -k1,1 > ${i}.score

score=${i}.score 
cut -f 2- -d ' ' ${trip} | sort -k1,1 | join -1 1 - -2 1 <(sort -k1,1 $score) |\
    sort -k2,2 | join -1 2 - -2 1 <(sort -k1,1 $score) | sort -k3,3 |\
    join -1 3 - -2 1 <(sort -k1,1 $score) | awk '{print $3,$2,$1,$4,$5,$6}' > ${i}.score.1-1-1.txt 
awk '$4>0.5 || $5>0.5 || $6>0.5' ${i}.score.1-1-1.txt > ${i}.score.1-1-1.filter.txt
n=`cat ${i}.score.1-1-1.filter.txt | wc -l`
awk '{print $1"-"$2"-"$3,$4/($4+$5+$6),$5/($4+$5+$6),$6/($4+$5+$6)}' ${i}.score.1-1-1.filter.txt > ${i}.score.1-1-1.filter.norm.txt
Rscript get_distance.r DHS.score.1-1-1.filter.norm.txt

cut -f 1,5 -d ' ' DHS.score.1-1-1.filter.norm.txt.dist > DHS.score.dist 
Rscript work_get_quantile.r DHS.score.dist 3 up 

###### DH divergence vs. expression divergence
rm plot.distance-box.txt
for x in {1..3}
do 
    n1=`awk 'NR=="'$x'"{print $2*1}' DHS.score.dist.quantile`
    n2=`awk 'NR=="'$x'"{print $3*1}' DHS.score.dist.quantile`

    sort -k1,1 DHS.score.dist | join -1 1 - -2 1 <(sort -k1,1 expr.1-1-1.fpkm.filter.norm.txt.dist | cut -f 1,5 -d ' ') | \
        awk -v n1=$n1 -v n2=$n2 '$2>=n1&&$2<=n2{print $0" level'$x'"}' >> plot.distance-box.txt 
done 

cut -f 1 DHS.score.1-1-1.filter.norm.txt.dist -d ' ' | cat - <(cut -f 1 expr.1-1-1.fpkm.filter.norm.txt.dist -d ' ') | sort | uniq -d | fgrep -f - DHS.score.1-1-1.filter.norm.txt.dist > plot.DHS.trip.txt
sort -k1,1 DHS.score.1-1-1.filter.norm.txt.dist | cut -f 1,5 -d ' ' | join -1 1 <(sort -k1,1 expr.1-1-1.fpkm.filter.norm.txt.dist | cut -f 1-4 -d ' ') -2 1 - | sort -k5,5nr  > plot.expr.trip.txt
