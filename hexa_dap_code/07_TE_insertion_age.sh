###################################################
########## Definition of full-length TE ###########
###################################################
i=cs 
grep "RL" ${i}_TE.bed > ${i}_TE_LTR.bed    
bedtools getfasta -fi ${i}_genome.fa -bed ${i}_TE_LTR.bed > ${i}_TE_LTR.fa
~/miniconda3/envs/yuyun_py2.7/bin/gt suffixerator -db ${i}_TE_LTR.fa -indexname $i -dna yes -md5 no -suf -lcp
~/miniconda3/envs/yuyun_py2.7/bin/gt ltrharvest -index ${i} -overlaps best -seed 30 -minlenltr 100 -maxlenltr 2000 \
    -mindistltr 3000 -maxdistltr 25000 -similar 85 -mintsd 4 -maxtsd 20 -motif tgca -motifmis 1 -vic 60 -xdrop 5 \
    -mat 2 -mis -2 -ins -3 -del -3 -out ${i}.lh.out -gff3 ${i}.lh.gff3 
echo finished_LTRharvest_${i}

grep '>' ${i}.lh.out | sed 's/ .*//;s/>//;s/:/\t/;s/-/\t/' | bedtools intersect -a ${i}_TE_LTR.bed -b - -wa > ${i}_TE_LTR.fl.bed 
grep '>' ${i}.lh.out | sed 's/ \[.*//;s/>//;s/(dbseq-nr /seq/;s/)//' | paste - <(grep repeat_region ${i}.lh.gff3 -w) | \
    awk '{print $1"\t"$6"\t"$7}' | sed 's/:/\t/;s/-/\t/' | awk '{print $1"\t"$2+$4"\t"$2+$5}' | bedtools intersect -a - -b ${i}_TE_LTR.bed -wa -wb |\
    cut -f 1,2,3,7- > ${i}_TE_LTR.fl.repeat_region.bed
grep '>' ${i}.lh.out | sed 's/ \[.*//;s/>//;s/(dbseq-nr /seq/;s/)//' | sort -k2,2 | \
    join -1 2 - -2 1 <(grep long_terminal_repeat ${i}.lh.gff3 -w | sort -k1,1) | awk '{print $2,$5,$6,$1}' | \
    sed 's/:/\t/;s/-/\t/' | awk '{print $1"\t"$2+$4"\t"$2+$5"\t"$6}' | \
    bedtools intersect -a - -b ${i}_TE_LTR.bed -wa -wb | \
    awk '{if(NR%2==1&&$10=="+")print $1,$2,$3,$8"_LTR"(NR+1)/2"%1",$9,$10;else if(NR%2!=1&&$10=="+") print $1,$2,$3,$8"_LTR"NR/2"%2",$9,$10;\
            else if (NR%2==1&&$10=="-")print $1,$2,$3,$8"_LTR"(NR+1)/2"%2",$9,$10;else print $1,$2,$3,$8"_LTR"NR/2"%1",$9,$10}' |\
    tr ' ' '\t' > ${i}_TE_LTR.fl.long_terminal_repeat.bed

###################################################
######### LTR distance of full-length TE ##########
###################################################
flLTR=../../cs_TE_LTR.fl.long_terminal_repeat.bed
flTE=../../cs_TE_LTR.fl.repeat_region.bed

n=`cat ../../cs_TE_LTR.fl.repeat_region.bed | wc -l` 
mkdir tmp
x=0
for i in `seq $n`
do 
x=`expr $x + 1`
{
    awk 'NR=="'$i'"*2-1||NR=="'$i'"*2' ${flLTR} | bedtools getfasta -fi $genome -bed - -nameOnly > tmp/${i}.fa 
    TE=`awk 'NR=="'$i'"*2-1' ${flLTR} | cut -f 4 | sed 's/%.*//'`
    muscle -quiet -in tmp/${i}.fa -out tmp/${i}.mus
    distmat tmp/${i}.mus -nucmethod 2 -outfile tmp/${i}.dist
    dist=`tail -n 2 tmp/${i}.dist | head -n 1 | awk '{print $2}'`
    echo -e "${TE}""\t""$dist" > tmp/${i}.dist.txt
    echo finished_${i}
}&
if [ $x == 20 ]
then 
    x=0
    wait 
fi 
done 
wait 

find ./tmp/ -name "*.dist.txt" -exec cat '{}' ';' > dist.out
grep '%1' ${flLTR} | bedtools intersect -a - -b ${flTE} -wa -wb -f 1 | awk '{print $7,$8,$9,$4}' | \
    sed 's/%1//;s/ /\t/g' | sort -k4,4 | join -1 4 - -2 1 <(sort -k1,1 dist.out) | \
    awk '{print $2,$3,$4,$1,$5}' | tr ' ' '\t' | bedtools sort -i - > ${flTE/bed/}LTRdist.bed 



###################################################
############## TE Phylogenetic Trees ##############
###################################################
genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
mp=mergeTFBS_PE_summits.filter.p10.bed

#### full length 
TE=~/yuyun/genome_data/TE_annotation/LTRharvest/CS_TE_LTR.fl.repeat_region.bed
for i in RLG_famc13 RLC_famc1.4 RLG_famc7.3 
do 
{
    grep $i -w $TE | bedtools intersect -a - -b merge_peaks_w_motifs.bed -wa | sort -u | bedtools sort -i - |\
        awk '{if($1~/A/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_A_peak";\
              else if($1~/B/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_B_peak";\
              else if($1~/D/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_D_peak"}' > ${i}.peak.fl.bed 
    grep $i -w $TE | bedtools intersect -a - -b $mp -wa -v | sort -u | bedtools sort -i - |\
        awk '{if($1~/A/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_A_nonpeak";\
              else if($1~/B/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_B_nonpeak";\
              else if($1~/D/)print $1"\t"$2"\t"$3"\t"$4"_"NR"_D_nonpeak"}' > ${i}.nonpeak.fl.bed 

    shuf -n 200 ${i}.peak.fl.bed | bedtools getfasta -fi $genome -bed - -nameOnly > ${i}.peak.fl.fa
    shuf -n 200 ${i}.nonpeak.fl.bed | bedtools getfasta -fi $genome -bed - -nameOnly > ${i}.nonpeak.fl.fa
    grep '>' ${i}.peak.fl.fa | sed 's/>//' | awk '{print $1,$1}' | sed 's/ .*_[0-9]*_/ /' > ${i}.fl.anno
    grep '>' ${i}.nonpeak.fl.fa | sed 's/>//' | awk '{print $1,$1}' | sed 's/ .*_[0-9]*_/ /'  >> ${i}.fl.anno
    grep '>' ${i}.peak.fl.fa | sed 's/>//' | awk '{print $1,$1}' | sed 's/ .*_[0-9]*_/ /;s/_peak$//' > ${i}.fl.anno2
    grep '>' ${i}.nonpeak.fl.fa | sed 's/>//' | awk '{print $1,$1}' | sed 's/ .*_[0-9]*_/ /;s/_nonpeak$//'  >> ${i}.fl.anno2

    cat ${i}.*.fl.fa > ${i}.fl.fa 
    mafft --thread 10 --quiet ${i}.fl.fa > ${i}.mafft.fl.fa
    fasttree -nt -quiet ${i}.mafft.fl.fa > ${i}.mafft.fl.tree
    echo finished_${i} 
}&
wait 
awk '{if($2=="A_peak")print $1,"1 0 0 4";else if($2=="B_peak")print $1,"0 2 0 4";else if($2=="D_peak")print $1,"0 0 3 4";else if($2=="A_nonpeak")print $1,"1 0 0 0";else if($2=="B_nonpeak")print $1,"0 2 0 0";else if($2=="D_nonpeak")print $1,"0 0 3 0"}' ${i}.fl.anno > ${i}.fl.heatmap 
done 

