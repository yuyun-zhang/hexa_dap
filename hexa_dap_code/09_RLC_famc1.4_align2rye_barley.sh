wheatfl=/public/home/zhangyuyun/yuyun/genome_data/TE_annotation/LTRharvest/CS_TE_LTR.fl.long_terminal_repeat.bed
ryefl=/public/home/zhangyuyun/yuyun/genome_data/TE_annotation/LTRharvest/rye_TE_LTR.fl.long_terminal_repeat.bed
barleyfl=/public/home/zhangyuyun/yuyun/genome_data/TE_annotation/LTRharvest/Barley_morexv2_TE_LTR.fl.long_terminal_repeat.bed

wheat=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
rye=~/yuyun/genome_data/rye/Weining.v1.0.fa
barley=~/yuyun/genome_data/barley_morex/Barley_Morex_V2_pseudomolecules.fasta

for i in RLC_famc1.4 
do 
    for sub in A B D 
    do  
        awk '$1~"'$sub'"' $wheatfl | grep ${i}"_" | awk '{if(($6=="+"&&$4~/%1/)||($6=="-"&&$4~/%2/))print $0>"'${i}'.wheat'$sub'.5prime.bed";else print $0>"'${i}'.wheat'$sub'.3prime.bed"}'
        sed 's/%[12]//' ${i}.wheat${sub}.5prime.bed | bedtools getfasta -fi $wheat -bed - -s -nameOnly > ${i}.wheat${sub}.5prime.fa 
        sed 's/%[12]//' ${i}.wheat${sub}.3prime.bed | bedtools getfasta -fi $wheat -bed - -s -nameOnly > ${i}.wheat${sub}.3prime.fa 
        seqkit concat -j 10 --quiet -w 0  ${i}.wheat${sub}.5prime.fa ${i}.wheat${sub}.3prime.fa > ${i}.wheat${sub}.5+3.fa 
    done 
    grep ${i}"_" $ryefl | awk '{if(($6=="+"&&$4~/%1/)||($6=="-"&&$4~/%2/))print $0>"'${i}'.rye.5prime.bed";else print $0>"'${i}'.rye.3prime.bed"}'
    sed 's/%[12]//' ${i}.rye.5prime.bed | bedtools getfasta -fi $rye -bed - -s -nameOnly > ${i}.rye.5prime.fa 
    sed 's/%[12]//' ${i}.rye.3prime.bed | bedtools getfasta -fi $rye -bed - -s -nameOnly > ${i}.rye.3prime.fa 
    seqkit concat -j 10 --quiet -w 0  ${i}.rye.5prime.fa ${i}.rye.3prime.fa > ${i}.rye.5+3.fa 

    grep ${i}"_" $barleyfl | awk '{if(($6=="+"&&$4~/%1/)||($6=="-"&&$4~/%2/))print $0>"'${i}'.barley.5prime.bed";else print $0>"'${i}'.barley.3prime.bed"}'
    sed 's/%[12]//' ${i}.barley.5prime.bed | bedtools getfasta -fi $barley -bed - -s -nameOnly > ${i}.barley.5prime.fa 
    sed 's/%[12]//' ${i}.barley.3prime.bed | bedtools getfasta -fi $barley -bed - -s -nameOnly > ${i}.barley.3prime.fa 
    seqkit concat -j 10 --quiet -w 0  ${i}.barley.5prime.fa ${i}.barley.3prime.fa > ${i}.barley.5+3.fa 
done 
rm *prime.fa *prime.bed 

###### align 
rm align.random.number 
for i in RLC_famc1.4 
do 
    for j in wheatA wheatB wheatD rye barley 
    do 
        n=`grep '>' ${i}.${j}.5+3.fa | wc -l`
        m=`grep '>' ${i}.wheatA.5+3.fa | wc -l`
        x=`echo ${n}*500/$m | bc`
        echo $i $j $x >> align.random.number 
        cat ${i}.${j}.5+3.fa | paste - - | shuf -n $x | sed "s/\t/\n/;s/>${i}_/>${j}%/;s/(.*)//" > ${i}.${j}.random.fa
    done 

    for line in wheatB-barley wheatA-barley wheatD-barley wheatA-wheatA wheatB-wheatB wheatD-wheatD wheatA-wheatD  wheatA-wheatB wheatB-wheatD wheatA-rye wheatB-rye wheatD-rye  
    do 
    {
        sp1=`basename $line | sed 's/-.*//'`
        sp2=`basename $line | sed 's/.*-//'`
        p=`cat ${i}.${sp2}.random.fa | wc -l`

        cat ${i}.$sp1.random.fa ${i}.$sp2.random.fa > ${i}.${sp1}-${sp2}.fa 
        if [ ! -s ${i}.${sp1}-${sp2}.mafft.fa ]
        then 
        echo start $i $sp1 $sp2 
        mafft --thread 15 --quiet --adjustdirection ${i}.${sp1}-${sp2}.fa > ${i}.${sp1}-${sp2}.mafft.fa
        distmat ${i}.${sp1}-${sp2}.mafft.fa -nucmethod 2 ${i}.${sp1}-${sp2}.mafft.distmat.out 
        tail -n +9 ${i}.${sp1}-${sp2}.mafft.distmat.out | cut -f 2- |\
            sed 's/\t\t[0-9]*-[0-9]* [0-9]*//' | sed 's/  //' > ${i}.${sp1}-${sp2}.mafft.distmat.max
        Rscript work_get_dist.r ${i}.${sp1}-${sp2}.mafft.distmat.max
        n1=`grep $i align.random.number | grep $sp1 | awk '{print $3}'`
        sed 's/V//;s/V//' ${i}.${sp1}-${sp2}.mafft.distmat.max.table | \
            awk '$1<="'$n1'"*1&&$2>"'$n1'"*1{print "'$i' '$sp1'-'$sp2'",$3}' > plot.${i}.${sp1}-${sp2}.align.txt   
        echo finished $i $sp1 $sp2 
        fi 
    }&
    done
    wait 
done 
wait
echo finished alignment 

cat plot.*.align.txt | grep -v % | cat - <(echo -e "RLC_famc7.3 wheatA-barley 49\nRLC_famc7.3 wheatB-barley 49\nRLC_famc7.3 wheatD-barley 49") > plot.align.txt 

##################################
#### Rscript: work_get_dist.r ####
##################################
file<-commandArgs(T)

max<-as.matrix(read.table(file,sep = "\t"))
rownames(max)<-paste0("V",1:nrow(max))
table<-na.omit(reshape2::melt(max))
write.table(table,paste0(file,".table"),quote=F,row.names=F,col.names=F)


