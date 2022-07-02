TEinfo=~/yuyun/genome_data/TE_annotation/clariTE/wheat6_TE.ingenome.txt #### the ratio of length of each TE subfamily in genome 
TE=~/yuyun/genome_data/cs_TE_split.bed
rm enrichment-length.TE.txt 
for i in mergeTFBS 
do 
    for sub in A B D 
    do 
        for j in homo nohomo 
        do             
            ####### overlap length enrichment 
            TEsummit=../04_TE_ratio/homo_bed_DHS/${i}.${sub}.summits.${j}.DHS.TE.bedpe ##### TE-embedded mergeTFBS with DHS
            n=`cut -f 4 $TEsummit | sort -u | wc -l`

            cut -f 4 ../04_TE_ratio/homo_bed_DHS/${i}.${sub}.summits.${j}.DHS.bed | \
                fgrep -f - ../01_stats/homo_bed/${i}.${sub}.peaks.${j}.bed -w | bedtools sort -i - | tee ../04_TE_ratio/homo_bed_DHS/${i}.${sub}.peaks.${j}.DHS.bed |\
                bedtools intersect -a - -b $TE -wo > ../04_TE_ratio/homo_bed_DHS/${i}.${sub}.peaks.${j}.DHS.TE.bedpe
            TEpeak=../04_TE_ratio/homo_bed_DHS/${i}.${sub}.peaks.${j}.DHS.TE.bedpe
            a=`awk '{s=s+$3-$2}END{print s}' $TEpeak`
            cut -f 4,8 $TEsummit | sort -u | cut -f 2 | sort | uniq -c | awk '{print $2,$1,"'$n'"}' | sort -k1,1 |\
                join -1 1 - -2 1 <(awk '{a[$8]=a[$8]+$NF}END{for(i in a)print i,a[i],"'$a'"}' $TEpeak | sort -k1,1)| \
                join -1 1 - -2 1 <(sort -k1,1 $TEinfo) | awk '{print "'$i' '$sub' '$j'",$0,($4/$5)/$6}' >> enrichment-length.TE.txt 
        done 
    done 
done 






