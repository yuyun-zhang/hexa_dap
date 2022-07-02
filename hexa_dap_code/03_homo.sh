###################################################
######## 1. subgenome synteny regions #############
###################################################
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/catalog.py ortholog A B --cpus=16
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py screen --minspan=10 --simple A.B.lifted.anchors A.B.lifted.anchors.new
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py simple A.B.lifted.anchors.new --qbed=A.bed --sbed=B.bed --coords --noheader #A.B.lifted.anchors.simple

python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/catalog.py ortholog A D --cpus=16
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py screen --minspan=10 --simple A.D.lifted.anchors A.D.lifted.anchors.new
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py simple A.D.lifted.anchors.new --qbed=A.bed --sbed=D.bed --coords --noheader

python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/catalog.py ortholog B D --cpus=16
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py screen --minspan=10 --simple B.D.lifted.anchors B.D.lifted.anchors.new
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/synteny.py simple B.D.lifted.anchors.new --qbed=B.bed --sbed=D.bed --coords --noheader

for i in A.B A.D B.D
do
    awk '{print $2"\t"$3"\t"$4"\t"$1}' ${i}.lifted.anchors.simple | bedtools sort -i - > ${i}.lifted.anchors.simple.bed
done
cat *.lifted.anchors.simple.bed | bedtools sort -i - > cs_ABD.synteny.bed

###################################################
######## 2. subgenome homologous regions ##########
###################################################
nucmer -t 3 --mum -p ${ref1}_vs_${que1} $ref  $que

rm samechr.filelist
for i in {1..7};do ls chrABD_${i}[ABD]_vs_chrABD_${i}[ABD].delta >> samechr.filelist ;done

genome_dir=/public/home/zhangyuyun/yuyun/dap-seq/07_conservation_score/00_genome/
synteny_dir=/public/home/zhangyuyun/yuyun/dap-seq/08_MCScanX/02_ABD

x=0
for i in `cat samechr.filelist`
do
x=`expr $x + 1`
{
    sub1=`basename $i | sed 's/chrABD_[0-9]//;s/_vs.*//'`
    sub2=`basename $i | sed 's/.delta//;s/.*_[0-9]//'`
    chr1=`basename $i | sed 's/chrABD_//;s/[ABD]_vs.*//'`
    chr2=`basename $i | sed 's/.*chrABD_//;s/[ABD]\.delta//'`
    /public/home/xieyilin/software/mummer-4.0.0beta2/show-coords -T ${i} > ${i/.delta/}.coords
    cat ${i/.delta/}.coords | tail -n +5 |awk 'BEGIN{OFS="\t"}{print $8,$1,$2,$9,$3,$4}' | \
        awk 'BEGIN{OFS="\t"}{if($2>$3)print $1,$3,$2,$4,$5,$6;else print $0}' |\
        awk 'BEGIN{OFS="\t"}{if($5>$6)print $1,$2,$3,$4,$6,$5;else print $0}' > ${i/.delta/}.bedpe
    awk '$3-$2>400&&$6-$5>400' ${i/.delta/}.bedpe > ${i/.delta/}.bedpe.filtertmp
    bedtools intersect -a ${i/.delta/}.bedpe.filtertmp -b ${genome_dir}/cs_marked.no-n.bed -wa -f 0.9 |\
        awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' | bedtools intersect -a - -b ${genome_dir}/cs_marked.no-n.bed -wa -f 0.9 |\
        awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' | sort -u | bedtools sort -i - |\
        sed 's/:1-[0-9]*//g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t'${chr1}$sub1'.vs.'${chr2}${sub2}'_"NR}' |\
        tee ${i/.delta/}.filtern.bedpe | awk '{print $1"\t"$2"\t"$3"\t"$7"\n"$4"\t"$5"\t"$6"\t"$7}' > ${i/.delta/}.filtern.bed
    bedtools intersect -a ${i/.delta/}.filtern.bedpe -b ${synteny_dir}/cs_ABD.synteny.bed -wa -f 0.9 -wb | \
        awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"\t"$11}' | \
        bedtools intersect -a - -b ${synteny_dir}/cs_ABD.synteny.bed -wa -wb -f 0.9 | \
        awk '$8==$12{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"-syn"}' | \
        sort -u | awk '{print $1"\t"$2"\t"$3"\t"$7"\n"$4"\t"$5"\t"$6"\t"$7}' | \
        bedtools sort -i - > ${i/.delta/}.filtern.synteny.bed
    cut -f 4 ${i/.delta/}.filtern.synteny.bed | sed 's/-syn//' | sort -u | fgrep -f - ${i/.delta/}.filtern.bed -v -w |\
        awk '{print $0"-nonsyn"}' > ${i/.delta/}.filtern.nonsynteny.bed
    cat ${i/.delta/}.filtern.synteny.bed ${i/.delta/}.filtern.nonsynteny.bed | bedtools sort -i - > ${i/.delta/}.filtern.bed
    rm ${i/.delta/}.bedpe.filtertmp ${i/.delta/}.coords ${i/.delta/}.filtern.bedpe ${i/.delta/}.filtern.nonsynteny.bed ${i/.delta/}.filtern.synteny.bed
    echo finished_${i}
}&
    if [ $x == 10 ]
    then
        wait
        x=0
    fi
done
wait

for i in {1..7}
do
{
    bedtools intersect -a chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -b chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -wa -wb | tee chrABD_${i}A_vs_chrABD_${i}BD.filtern.RH.temp | \
        grep -v nonsyn | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}A_vs_chrABD_${i}BD.filtern.RH.syn.bed
    grep nonsyn chrABD_${i}A_vs_chrABD_${i}BD.filtern.RH.temp | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}A_vs_chrABD_${i}BD.filtern.RH.nonsyn.bed
    bedtools intersect -a chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -b chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -v |\
        cat - <(bedtools intersect -b chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -a chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -v) |\
        cut -f 1-3 | awk '$1~"A"' | sort-bed - | mergeBed -i - > chrABD_${i}A_vs_chrABD_${i}BD.filtern.1sub.bed
    bedtools intersect -a chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -b chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -wa -wb | tee chrABD_${i}B_vs_chrABD_${i}AD.filtern.RH.temp | \
        grep -v nonsyn | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}B_vs_chrABD_${i}AD.filtern.RH.syn.bed
    grep nonsyn chrABD_${i}B_vs_chrABD_${i}AD.filtern.RH.temp | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}B_vs_chrABD_${i}AD.filtern.RH.nonsyn.bed
    bedtools intersect -a chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -b chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -v |\
        cat - <(bedtools intersect -b chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -a chrABD_${i}A_vs_chrABD_${i}B.filtern.RH.bed -v) |\
        cut -f 1-3 | awk '$1~"B"' | sort-bed - | mergeBed -i - > chrABD_${i}B_vs_chrABD_${i}AD.filtern.1sub.bed
    bedtools intersect -a chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -b chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -wa -wb | tee chrABD_${i}D_vs_chrABD_${i}AB.filtern.RH.temp | \
        grep -v nonsyn | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}D_vs_chrABD_${i}AB.filtern.RH.syn.bed
    grep nonsyn chrABD_${i}D_vs_chrABD_${i}AB.filtern.RH.temp | cut -f 1-3 | sort-bed - | mergeBed -i - | \
        awk '$3-$2>400' > chrABD_${i}D_vs_chrABD_${i}AB.filtern.RH.nonsyn.bed
    bedtools intersect -a chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -b chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -v |\
        cat - <(bedtools intersect -b chrABD_${i}A_vs_chrABD_${i}D.filtern.RH.bed -a chrABD_${i}B_vs_chrABD_${i}D.filtern.RH.bed -v) |\
        cut -f 1-3 | awk '$1~"D"' | sort-bed - | mergeBed -i - > chrABD_${i}D_vs_chrABD_${i}AB.filtern.1sub.bed
}&
done
wait
cat `ls chrABD_[1-7]A_vs_chrABD_[1-7]BD.filtern.RH.syn.bed chrABD_[1-7]B_vs_chrABD_[1-7]AD.filtern.RH.syn.bed chrABD_[1-7]D_vs_chrABD_[1-7]AB.filtern.RH.syn.bed` | \
    sort-bed - | awk '{print $0"\thomo3-syn_"NR}' > cs_ABDhomo.multiple.RH.homo3-syn.bed
cat `ls chrABD_[1-7]A_vs_chrABD_[1-7]BD.filtern.RH.nonsyn.bed chrABD_[1-7]B_vs_chrABD_[1-7]AD.filtern.RH.nonsyn.bed chrABD_[1-7]D_vs_chrABD_[1-7]AB.filtern.RH.nonsyn.bed` | \
    sort-bed - | bedtools subtract -a - -b cs_ABDhomo.multiple.RH.homo3-syn.bed | \
    awk '$3-$2>400{print $0"\thomo3-nonsyn_"NR}' > cs_ABDhomo.multiple.RH.homo3-nonsyn.bed
cat `ls chrABD_[1-7]A_vs_chrABD_[1-7]BD.filtern.1sub.bed chrABD_[1-7]B_vs_chrABD_[1-7]AD.filtern.1sub.bed chrABD_[1-7]D_vs_chrABD_[1-7]AB.filtern.1sub.bed` | \
    sort -u | sort-bed - | bedtools subtract -a - -b <(cat cs_ABDhomo.multiple.RH.homo3-syn.bed cs_ABDhomo.multiple.RH.homo3-nonsyn.bed) |\
    awk '$3-$2>400{print $0"\thomo2_"NR}' > cs_ABDhomo.multiple.RH.homo2.bed

cat cs_ABDhomo.multiple.RH.homo2.bed cs_ABDhomo.multiple.RH.homo3-syn.bed cs_ABDhomo.multiple.RH.homo3-nonsyn.bed |\
    bedtools subtract -a ~/yuyun/genome_data/wheat6_bed/IWGSC_v1.1_HC_20170706_parts.genome.bed -b - | grep -v Un  > cs_ABDhomo.multiple.RH.nohomo.bed

cat cs_ABDhomo.multiple.RH.homo3-syn.bed cs_ABDhomo.multiple.RH.homo3-nonsyn.bed | bedtools sort -i - | awk '{print $1"\t"$2"\t"$3"\thomo3"}' > cs_ABDhomo.multiple.RH.homo3.bed
