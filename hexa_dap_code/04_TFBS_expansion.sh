###################################################
################### fig. 3d #######################
###################################################
python /public/home/zhangyuyun/yuyun/software/jcvi/jcvi/jcvi/compara/catalog.py ortholog A B --cpus=16
genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
mkdir blast_out_bHLH-1A-1/

for sub in A B D
do
{
    nohomo=../02_ABD_seq_homo/01_stats/homo_bed/bHLH-1A-1.${sub}.summits.nohomo.bed
    bedtools slop -g $genomesize -i $nohomo -b 300 | awk '{print $1"\t"$2"\t"$3"\t"NR}' |\
        bedtools getfasta -fi $genome -bed - -name > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa
    makeblastdb -in blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa -out blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa -dbtype nucl
    blastn -query blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa -db blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa -outfmt 6 \
        -out blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.out6 -num_threads 10 -evalue 1e-30
    awk '$1!=$2&&$3>70&&($4-$5-$6)/600>0.7' blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.out6 > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.filter.out6
    awk '{print $1,$2}' blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.filter.out6 | sed 's/::/ /g' | \
        awk '{if($1>$3)print $4,$2;else print $2,$4}' | sort -u > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.pair
    sed 's/ /\n/' blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.pair | sort -u | tee blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.ortho.list | \
        sed 's/:/\t/;s/-/\t/' > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.ortho.bed
    grep '>' blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa | sed 's/>[0-9]*:://' | cat - blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.ortho.list |\
        sort | uniq -u | sed 's/:/\t/;s/-/\t/' > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.non-ortho.bed
    rm blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.ortho.list
    echo finished_wheat_${sub}
}&
done

### arab
genome=~/yuyun/genome_data/arab_allgenome.fa
genomesize=~/yuyun/genome_data/arab_allgenome.genomesize

peak=dist_out_arab/bHLH28_col_a.peaks.bed
bedtools sort -i $peak | bedtools merge -i - | sed 's/chr/Chr/' | awk '{print $1"\t"$2"\t"$3"\t"NR}' |\
    bedtools getfasta -fi $genome -bed - -name > blast_out_bHLH-1A-1/arab_bHLH28.fa
makeblastdb -in blast_out_bHLH-1A-1/arab_bHLH28.fa -out blast_out_bHLH-1A-1/arab_bHLH28.fa -dbtype nucl
blastn -query blast_out_bHLH-1A-1/arab_bHLH28.fa -db blast_out_bHLH-1A-1/arab_bHLH28.fa -outfmt 6 \
    -out blast_out_bHLH-1A-1/arab_bHLH28.out6 -num_threads 10 -evalue 1e-30
awk '{print $1,$0}' blast_out_bHLH-1A-1/arab_bHLH28.out6 | sed 's/:/\t/3;s/-/\t/' | \
    awk '$4!=$5&&$6>70&&($7-$8-$9)/($3-$2)>0.7{print $4,$5}' | sed 's/::/ /g' | \
    awk '{if($1>$3)print $4,$2;else print $2,$4}' | sort -u > blast_out_bHLH-1A-1/arab_bHLH28.pair
sed 's/ /\n/' blast_out_bHLH-1A-1/arab_bHLH28.pair | sort -u | tee blast_out_bHLH-1A-1/arab_bHLH28.ortho.list | \
    sed 's/:/\t/;s/-/\t/' > blast_out_bHLH-1A-1/arab_bHLH28.ortho.bed
grep '>' blast_out_bHLH-1A-1/arab_bHLH28.fa | sed 's/>[0-9]*:://' | cat - blast_out_bHLH-1A-1/arab_bHLH28.ortho.list |\
    sort | uniq -u | sed 's/:/\t/;s/-/\t/' > blast_out_bHLH-1A-1/arab_bHLH28.non-ortho.bed
rm blast_out_bHLH-1A-1/arab_bHLH28.ortho.list
echo finished_arab
wait


###### shuf 400 wheat bHLH-1A-1
for sub in A B D
do
    grep '>' blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.fa | shuf -n 400 | \
        sed 's/[0-9]*:://;s/>//' | sort > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.shuf.list
    fgrep -f blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.shuf.list blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.pair -w > blast_out_bHLH-1A-1/bHLH-1A-1.${sub}.nohomo.shuf.pair
done

grep '>' blast_out_bHLH-1A-1/arab_bHLH28.fa | shuf -n 400 | \
    sed 's/[0-9]*:://;s/>//' | sort > blast_out_bHLH-1A-1/arab_bHLH28.shuf.list
fgrep -f blast_out_bHLH-1A-1/arab_bHLH28.shuf.list blast_out_bHLH-1A-1/arab_bHLH28.pair -w > blast_out_bHLH-1A-1/arab_bHLH28.shuf.pair


###################################################
################### fig. 3e #######################
###################################################
genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
dir=/public/home/zhangyuyun/yuyun/dap-seq/02_TF_analysis/28_CS_DAP/06_summary/
mkdir dist_out/
x=0
for tf in Trihelix-2B-1 `cat ${dir}/HC.tf.list`
do
{
    n=`cat ../01_statsv2/homo_bed/${tf}.allsub.summits.nohomo.bed | wc -l`
    if [ $n -gt 500 ]
    then
        x=`expr $x + 1`
        for sub1 in A B D
        do
        if [ ! -s dist_out/${tf}.${sub1}.dist ];then
            echo start $tf ${sub1}
            summit=../01_statsv2/homo_bed/${tf}.${sub1}.summits.nohomo.bed
            grep -v Un $summit | shuf -n 500 | bedtools sort -i - | bedtools slop -i - -g $genomesize -b 300 |\
                bedtools getfasta -fi $genome -bed - > dist_out/${tf}.${sub1}.nohomo.shuf500.peak.fa

            mafft --adjustdirection --thread 10 --quiet dist_out/${tf}.${sub1}.nohomo.shuf500.peak.fa > dist_out/${tf}.${sub1}.mafft.fa
            distmat dist_out/${tf}.${sub1}.mafft.fa -nucmethod 2 dist_out/${tf}.${sub1}.distmat.out

            tail -n +9 dist_out/${tf}.${sub1}.distmat.out | cut -f 2- |\
                sed 's/\t\t[0-9]*-[0-9]* [0-9]*//' | sed 's/  //' > dist_out/${tf}.${sub1}.distmat.max
            Rscript work_get_dist.r dist_out/${tf}.${sub1}.distmat.max

            awk '$1!=$2{print "'$tf' '$sub1'-'$sub1' "$3}' dist_out/${tf}.${sub1}.distmat.max.table > dist_out/${tf}.${sub1}.dist
            echo "finished_"${tf} ${sub1}
        fi
        done
    fi
}&
if [ $x == 10 ]
then
    x=0
    wait
fi
done
wait
cat dist_out/*.*.dist > plot.intra.nohomo.txt


