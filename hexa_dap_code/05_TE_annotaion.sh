genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
output=cs

######## 1. repeatmasker
x=0
for i in chr{1..7}A chr{1..7}B chr{1..7}D
do 
x=`expr $x + 1`
{
    echo ${i} | ~/miniconda3/envs/yuyun/bin/seqtk subseq ${genome} - > ${i}.fa
    ~/miniconda3/envs/yuyun/bin/RepeatMasker -pa 20 ${i}.fa -gff -nolow -xm -lib ~/yuyun/software/CLARI-TE/clariTeRep.fna -e ncbi
    echo "RepeatMasker ${i}_finished!"
}&
if [ $x == 4 ]
then 
    x=0
    wait 
fi 
done
wait 

######## 2. get required file format for clari_te
genomesize=$genome.fai
x=0
for i in chr{1..7}A chr{1..7}B chr{1..7}D
do
x=`expr $x + 1`
{
    ### get gene annotation in embl format
    echo "get embl"
    ~/miniconda3/envs/yuyun/bin/seqret -sequence ${i}.fa -outseq ${i}.embl -osformat embl
    embl_file=${i}_anno.embl
    head -n 2 ${i}.embl >${i}_anno.embl
    echo -e "AC   unknown"";""\n""XX""\n""XX""\n""FH   Key            Location/Qualifiers""\n""FH" >> $embl_file
    count=0
    repeatmasker_gff=${i}.fa.out.gff
    grep ${i} $repeatmasker_gff | while read line 
    do
        count=$(( $count + 1 ))
        location=`echo $line |awk 'BEGIN{OFS=".."}{print$4,$5}'`
        locustag=`echo $line |awk -v count=$count 'BEGIN{OFS="_"}{print$1,"REPEATMASKER",count}'`
        id=`echo $line |awk -v count=$count 'BEGIN{OFS="_"}{print$1,"REPEATMASKER",$4,$5,"Repeat_region",count}'`
        echo -e "FT   repeat_region   $location""\n""FT                   /locus_tag=\"$locustag\"""\n""FT                   /id=\"$id\"" >> $embl_file
    done
    echo "XX" >> $embl_file
    tail -n +3 ${i}.embl >> $embl_file
    # rm -f ${i}.embl
    ###repeatmasker
    echo "get xm_file"
    xm_file=${i}.fa.out.xm
    grep ${i} $xm_file | grep -v "Simple_repeat" |sed 's/Unspecified/Unknown/g' >${i}.fas.out.xm

    ###CLARI_TE
    echo "clari_te"
    perl ~/yuyun/software/CLARI-TE/clariTE.pl -LTR ~/yuyun/software/CLARI-TE/clariTeRep_LTRposition.txt -gene $embl_file -classi ~/yuyun/software/CLARI-TE/clariTeRep_classification.txt -fasta ${i}.fa -dir ./ ${i}.fas.out.xm
}&
if [ $x == 7 ]
then 
    x=0
    wait 
fi 
done
wait 

######## 3. embl2gff 
for i in chr{1..7}A chr{1..7}B chr{1..7}D
do
    seqret -sformat embl -sequence $i.fas_annoTE.embl  -feature -osformat fasta -offormat gff -auto -osname2 $i
    awk '$3!="repeat_region"' $i.gff |cut -f 9 |tr ';' '\t'|awk -F "\t" '{for(i=1;i<=NF;i++)if($i~/post/)print $1"\t"$i}'|sed 's/ID=//'|awk '{print $1"\t"$3}' >  $i.merge.id 
    awk '$3=="repeat_region"' $i.gff|tr ';' '\t'|awk -F "\t"  -v chr=$i 'BEGIN{OFS="\t"}{for(i=9;i<=NF;i++)if($i~/post/)print  chr,$4-1,$5,$7,$i}'|tr ' ' '\t'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,".",$4}' >  $i.TE.1.bed
    awk '$9~/^Parent/' $i.gff|awk -v chr=$i 'BEGIN{OFS="\t"}{print chr,$4-1,$5,$9,".",$7}'|sed 's/Parent=//' > $i.merge.bed
    Rscript f_merge.r $i
    cat $i.TE.1.bed $i.merge.1.bed > $i.TE.bed
done

cat *[1-7]*TE.bed > /public/home/zhangyuyun/yuyun/genome_data/TE_annotation/clariTE/${output}_TE_clariTE.bed 

rm *bed *embl *fa *cat.gz *masked *out *xm *fasta *tbl *id 
mkdir ${output}_repeatmasker_out && mv *gff ${output}_repeatmasker_out/
