function getmotif(){
    i=$1
    jaspar=/opt/yuyun/03_dap-seq/meme/JASPAR2020_CORE_plants_non-redundant_pfms.meme
    meme-chip -meme-minw 5 -meme-maxw 12 -meme-nmotifs 3 -meme-p 8 -dreme-m 0 -spamo-skip -fimo-skip -norc -oc meme-chip/${file/.top6000.summitupd300.fa/} ${file} -db ${jaspar}

    genome=/mnt/zhaofei/wheat_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
    genomesize=/mnt/zhaofei/wheat_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.size
    mkdir ame_out-${date}
    mkdir ame_out-${date}/${i/.top6000.summitupd300.fa/}_jaspar
    ame --control --shuffle-- -oc ame_out-${date}/${i/.top6000.summitupd300.fa/}_jaspar $i ${jaspar}
}

j=0
for file in `ls *.top6000.summitupd300.fa`
do
    j=`expr $j + 1`
{
    getmotif ${file}
    echo "finished $file"
}&
    if [ $j == 4 ]
    then
        echo 'wait'
        j=0
        wait
    fi
done

wait
