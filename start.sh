#/bin/bash/
set -e

trap 'echo "Error occurred at line $LINENO"' ERR

cmake -B build > cmake_log.txt 2>&1
make -C build > build_log.txt 2>&1

# ./data/test.fas

# ./data/Variola.fasta
# data/sars2_similarity_70_1.fas
# data/23s_rRNA.fasta


dataNmae="./data/23s_rRNA.fasta"

alignName="align.fasta"

echo "======================"
echo "======================"

echo " GLS: "

time ./build/GLS -i ${dataNmae} -k 4 > log.txt

echo "======================"
echo "======================"

echo " mafft: "
time mafft --thread 16 ${dataNmae} > ${alignName}

echo "======================"
echo "======================"

outName="./IqOut/"
echo " iqtree: "
./iqtree-3.0.1-Linux-intel/bin/iqtree3 -s ./${alignName} -m GTR+G -B 1000 -T AUTO --prefix ${outName}/result


echo "======================"
echo "======================"

outName="./IqOut/"

echo " TreeDist: "
Rscript compare_rf.R ${outName}/result.bionj ./out.tree 
rm -r ${outName}/*

