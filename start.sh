#/bin/bash/
set -e

trap 'echo "Error occurred at line $LINENO"' ERR

cmake -B build > cmake_log.txt 2>&1
make -C build > build_log.txt 2>&1

# ./data/test.fas

# ./data/Variola.fasta
# data/sars2_similarity_70_1.fas
# data/23s_rRNA.fasta
# mt1x.fasta
# steps-1024_0

dataNmae="./data/steps-1024_0.fasta"

alignName="align.fasta"

echo "======================"
echo "======================"

echo " GLS: "

/usr/bin/time -v ./build/GLS -i ${dataNmae} -k 6 -t 1 > log.txt

# echo "======================"
# echo "======================"

# echo " mafft: "
# time mafft --thread 16 ${dataNmae} > ${alignName}

# echo "======================"
# echo "======================"

# outName="./IqOut/"
# echo " iqtree: "
# ./iqtree-3.0.1-Linux-intel/bin/iqtree3 -s ./${alignName} -m GTR+G -B 1000 -T AUTO --prefix ${outName}/result


# echo "======================"
# echo "======================"

# outName="./IqOut/"

# echo " TreeDist: "
# Rscript compare_rf.R ${outName}/result.bionj ./out.tree 
# rm -r ${outName}/*

