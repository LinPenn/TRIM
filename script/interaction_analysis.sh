#!/bin/bash

pdb="$1"
out="$2"
rec_chains="$3"
lig_chains="$4"
base="$5"
BASE_DIR="$6"
CONDA_BASE="$7"
result="$8"

source $CONDA_BASE/etc/profile.d/conda.sh
source $BASE_DIR/utils.sh
conda activate trim

echo "==============================================================================="
echo "=================Interaction Interface Identification Start===================="
echo -e "===============================================================================\n"

#将链信息转化为PLIP的输入格式
rec=$(echo "$rec_chains" | sed -E "s/(.)/'\1', /g" | sed 's/, $//')
lig=$(echo "$lig_chains" | sed -E "s/(.)/'\1', /g" | sed 's/, $//')
chain="[[$rec], [$lig]]"
echo "受体链：${rec}"
echo "配体链：${lig}"

#PLIP分析
plip -f "$pdb" -o "$out" --chains "$chain" -qxy --name $base

# 提取PLIP挖掘得到的互作残基
for xml in $out/*.xml; do
    python $BASE_DIR/tools/plip_extract.py -i "$xml" -o "$out"
done

mv $out/*.pse result
rm $out/*.pdb

#蛋白互作残基可视化
python $BASE_DIR/tools/interaction_plot.py -p ${out}/PLIP_${base} -o $out/result/Chord_Interaction.pdf -e $lig_chains

conda deactivate

echo "==============================================================================="
echo "==================Interaction Interface Identification End====================="
echo -e "===============================================================================\n"