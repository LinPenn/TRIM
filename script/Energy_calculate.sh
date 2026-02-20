#!/bin/bash

pdb="$1"
out="$2"
rec_chains="$3"
lig_chains="$4"
base="$5"
mutout="$6"
energy_out="$7"
result="$8"
foldx="$9"
BASE_DIR="${10}"
CONDA_BASE="${11}"
THREAD="${12}"
cd $BASE_DIR

source $CONDA_BASE/etc/profile.d/conda.sh
source $BASE_DIR/utils.sh

echo "========================================================================="
echo "=================Saturation Mutagenesis Simulation Start================="
echo -e "=========================================================================\n"

# 提取突变位点信息
PLIP_MUT=$(awk -F',' 'FNR==1 {col=0; for (i=1;i<=NF;i++) if ($i=="lig") col=i; next} col>0 {if (match($col, /[0-9]+/, a)) print a[0], $col}' $out/PLIP*.csv | sort -n -k1,1 | awk '{print $2 "a"}' | uniq | paste -sd, -)
echo "所有突变位点：${PLIP_MUT}"

echo "[1] Repair PDB Strat"
#修复蛋白结构
$foldx/foldx \
    --command=RepairPDB \
    --pdb-dir=$(dirname "$pdb") \
    --pdb=$(basename "$pdb") \
    --output-dir=$mutout \
    --screen=false \
    >>"$BASE_DIR/log/${base}_folx.out" 2>>"$BASE_DIR/log/${base}_foldx.err"
echo -e "[1] Repair PDB End\n"

cd $mutout

echo "[2] Build Mutants Library Start"
#构建突变体库
$foldx/foldx \
        --command=PositionScan \
        --pdb=${base}_Repair.pdb \
        --positions=$PLIP_MUT \
        --output-dir=$mutout \
        --screen=false \
        >>"$BASE_DIR/log/${base}_folx.out" 2>>"$BASE_DIR/log/${base}_foldx.err"
echo -e "[2] Build Mutants Library End\n"

echo "[3] Calculate Mutants Energy Start"
#计算突变体能量变化
COMPLEX="${rec_chains},${lig_chains}"

find . -maxdepth 1 -name "*.pdb" -print0 |
xargs -0 -P $THREAD -I {} bash -c '
mutpdb="$1"
complex="$2"
foldx="$3"
energy_out="$4"

name=$(basename "$mutpdb" .pdb)
echo "--> Calculating $name"

"$foldx" --command=AnalyseComplex \
  --pdb="$mutpdb" \
  --analyseComplexChains="$complex" \
  --output-dir="$energy_out" \
' _ {} "$COMPLEX" "$foldx/foldx" "$energy_out"
>>"$BASE_DIR/log/${base}_folx.out" 2>>"$BASE_DIR/log/${base}_foldx.err"
echo -e "[3] Calculate Mutants Energy End\n"

cd $BASE_DIR

# 将单字母格式转化为三字母格式
CONVERTED=$(convert_to_three_letter "$PLIP_MUT")

# 移除原氨基酸的能量计算结果
for res in $(echo "$CONVERTED" | tr ',' ' '); do
    rm -f $energy_out/*$res*
done

echo "[4] Calculate DDG of Mutants Start"
conda activate trim
#计算突变体的能量变化
python $BASE_DIR/tools/calculate_ddg_by_position.py \
        --wt $energy_out/Summary_${base}_Repair_AC.fxout \
        --indir $energy_out \
        --outdir $result \
        --name $base
echo -e "[4] Calculate DDG of Mutants End\n"
echo "[5] Screen mutants Start"
#依据阈值筛选突变体
python $BASE_DIR/tools/filter_high_ddg_mutations.py \
        --dir $result \
        --chain "$lig_chains" 
#筛选结果可视化
python $BASE_DIR/tools/bubble_heatmap.py \
        -d $result \
        -o $result/Bubble_Heatmap.pdf 
conda deactivate
echo -e "[5] Screen mutants End"

echo "========================================================================="
echo "=================Saturation Mutagenesis Simulation End==================="
echo -e "=========================================================================\n"


