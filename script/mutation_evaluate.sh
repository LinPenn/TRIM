#!/bin/bash

pdb="$1"
out="$2"
rec_chains="$3"
lig_chains="$4"
result="$5"
ROSETTA_DIR="$6"
BASE_DIR="$7"
CONDA_BASE="$8"
THREAD="$9"

cd $BASE_DIR

source $CONDA_BASE/etc/profile.d/conda.sh
source $BASE_DIR/utils.sh

conda activate trim

pdb_name=$(basename "$pdb" .pdb)
rec=$(echo "$rec_chains" | sed -E "s/(.)/'\1', /g" | sed 's/, $//')
lig=$(echo "$lig_chains" | sed -E "s/(.)/'\1', /g" | sed 's/, $//')
chain="[[$rec], [$lig]]"

echo "========================================================================="
echo "=================Mutant Interaction Assessment Start====================="
echo -e "=========================================================================\n"

echo "[1] Extract Protein Chains Start"
#提取配体蛋白链
pymol -cq -r $BASE_DIR/tools/pymol_chains.py -- extract \
    -i $pdb \
    -c $lig_chains \
    -o "$out/${pdb_name}_${lig_chains}.pdb"
echo "  [+] 已提取配体蛋白链->${out}/${pdb_name}_${lig_chains}.pdb"

#提取受体蛋白链
pymol -cq -r $BASE_DIR/tools/pymol_chains.py -- extract \
    -i $pdb \
    -c $rec_chains \
    -o "$out/${pdb_name}_${rec_chains}.pdb"
echo "  [+] 已提取配体蛋白链->${out}/${pdb_name}_${rec_chains}.pdb"
echo -e "[1] Extract Protein Chains End\n"

host_pdb="$out/${pdb_name}_${lig_chains}.pdb"
virus_pdb="$out/${pdb_name}_${rec_chains}.pdb"
RES="$result/filtered_ddg_mutations.csv"
docking="$out/docking"

mkdir -p $out/resfiles
mkdir -p $docking/best/plip_result
mkdir -p $docking/best/analysis

echo "[2] Create resfiles for mutation Start:"
joblist="$out/joblist.tsv"
: > "$joblist"

tail -n +2 "$RES" | while IFS=',' read -r chain resi mut_aa _; do
    aa1=$(convert_to_one_letter "$mut_aa")
    if [[ -z "$aa1" ]]; then
        echo "  [WARNING] Unknown AA: $mut_aa"
        continue
    fi
    mut_name="${chain}${resi}${aa1}"
    workdir="$out/resfiles/$mut_name"
    mkdir -p "$workdir"
    create_resfile "$resi" "$chain" "$aa1" > "$workdir/resfile.txt"
    echo -e "$mut_name\t$workdir" >> "$joblist"
    echo "  [PREPARED] $mut_name"
done
echo -e "[2] Create resfiles for mutation End\n"

echo "[3] Mutate protein ${pdb_name} Start"
cat "$joblist" |
xargs -P "$THREAD" -n 2 bash -c '
mut_name="$1"
workdir="$2"

echo "  [INFO] FixBB $mut_name" >> "'"$BASE_DIR/log/${pdb_name}.out"'"

'"$ROSETTA_DIR"'/fixbb.linuxgccrelease \
    -s "'"$host_pdb"'" \
    -resfile "$workdir/resfile.txt" \
    -ex1 -ex2 -use_input_sc \
    -nstruct 1 \
    -mute all \
    -out:path:all "$workdir" \
    -out:suffix "_fixbb" \
    -overwrite \
    >> "'"$BASE_DIR/log/${pdb_name}_rosetta.out"'" \
    2>> "'"$BASE_DIR/log/${pdb_name}_rosetta.err"'" || exit 1

fixbb_pdb=$(ls "$workdir"/*_fixbb*.pdb 2>/dev/null | head -1)
[[ -f "$fixbb_pdb" ]] || exit 1

echo "  [INFO] Relax $mut_name" >> "'"$BASE_DIR/log/${pdb_name}.out"'"

'"$ROSETTA_DIR"'/relax.linuxgccrelease \
    -s "$fixbb_pdb" \
    -relax:fast \
    -relax:constrain_relax_to_start_coords \
    -use_input_sc \
    -nstruct 1 \
    -ex1 -ex2 \
    -score:weights ref2015 \
    -mute all \
    -out:path:all "$workdir" \
    -out:suffix "_relax" \
    -overwrite \
    >> "'"$BASE_DIR/log/${pdb_name}_rosetta.out"'" \
    2>> "'"$BASE_DIR/log/${pdb_name}_rosetta.err"'"

echo "  [DONE] $mut_name" >> "'"$BASE_DIR/log/${pdb_name}.out"'"
' _

shopt -s nullglob
echo -e "[3] Mutate protein ${pdb_name} Done\n"

echo "[4] Merge receptor-Ligand Chains Strat:"
# 拼接配体-受体蛋白链
for resdir in "$out/resfiles"/*/; do
    res_name=$(basename "$resdir")
    pymol -cq -r "$BASE_DIR/tools/pymol_chains.py" -- merge \
        -i $resdir/${pdb_name}_${lig_chains}_fixbb_0001_relax_0001.pdb $out/${pdb_name}_${rec_chains}.pdb \
        -o "$out/resfiles/${pdb_name}_${res_name}.pdb"
    echo "Succeed in merging ${pdb_name}_${res_name}"
done
echo -e "[4] Merge receptor-Ligand Chains Done\n"

echo "[5] Local_Docking Start:"
# 对拼接蛋白进行局部对接
find "$out/resfiles" -maxdepth 1 -name "*.pdb" -print0 |
xargs -0 -P "$THREAD" -I {} bash -c '
mut="$1"; partners="$2"; rosetta="$3"; docking="$4"; logdir="$5"

pdbname=$(basename "$mut" .pdb)
outlog="$logdir/'"$pdb_name"'_rosetta.out"
errlog="$logdir/'"$pdb_name"'_rosetta.err"
scriptoutlog=""$logdir/'"$pdb_name"'.out""
scripterrlog=""$logdir/'"$pdb_name"'.err""
{
  echo "  [$(date "+%F %T")] [START docking] $pdbname  file=$mut"
} >>"$scriptoutlog" 2>>"$scripterrlog"
  "$rosetta"/docking_protocol.linuxgccrelease \
    -s "$mut" \
    -partners "$partners" \
    -docking_local_refine \
    -use_input_sc \
    -docking:sc_min \
    -ex1 -ex2aro -spin \
    -no_optH false \
    -flip_HNQ true \
    -nstruct 1 \
    -score:weights ref2015 \
    -mute all \
    -out:file:silent "${pdbname}_dock.out" \
    -out:path:all "$docking" \
    >>"$outlog" 2>>"$errlog"
  rc=$?
{  
  echo "  [$(date "+%F %T")] [END docking] $pdbname  rc=$rc"
} >>"$scriptoutlog" 2>>"$scripterrlog"
  exit $rc
' _ {} "${rec_chains}_${lig_chains}" "$ROSETTA_DIR" "$docking" "$BASE_DIR/log"
echo -e "[5] Local_Docking Done!\n"

echo "[6] Extract Protein Structure Start:"
#从静默文件中提取蛋白结构
for silent in $docking/*_dock.out; do
    pdbname=$(basename "$silent" _dock.out)
    best_tag=$(
        awk '
        $1=="SCORE:" {
            if (!header_found) {
                for (i=1;i<=NF;i++)
                    if ($i=="score") col=i
                header_found=1
                next
            }
            print $col, $NF
        }
        ' "$silent" | sort -n | head -1 | awk '{print $2}'
    )
    # 提取该构象
    (cd $docking/best
    $ROSETTA_DIR/extract_pdbs.linuxgccrelease \
        -mute all \
        -in:file:silent "$silent" \
        -in:file:tags "$best_tag") \
        >>"$BASE_DIR/log/${pdb_name}_rosetta.out" 2>>"$BASE_DIR/log/${pdb_name}_rosetta.err"
    echo "  [+] ${best_tag} Done"
done
echo -e "[6] Extract Protein Structure End\n"

cp $pdb $docking/best

echo "[7] Assess Protein Complex Interface Start"
for pdb_best in $docking/best/*.pdb; do
    basename=$(basename "$pdb_best" .pdb)

    #对精修结构的互作界面进行评分
    $ROSETTA_DIR/InterfaceAnalyzer.linuxgccrelease \
        -s $pdb_best \
        -interface ${rec_chains}_${lig_chains} \
        -scorefxn ref2015 \
        -pack_input false \
        -pack_separated false \
        -mute all \
        -out:file:score_only $docking/best/score.sc \
        >>"$BASE_DIR/log/${pdb_name}_rosetta.out" 2>>"$BASE_DIR/log/${pdb_name}_rosetta.err"
    # 分析互作残基信息
    plip -f $pdb_best -o $docking/best/plip_result --chains "$chain" -qx --name $basename >>"$BASE_DIR/log/${pdb_name}_rosetta.out" 2>>"$BASE_DIR/log/${pdb_name}_rosetta.err"
    echo "[✓] $basename Done！"
done

# 提取互作残基信息
for xml in $docking/best/plip_result/*.xml; do
    python $BASE_DIR/tools/plip_extract.py -i $xml -o $docking/best/analysis
done
echo -e "[7] Assess Protein Complex Interface End\n"

rm $docking/best/plip_result/*.pdb

echo "Summary Result Start"
# 汇总残基信息以及界面评分
python $BASE_DIR/tools/summary.py \
        -p $docking/best/analysis \
        -s $docking/best/score.sc \
        -n $pdb_name \
        -m $out/PLIP_${pdb_name}_chain${lig_chains}_residues.txt \
        -o $result/interaction_summary.csv

# 突变体评估结果可视化
python $BASE_DIR/tools/evaluate_plot.py \
        -i $result/interaction_summary.csv \
        -o $result \
        -n $pdb_name
echo -e "Summary Result End\n"
conda deactivate

echo "========================================================================="
echo "==================Mutant Interaction Assessment End======================"
echo -e "=========================================================================\n"
