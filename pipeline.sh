#!/bin/bash
#SBATCH --job-name=TRIM                  # 作业名称
#SBATCH --cpus-per-task=20               # 申请的CPU数（需大于等于并行线程数 ）
#SBATCH --mem=64G                        # 申请内存大小
#SBATCH --time=15-00:00:00               # 最大运行时间
#SBATCH -D /home/hzclab/TRIM-test        # 工作目录（TRIM所在目录）

#=================================软件配置=====================================
#CONDA环境
CONDA_BASE=$(conda info --base)
#软件目录
FOLDX_DIR="/home/hzclab/software/foldx5.1"
ROSETTA_DIR="/home/hzclab/software/rosetta/source/bin"
BASE_DIR="/home/hzclab/TRIM-test"
#复合物结构路径
PDB="$BASE_DIR/pdb/2OOB.pdb"
# 受体蛋白链(若蛋白为多条链："ABCD")
rec_chains="B"
# 配体蛋白链（改造目标）
lig_chains="A"
#并行线程数
THREAD=20
#======================================================================================


OUTDIR="$BASE_DIR/out"
pdb_name=$(basename "$PDB" .pdb)

echo "============== Processing $pdb_name =============="
exec 1>>"$BASE_DIR/log/${pdb_name}.out" 2>>"$BASE_DIR/log/${pdb_name}.err"
out="$OUTDIR/${pdb_name}_out"
mkdir -p "$out"

mutout="$out/mutation"
energy_out="$out/energy"
result="$out/result"
mkdir -p "$mutout" "$energy_out" "$result"

# 蛋白互作分析模块
bash $BASE_DIR/script/interaction_analysis.sh "$PDB" "$out" "$rec_chains" "$lig_chains" "$pdb_name" "$BASE_DIR" "$CONDA_BASE" "$result"

# 饱和突变模拟模块
bash $BASE_DIR/script/Energy_calculate.sh "$PDB" "$out" "$rec_chains" "$lig_chains" "$pdb_name" "$mutout" "$energy_out" "$result" "$FOLDX_DIR" "$BASE_DIR" "$CONDA_BASE" "$THREAD"

# 突变体评估模块 
bash $BASE_DIR/script/mutation_evaluate.sh "$PDB" "$out" "$rec_chains" "$lig_chains" "$result"  "$ROSETTA_DIR" "$BASE_DIR" "$CONDA_BASE" "$THREAD"
echo "============== Processing $pdb_name Done=============="
