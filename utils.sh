#!/usr/bin/env bash

# --- 定义1->3映射表 --- 
declare -A AA_MAP=( 
    [A]=ALA 
    [R]=ARG 
    [N]=ASN 
    [D]=ASP 
    [C]=CYS 
    [E]=GLU 
    [Q]=GLN 
    [G]=GLY 
    [H]=HIS 
    [I]=ILE 
    [L]=LEU 
    [K]=LYS 
    [M]=MET 
    [F]=PHE 
    [P]=PRO 
    [S]=SER 
    [T]=THR 
    [W]=TRP 
    [Y]=TYR 
    [V]=VAL ) 

# --- 定义3->1映射表 --- 
declare -A AA_MAP_REV
for k in "${!AA_MAP[@]}"; do
    AA_MAP_REV["${AA_MAP[$k]}"]=$k
done

# 单字母转化为三字母
convert_to_three_letter() { 
    local input=$1 
    local result="" 
        
    IFS=',' read -ra arr <<< "$input" 
    for item in "${arr[@]}"; do 
        local one=${item:0:1} 
        local num=$(echo "$item" | grep -o '[0-9]\+') 
        local three=${AA_MAP[$one]} 
        if [[ -n "$three" && -n "$num" ]]; then 
            result+="${three}${num}," 
        fi 
    done 
    echo "${result%,}" 
    }

# 三字母转化为单字母
convert_to_one_letter() {
    local aa3=$1
    aa3=$(echo "$aa3" | tr 'a-z' 'A-Z') 
    echo "${AA_MAP_REV[$aa3]}"
}

# 创建突变信息文件
create_resfile() {
    local resi=$1
    local chain=$2
    local aa=$3

    cat <<EOF
NATRO
start
${resi} ${chain} PIKAA ${aa}
EOF
}