#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import csv
import os
import argparse

def write_csv(filename, rows, header):
    """写入 CSV 文件"""
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def convert_three_to_one(aa_three, aa_dict):

    # 单氨基酸
    if isinstance(aa_three, str):
    
        if aa_three in aa_dict:
            return aa_dict[aa_three]
        else:
            raise ValueError(f"未知的氨基酸三字母代码: {aa_three}")

    # 多氨基酸
    elif isinstance(aa_three, list):
        result = []
        for code in aa_three:
            if code in aa_dict:
                result.append(aa_dict[code])
            else:
                raise ValueError(f"未知的氨基酸三字母代码: {code}")
        # 连接成一个字符串返回
        return ''.join(result)

    else:
        raise TypeError("输入数据必须是字符串或列表")

def parse_plip_xml(xml_file, outdir):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    res_dict = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }

    # 存储不同相互作用类型的残基对与相关参数
    interactions = {
        "hydrophobic_interactions": {
            "rows": [],
            "header": ["rec", "lig", "distance"]
        },
        "hydrogen_bonds": {
            "rows": [],
            "header": ["rec", "lig", "distance", "angle"]
        },
        "salt_bridges": {
            "rows": [],
            "header": ["rec", "lig", "distance"]
        },
        "pi_stacks": {
            "rows": [],
            "header": ["rec", "lig", "distance", "angle", "offset"]
        },
        "pi_cation_interactions": {
            "rows": [],
            "header": ["rec", "lig", "distance"]
        },
        "halogen_bonds": {
            "rows": [],
            "header": ["rec", "lig", "distance", "angle"]
        },
        "metal_complexes": {
            "rows": [],
            "header": ["rec", "lig", "distance"]
        }
    }

    # 遍历不同的 interaction 类型
    for itype, data in interactions.items():
        for inter in root.findall(f".//{itype}/*"):
            # 受体残基
            rec_resnr = inter.findtext("resnr")
            rec_restype = convert_three_to_one(inter.findtext("restype"), res_dict)
            rec_reschain = inter.findtext("reschain")
            rec = f"{rec_restype}{rec_reschain}{rec_resnr}" if rec_resnr and rec_restype and rec_reschain else "NA"

            # 配体残基
            lig_resnr = inter.findtext("resnr_lig")
            lig_restype = convert_three_to_one(inter.findtext("restype_lig"), res_dict)
            lig_reschain = inter.findtext("reschain_lig")
            lig = f"{lig_restype}{lig_reschain}{lig_resnr}" if lig_resnr and lig_restype and lig_reschain else "NA"

            # 公共参数
            hbond_dict = inter.findtext("dist_h-a")
            hbond_angle = inter.findtext("don_angle")
            distance = inter.findtext("distance") or inter.findtext("dist") or ""
            angle = inter.findtext("angle") or ""
            offset = inter.findtext("offset") or ""

            # 按类型保存
            if itype == "hydrogen_bonds":
                data["rows"].append([rec, lig, hbond_dict, hbond_angle])
            elif itype == "salt_bridges":
                data["rows"].append([rec, lig, distance])
            elif itype == "hydrophobic_interactions":
                data["rows"].append([rec, lig, distance])
            elif itype == "pi_stacks":
                data["rows"].append([rec, lig, distance, angle, offset])
            elif itype == "pi_cation_interactions":
                data["rows"].append([rec, lig, distance])
            elif itype == "halogen_bonds":
                data["rows"].append([rec, lig, distance, angle])
            elif itype == "metal_complexes":
                data["rows"].append([rec, lig, distance])

    # 写入文件
    prefix = os.path.splitext(os.path.basename(xml_file))[0]
    os.makedirs(outdir, exist_ok=True)
    for itype, data in interactions.items():
        if data["rows"]: 
            fname = os.path.join(outdir, f"PLIP_{prefix}_{itype}.csv")
            write_csv(fname, data["rows"], data["header"])
            print(f"[+] 写入 {fname}，共 {len(data['rows'])} 条")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="解析 PLIP XML 并提取残基对及互作信息")
    parser.add_argument("-i", "--input", required=True, help="PLIP 生成的 XML 文件")
    parser.add_argument("-o", "--outdir", default="plip_csv", help="输出目录")
    args = parser.parse_args()

    parse_plip_xml(args.input, args.outdir)
