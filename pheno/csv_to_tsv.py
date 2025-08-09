r"""
将 CSV 文件转换为 TSV 文件的简单脚本。

功能：
- 使用 Python 内置 csv 库按流式读取 CSV，写出为制表符分隔的 TSV；
- 默认输入/输出路径已设置，可用命令行参数覆盖；
- 自动处理 UTF-8 BOM（默认使用 utf-8-sig 读取）。

示例：
    python csv_to_tsv.py \
        --input C:\Users\Cloud\Desktop\showANDdnngp\pheno\maize.hybrid.train_phe.csv \
        --output C:\Users\Cloud\Desktop\showANDdnngp\pheno\maize.hybrid.train_phe.tsv
"""

import argparse
import os
import sys
import csv


def convert_csv_to_tsv(input_csv: str, output_tsv: str, encoding: str = "utf-8-sig") -> None:
    """将 CSV 转换为 TSV（流式处理，适合大文件）。"""
    with open(input_csv, "r", encoding=encoding, newline="") as fin, \
         open(output_tsv, "w", encoding="utf-8", newline="") as fout:
        reader = csv.reader(fin, delimiter=",", quotechar='"')
        writer = csv.writer(fout, delimiter="\t", quotechar='"', lineterminator="\n")
        for row in reader:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="将 CSV 文件转换为 TSV 文件")
    parser.add_argument(
        "--input",
        dest="input_csv",
        default=r"C:\\Users\\Cloud\\Desktop\\showANDdnngp\\pheno\\maize.hybrid.train_phe.csv",
        help="输入 CSV 文件路径（默认指向项目内示例文件）",
    )
    parser.add_argument(
        "--output",
        dest="output_tsv",
        default=r"C:\\Users\\Cloud\\Desktop\\showANDdnngp\\pheno\\maize.hybrid.train_phe.tsv",
        help="输出 TSV 文件路径（默认与输入同目录、同名不同后缀）",
    )
    parser.add_argument(
        "--encoding",
        dest="encoding",
        default="utf-8-sig",
        help="输入文件编码（默认 utf-8-sig，可处理 UTF-8 BOM）",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"[Error] 输入文件不存在: {args.input_csv}")
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)

    print(f"[Info] 读取 CSV: {args.input_csv}")
    print(f"[Info] 写入 TSV: {args.output_tsv}")
    try:
        convert_csv_to_tsv(args.input_csv, args.output_tsv, encoding=args.encoding)
    except Exception as e:
        print(f"[Error] 转换失败: {e}")
        sys.exit(1)
    print("[Done] 转换完成")


if __name__ == "__main__":
    main()


