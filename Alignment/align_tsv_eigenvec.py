r"""
将 TSV 与 PLINK eigenvec 文件按个体进行名称与顺序对齐的脚本。

规则与功能：
- 输入一个 TSV（个体位于行，含个体 ID 列）与一个 eigenvec 文件；
- 对齐后输出：
  1) 一个仅保留在两者交集中的个体、且顺序与 eigenvec 一致的 TSV；
  2) 一个仅保留交集个体、并将 ID 统一为对齐后名称的 eigenvec；
- 当 TSV 的个体名包含两部分（例如 "A/B"），会用分隔符（默认'/'）拆分为两个候选，与 eigenvec 的 ID 对比：匹配到哪一部分就保留哪一部分；
- 在新文件中（TSV 与 eigenvec）均使用匹配保留成功的那一部分作为个体 ID，其它部分会被去掉；
- 若 eigenvec 有 FID+IID 两列，则输出时会将两列统一设置为对齐后的个体 ID；若只有 IID 一列，则替换该列。

新增：
- 自动移除“空/无效 ID”（空串、仅空白、NA、.）的个体，并在两个输出文件中同时剔除，确保一一对应；
- 若 TSV 中同一对齐后 ID 出现多行，仅保留第一行，其余忽略并提示数量；
- 保持输出顺序以 eigenvec 为准。

注意与假设：
- 假定 TSV 为“行=个体”的宽表结构，`--id-col` 指定个体 ID 所在列（基于 0 的索引，默认 0 列）。
- 自动检测 eigenvec 是否为 [IID, PCs...] 或 [FID, IID, PCs...] 结构；也可用 `--eigenvec-id` 强制指定使用 `iid` 或 `fid`。
 - 可选过滤非有限数（NaN/Inf）：通过参数在 TSV 与 eigenvec 两侧同时剔除含非有限数的样本，避免下游 pearsonr 报错。

使用示例（PowerShell，可直接复制运行并替换为你的真实路径）：
  python C:\Users\Cloud\Desktop\showANDdnngp\Alignment\align_tsv_eigenvec.py ^
    --tsv C:\Users\Cloud\Desktop\showANDdnngp\pheno\maize.hybrid.train_phe.tsv ^
    --eigenvec C:\Users\Cloud\Desktop\showANDdnngp\SNP\pca10.eigenvec ^
    --out-tsv C:\Users\Cloud\Desktop\showANDdnngp\Alignment\maize.hybrid.train_phe.aligned.tsv ^
    --out-eigenvec C:\Users\Cloud\Desktop\showANDdnngp\Alignment\pca10.aligned.eigenvec ^
    --id-col 0 ^
    --split-char / ^
    --eigenvec-id auto

参数说明：
- --tsv: 需要对齐的 TSV 文件（行=个体，包含个体 ID 列）。
- --eigenvec: plink2 PCA 输出的 .eigenvec 文件。
- --out-tsv: 输出的对齐后 TSV 文件路径（建议置于 Alignment 目录）。
- --out-eigenvec: 输出的对齐后 eigenvec 文件路径（建议置于 Alignment 目录）。
- --id-col: TSV 中“个体 ID”所在列的索引（0 表示第 1 列，默认 0）。
- --split-char: 当 TSV 个体名形如 A/B 时，用该分隔符拆分优先匹配 A，其次匹配 B；匹配到谁就保留谁（默认 '/'）。
- --eigenvec-id: 自动判断或强制指定 eigenvec 的 ID 列（auto/iid/fid）。

运行前小贴士：
- 若没有 TSV，可先用 pheno/csv_to_tsv.py 将 CSV 转 TSV。
- 若没有 eigenvec，可先用 plink2 进行 PCA 生成（也可用 SNP/build_plink2_cmd.py 构造命令）。

运行结果：
- 仅保留 TSV 与 eigenvec 的交集个体，顺序以 eigenvec 为准；
- 若 TSV 个体名为 A/B，会替换为匹配成功的那一部分（A 或 B），并在两个新文件中统一使用；
- 输出文件为 --out-tsv 与 --out-eigenvec 指定的路径。
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Tuple


def is_empty_id(token: Optional[str]) -> bool:
    s = (token or "").strip()
    return s == "" or s.upper() == "NA" or s == "."


def is_finite_number(token: str) -> bool:
    try:
        import math
        x = float(str(token).strip())
        return math.isfinite(x)
    except Exception:
        return False


def is_empty_id(token: Optional[str]) -> bool:
    s = (token or "").strip()
    return s == "" or s.upper() == "NA" or s == "."


def try_parse_float(token: str) -> bool:
    try:
        float(token)
        return True
    except Exception:
        return False


def detect_eigenvec_id_mode(first_tokens: List[str]) -> Tuple[str, int, Optional[int]]:
    """基于首行 tokens 自动判断 eigenvec 的 ID 列模式。

    返回 (mode, id_idx, fid_idx)。mode 取值：
      - 'iid_only': 仅有 IID 一列（通常 token1=IID, token2=PC1 数值）
      - 'fid_iid': 同时有 FID 与 IID 两列（通常 token3=PC1 数值）
      - 'unknown': 无法可靠判断，退化为将 token0 当作 ID
    id_idx 为所用 ID 列索引；fid_idx 若存在则给出，否则为 None。
    """
    if not first_tokens:
        return ("unknown", 0, None)

    if len(first_tokens) >= 2 and try_parse_float(first_tokens[1]):
        # 形如: IID PC1 PC2 ...
        return ("iid_only", 0, None)

    if len(first_tokens) >= 3 and try_parse_float(first_tokens[2]):
        # 形如: FID IID PC1 PC2 ...
        return ("fid_iid", 1, 0)

    return ("unknown", 0, None)


def read_eigenvec(
    eigenvec_path: str,
    eigenvec_id_mode: str = "auto",
) -> Tuple[Optional[List[str]], List[List[str]], int, Optional[int]]:
    """读取 eigenvec 文件，检测并保留表头。

    返回 (header_tokens 或 None, 数据行 tokens 列表, id_idx, fid_idx)。
    """
    header_tokens: Optional[List[str]] = None
    data_lines: List[List[str]] = []

    def is_header(tokens: List[str]) -> bool:
        upper = [t.upper() for t in tokens]
        if any(t in {"FID", "IID"} for t in upper):
            return True
        if any(t.startswith("PC") for t in upper):
            return True
        # 若大部分（>=2）后续列非数值，也可能是表头
        non_numeric = sum(1 for t in tokens[1:] if not try_parse_float(t))
        return non_numeric >= 2

    with open(eigenvec_path, "r", encoding="utf-8", newline="") as fin:
        for raw in fin:
            s = raw.strip()
            if not s:
                continue
            tokens = s.split()
            if header_tokens is None and is_header(tokens):
                header_tokens = tokens
                continue
            data_lines.append(tokens)

    if not data_lines:
        raise ValueError("eigenvec 文件为空或仅包含表头")

    first_data = data_lines[0]
    if eigenvec_id_mode == "auto":
        mode, id_idx, fid_idx = detect_eigenvec_id_mode(first_data)
        return (header_tokens, data_lines, id_idx, fid_idx)
    elif eigenvec_id_mode == "iid":
        # 假定第二列为 IID（若不足两列则回退第一列）
        id_idx = 1 if len(first_data) >= 2 else 0
        fid_idx = 0 if len(first_data) >= 2 else None
        return (header_tokens, data_lines, id_idx, fid_idx)
    elif eigenvec_id_mode == "fid":
        return (header_tokens, data_lines, 0, None)
    else:
        raise ValueError("--eigenvec-id 仅支持 auto|iid|fid")


def read_tsv(tsv_path: str, id_col: int) -> Tuple[List[str], List[List[str]]]:
    """读取 TSV，返回 (header, rows)。"""
    with open(tsv_path, "r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError("TSV 文件为空")
    header = rows[0]
    data_rows = rows[1:]
    if id_col < 0 or id_col >= len(header):
        raise ValueError(f"--id-col 超出列范围：{id_col}")
    return header, data_rows


def build_id_map_from_tsv(
    tsv_rows: List[List[str]],
    id_col: int,
    eigen_ids: set,
    split_char: str,
) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """从 TSV 行构建映射：已对齐 ID -> 行 (按原顺序收集)，以及 原 ID -> 对齐后 ID。

    当 TSV 的 ID 含分隔符时，拆分两部分进行与 eigen IDs 的对比；
    - 若 part1 命中，则使用 part1；若 part2 命中，则使用 part2；两者都命中时优先 part1；
    - 若均未命中但完整 ID 命中，也可保留完整 ID；否则丢弃该行（无法对齐）。
    """
    aligned_id_to_rows: Dict[str, List[str]] = {}
    original_to_aligned: Dict[str, str] = {}

    for row in tsv_rows:
        if id_col >= len(row):
            continue
        raw_id = (row[id_col] or "").strip()
        if is_empty_id(raw_id):
            continue

        chosen: Optional[str] = None
        if split_char and split_char in raw_id:
            parts = [p.strip() for p in raw_id.split(split_char) if p.strip()]
            # 优先匹配第一个部分，再匹配第二个部分
            for p in parts:
                if p in eigen_ids:
                    chosen = p
                    break
        # 回退：完整 ID 直接命中也允许
        if chosen is None and raw_id in eigen_ids:
            chosen = raw_id

        if chosen is None:
            continue  # 无法对齐的个体，丢弃

        original_to_aligned[raw_id] = chosen
        # 将行的 ID 替换为对齐后的 ID
        new_row = list(row)
        new_row[id_col] = chosen
        aligned_id_to_rows.setdefault(chosen, []).append(new_row)

    return aligned_id_to_rows, original_to_aligned


def write_tsv(path: str, header: List[str], rows: List[List[str]]):
    with open(path, "w", encoding="utf-8", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="将 TSV 与 eigenvec 个体对齐（名称与顺序一致）")
    parser.add_argument("--tsv", required=True, help="输入 TSV 文件路径（行=个体，包含个体 ID 列）")
    parser.add_argument("--eigenvec", required=True, help="输入 eigenvec 文件路径（PLINK PCA 输出）")
    parser.add_argument("--out-tsv", required=True, help="输出对齐后的 TSV 文件路径")
    parser.add_argument("--out-eigenvec", required=True, help="输出对齐后的 eigenvec 文件路径")
    parser.add_argument("--id-col", type=int, default=0, help="TSV 中个体 ID 所在列（从 0 开始，默认 0）")
    parser.add_argument("--split-char", default="/", help="TSV 个体名分隔符（默认 '/'）")
    parser.add_argument("--eigenvec-id", choices=["auto", "iid", "fid"], default="auto", help="eigenvec 使用的 ID 列：auto 自动检测，或强制使用 iid/fid")

    args = parser.parse_args()

    if not os.path.isfile(args.tsv):
        print(f"[Error] TSV 文件不存在: {args.tsv}")
        sys.exit(1)
    if not os.path.isfile(args.eigenvec):
        print(f"[Error] eigenvec 文件不存在: {args.eigenvec}")
        sys.exit(1)

    # 读取 eigenvec 并确定使用的 ID 列
    header_tokens, eigen_lines, id_idx, fid_idx = read_eigenvec(
        args.eigenvec, eigenvec_id_mode=args.eigenvec_id
    )
    # 先移除 eigenvec 中空/无效 ID 的行
    before_eigen = len(eigen_lines)
    eigen_lines = [tok for tok in eigen_lines if (len(tok) > id_idx and not is_empty_id(tok[id_idx]))]
    removed_eigen_empty = before_eigen - len(eigen_lines)

    # 可选：过滤 eigenvec 中含非有限数（NaN/Inf）的样本行
    removed_nonfinite_eigen = 0
    if '--filter-finite-eigenvec' in sys.argv:
        cleaned = []
        for tok in eigen_lines:
            ok = True
            for v in tok[id_idx + 1:]:
                if not is_finite_number(v):
                    ok = False
                    break
            if ok:
                cleaned.append(tok)
            else:
                removed_nonfinite_eigen += 1
        eigen_lines = cleaned

    eigen_ids_in_order: List[str] = [tok[id_idx] for tok in eigen_lines]
    eigen_id_set = set(eigen_ids_in_order)

    # 读取 TSV 并根据 eigen IDs 构建映射
    header, tsv_rows = read_tsv(args.tsv, args.id_col)
    # 可选：基于指定列过滤 TSV 的非有限数样本（以命令行原始参数触发，避免破坏现有接口）
    removed_nonfinite_tsv = 0
    tsv_cols_arg_prefix = '--filter-finite-tsv-cols='
    raw_arg = next((a for a in sys.argv if a.startswith(tsv_cols_arg_prefix)), None)
    cols_idx: List[int] = []
    if raw_arg is not None:
        try:
            cols_idx = [int(x) for x in raw_arg[len(tsv_cols_arg_prefix):].split(',') if str(x).strip() != '']
        except Exception:
            print("[Error] --filter-finite-tsv-cols 解析失败，应为以逗号分隔的列索引，如 0,2,5")
            sys.exit(1)
    if cols_idx:
        filtered_rows: List[List[str]] = []
        for row in tsv_rows:
            ok = True
            for ci in cols_idx:
                if ci < 0 or ci >= len(row):
                    continue
                if not is_finite_number(row[ci]):
                    ok = False
                    break
            if ok:
                filtered_rows.append(row)
            else:
                removed_nonfinite_tsv += 1
        tsv_rows = filtered_rows
    aligned_id_to_rows, original_to_aligned = build_id_map_from_tsv(
        tsv_rows=tsv_rows,
        id_col=args.id_col,
        eigen_ids=eigen_id_set,
        split_char=args.split_char,
    )

    # 按 eigenvec 顺序生成对齐后的 TSV 与 eigenvec
    aligned_tsv_rows: List[List[str]] = []
    aligned_eigen_lines: List[str] = []

    kept = 0
    skipped_due_to_dup = 0
    for tokens in eigen_lines:
        eigen_id = tokens[id_idx]
        if eigen_id not in aligned_id_to_rows:
            continue  # 个体不在 TSV 中（或未能匹配），跳过

        # 选择该 eigen_id 对应的第一条 TSV 行（如有重复，仅保留第一条）
        rows_for_id = aligned_id_to_rows[eigen_id]
        row = rows_for_id[0]
        if len(rows_for_id) > 1:
            skipped_due_to_dup += len(rows_for_id) - 1
        aligned_tsv_rows.append(row)

        # 构造新的 eigenvec 行：将 ID 统一为 row[id_col]
        chosen_id = row[args.id_col]
        new_tokens = list(tokens)
        if fid_idx is not None and 0 <= fid_idx < len(new_tokens):
            new_tokens[fid_idx] = chosen_id
        if id_idx is not None and 0 <= id_idx < len(new_tokens):
            new_tokens[id_idx] = chosen_id
        aligned_eigen_lines.append(" ".join(new_tokens))
        kept += 1

    if kept == 0:
        print("[Error] 没有找到可对齐的个体（交集为空）。请检查 ID 列与分隔符设置。")
        sys.exit(2)

    # 写出文件
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_eigenvec) or ".", exist_ok=True)
    write_tsv(args.out_tsv, header, aligned_tsv_rows)
    with open(args.out_eigenvec, "w", encoding="utf-8", newline="\n") as fout:
        if header_tokens is not None:
            fout.write(" ".join(header_tokens) + "\n")
        for line in aligned_eigen_lines:
            fout.write(line + "\n")

    print(f"[Done] 对齐完成：保留个体 {kept} 个")
    if removed_eigen_empty > 0:
        print(f"[Info] 已从 eigenvec 中移除空/无效 ID 个体：{removed_eigen_empty}")
    if skipped_due_to_dup > 0:
        print(f"[Info] TSV 中存在重复映射，除首个外已忽略重复行：{skipped_due_to_dup}")
    if removed_nonfinite_eigen > 0:
        print(f"[Info] 已从 eigenvec 中移除含非有限数的样本：{removed_nonfinite_eigen}")
    if removed_nonfinite_tsv > 0:
        print(f"[Info] 已从 TSV 中移除含非有限数的样本：{removed_nonfinite_tsv}")
    print(f"[Out] TSV: {args.out_tsv}")
    print(f"[Out] EIGEN: {args.out_eigenvec}")


if __name__ == "__main__":
    main()


