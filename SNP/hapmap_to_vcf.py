r"""
HapMap(hmp.txt) -> VCF 转换脚本（中文使用说明）

功能：
- 解析 HapMap 格式（hmp.txt）的标记与样本基因型，输出符合 VCFv4.2 的文件；
- 识别等位基因字段（如 A/C、AC、A|C），支持 IUPAC 简并码（R,Y,W,S,K,M 等）；
- 缺失/不确定基因型输出为 ./.；单态位点 ALT 列输出为 .；格式列为 GT。

示例（Windows）：
  1) 在仓库根目录运行（使用默认示例文件，默认输出到 `SNP/origin_maize.geno.selected.vcf`）：
     PowerShell:
       python .\\SNP\\hapmap_to_vcf.py

  2) 指定输入/输出路径（单行）：
     PowerShell:
       python .\\SNP\\hapmap_to_vcf.py "C:\\path\\to\\input.hmp.txt" "C:\\path\\to\\output.vcf"
     CMD:
       python .\SNP\hapmap_to_vcf.py "C:\path\to\input.hmp.txt" "C:\path\to\output.vcf"

  3) 指定输入/输出路径（PowerShell 多行，使用反引号续行）：
       python .\\SNP\\hapmap_to_vcf.py `
         "C:\\path\\to\\input.hmp.txt" `
         "C:\\path\\to\\output.vcf"

说明：
- 输入必须是标准 HapMap 列表头（前 11 列为元数据，样本列从第 12 列开始）。
- 若 Alleles 列无法提供 REF/ALT，将从首个非缺失基因型中推断 REF，并收集 ALT；若仍无法判断，REF 用 N。
- 输出文件首部包含 VCF 头，FORMAT 仅含 GT（基因型）。
"""

import os
import sys
from datetime import datetime, timezone


def parse_alleles_field(alleles_raw: str) -> list:
    """Parse HapMap alleles field into an ordered list of unique bases.

    Accepts formats like "A/C", "AC", "A C", or "A|C". Returns only A,C,G,T.
    The first base is treated as REF, the rest as ALT(s).
    """
    if not alleles_raw:
        return []
    tokens = []
    for ch in alleles_raw.replace('/', ' ').replace('|', ' ').replace(',', ' ').split():
        for c in ch:
            tokens.append(c.upper())

    ordered = []
    for b in tokens:
        if b in {"A", "C", "G", "T"} and b not in ordered:
            ordered.append(b)
    return ordered


IUPAC_MAP = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": set(),
    "-": set(),
}


def normalize_genotype_raw(gt_raw: str):
    """Normalize a raw HapMap genotype cell to a tuple of alleles (a1, a2).

    Returns:
      (allele1, allele2): each in {A,C,G,T}, or None if missing/ambiguous.
    """
    if gt_raw is None:
        return None
    s = str(gt_raw).strip().upper()
    if s == "" or s == "NA":
        return None

    # Single-character IUPAC handling
    if len(s) == 1 and s.isalpha():
        if s in IUPAC_MAP:
            bases = IUPAC_MAP[s]
            if len(bases) == 1:
                b = next(iter(bases))
                return (b, b)
            elif len(bases) == 2:
                a1, a2 = sorted(list(bases))  # deterministic order
                return (a1, a2)
            else:
                return None  # 0 or >2 -> treat as missing/ambiguous
        # Unknown letter -> missing
        return None

    # Multi-character formats like "A/G", "AG", "N/N", "A|C"
    cleaned = []
    for ch in s:
        if ch in {"A", "C", "G", "T", "N"}:
            cleaned.append(ch)
        # ignore separators like '/', '|', '\\', '-', ':'
    if not cleaned:
        return None
    # Consider only first two base-like chars
    if len(cleaned) < 2:
        # If only one base and not N, duplicate it (homozygous)
        if cleaned[0] in {"A", "C", "G", "T"}:
            return (cleaned[0], cleaned[0])
        return None

    a1, a2 = cleaned[0], cleaned[1]
    if a1 == "N" or a2 == "N":
        return None
    if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"}:
        return None
    return (a1, a2)


def write_vcf_header(out_fh, sample_ids):
    out_fh.write("##fileformat=VCFv4.2\n")
    out_fh.write(f"##source=hapmap_to_vcf.py ({datetime.now(timezone.utc).isoformat()})\n")
    out_fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    header_cols.extend(sample_ids)
    out_fh.write("\t".join(header_cols) + "\n")


def convert_hapmap_to_vcf(hapmap_path: str, vcf_path: str):
    with open(hapmap_path, "r", encoding="utf-8") as hin, open(vcf_path, "w", encoding="utf-8") as hout:
        header = hin.readline()
        if not header:
            raise ValueError("Empty HapMap file.")
        header_fields = header.rstrip("\n\r").split("\t")
        if len(header_fields) < 12 or header_fields[0].lower().startswith("rs#") is False:
            # Some HapMap headers may not start with 'rs#' exactly; do a looser check
            # but still expect at least 12 columns
            if len(header_fields) < 12:
                raise ValueError("Unexpected HapMap header: fewer than 12 columns.")

        # Standard HapMap: first 11 metadata columns, samples start at index 11
        sample_ids = header_fields[11:]
        write_vcf_header(hout, sample_ids)

        line_num = 1
        for line in hin:
            line_num += 1
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 12:
                # Skip malformed lines
                continue

            rs_id = fields[0]
            alleles_raw = fields[1]
            chrom = fields[2]
            pos = fields[3]
            # fields[4]..fields[10] are not used for VCF minimal output
            sample_gts_raw = fields[11:]

            # Parse REF/ALT candidates from alleles field
            allele_order = parse_alleles_field(alleles_raw)
            if allele_order:
                ref_base = allele_order[0]
                alt_candidates = allele_order[1:]
            else:
                # Fallback: infer REF from the first non-missing genotype allele we see
                ref_base = None
                alt_candidates = []

            # First pass: normalize all genotypes and collect observed alleles
            normalized_gts = []
            observed_alts = set(alt_candidates)

            for gt_raw in sample_gts_raw:
                gt = normalize_genotype_raw(gt_raw)
                normalized_gts.append(gt)
                if gt is None:
                    continue
                a1, a2 = gt
                if ref_base is None and a1 in {"A", "C", "G", "T"}:
                    ref_base = a1
                if ref_base is None and a2 in {"A", "C", "G", "T"}:
                    ref_base = a2
                # Collect ALTs distinct from REF
                if a1 in {"A", "C", "G", "T"}:
                    if ref_base is not None and a1 != ref_base:
                        observed_alts.add(a1)
                if a2 in {"A", "C", "G", "T"}:
                    if ref_base is not None and a2 != ref_base:
                        observed_alts.add(a2)

            if ref_base is None:
                # If we still cannot decide REF, default to 'N' and mark site missing
                ref_base = "N"

            # Build ALT list in deterministic order
            alt_list = [a for a in allele_order[1:] if a in observed_alts]
            for a in sorted(observed_alts):
                if a not in alt_list:
                    alt_list.append(a)

            # Create allele index mapping for GT encoding
            allele_to_index = {ref_base: 0}
            for idx, alt in enumerate(alt_list, start=1):
                allele_to_index[alt] = idx

            # If no ALT (monomorphic), VCF expects '.' in ALT column
            alt_field = ",".join(alt_list) if alt_list else "."

            # Encode genotypes
            gt_strings = []
            for gt in normalized_gts:
                if gt is None:
                    gt_strings.append("./.")
                else:
                    a1, a2 = gt
                    if a1 not in allele_to_index or a2 not in allele_to_index:
                        gt_strings.append("./.")
                    else:
                        g = f"{allele_to_index[a1]}/{allele_to_index[a2]}"
                        gt_strings.append(g)

            # Compose VCF line
            chrom_field = str(chrom)
            pos_field = str(pos)
            id_field = rs_id if rs_id != "." else "."
            ref_field = ref_base
            qual_field = "."
            filter_field = "PASS"
            info_field = "."
            format_field = "GT"

            row = [
                chrom_field,
                pos_field,
                id_field,
                ref_field,
                alt_field,
                qual_field,
                filter_field,
                info_field,
                format_field,
            ] + gt_strings
            hout.write("\t".join(row) + "\n")


def main():
    # 默认输入/输出路径：相对于本脚本所在目录的示例文件
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input = os.path.join(script_dir, "origin_maize.geno.selected.hmp.txt")
    default_output = os.path.join(script_dir, "origin_maize.geno.selected.vcf")

    # 支持命令行参数覆盖：python hapmap_to_vcf.py <input_hmp.txt> <output.vcf>
    if len(sys.argv) >= 3:
        input_path = sys.argv[1]
        output_path = sys.argv[2]
    else:
        input_path = default_input
        output_path = default_output

    if not os.path.isfile(input_path):
        print(f"[Error] 输入文件不存在: {input_path}")
        sys.exit(1)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    print(f"[Info] 读取 HapMap: {input_path}")
    print(f"[Info] 写入 VCF:   {output_path}")
    try:
        convert_hapmap_to_vcf(input_path, output_path)
    except Exception as e:
        print(f"[Error] 转换失败: {e}")
        sys.exit(1)
    print("[Done] 转换完成")


if __name__ == "__main__":
    main()


