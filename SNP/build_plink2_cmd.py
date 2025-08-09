"""
构建可运行的 plink2 PCA 命令的辅助脚本。

功能：
- 接收 plink2 的可执行文件路径，以及 --threads、--vcf、--pca、--out、--out-dir 等参数；
- 生成一条可直接在命令行运行的 plink2 命令并打印到标准输出。

用法示例：
    python build_plink2_cmd.py --plink2 ./plink2 --threads 30 --vcf *.vcf --pca 10 --out pca10 --out-dir ./results

说明：
在--plink2 参数中填入plink2.exe的路径
在--vcf 参数中填入vcf文件的路径
在--pca 参数中填入主成分个数
在--out 参数中填入输出前缀。
在--out-dir 参数中填入输出目录（可选，不填则默认当前工作目录）。
在--threads 参数中填入线程数。
输出文件为 pca10.eigenvec 和 pca10.eigenval；可通过 --out-dir 指定输出目录。
"""
import argparse
from typing import List


def needs_quotes(value: str) -> bool:
    """判断参数值是否需要加引号。

    规则：
    - 如果包含空格或常见会打断命令的特殊字符，则需要加引号；
    - 对于纯通配模式（例如 *.vcf），默认不加引号，以便在支持的 shell 中让通配符展开。
    """
    if value is None:
        return False
    special = set(' <>|&()')
    if any(ch in special for ch in value):
        return True
    return ' ' in value


def quote(value: str) -> str:
    """根据需要为参数加引号。"""
    if value is None:
        return ''
    if needs_quotes(value):
        return f'"{value}"'
    return value


def build_plink2_command(plink2_path: str, threads: int, vcf_pattern: str, pca_k: int, out_prefix: str, out_dir: str | None) -> str:
    """根据输入参数构建 plink2 命令字符串。

    参数：
    - plink2_path: plink2 可执行文件路径（可绝对/相对）；
    - threads: 线程数；
    - vcf_pattern: VCF 文件路径或通配符；
    - pca_k: 主成分个数（--pca 的 K 值）；
    - out_prefix: 输出前缀（--out）。
    """
    tokens: List[str] = []
    # 第一个 token：plink2 可执行文件路径（若包含空格会自动加引号）
    tokens.append(quote(plink2_path))
    tokens.extend(["--threads", str(threads)])
    # 对于包含 * 或 ? 的通配模式，保持原样以便在支持的 shell 中进行展开；
    # 若传入的是具体路径（可能含空格），会自动加引号。
    vcf_token = vcf_pattern if ('*' in vcf_pattern or '?' in vcf_pattern) else quote(vcf_pattern)
    tokens.extend(["--vcf", vcf_token])
    tokens.extend(["--pca", str(pca_k)])
    # 组合输出路径（目录 + 前缀）
    if out_dir:
        # 不在此处添加扩展名，由 plink2 自行追加 .eigenvec/.eigenval
        full_out = quote(out_dir.rstrip('/\\') + '/' + out_prefix)
    else:
        full_out = quote(out_prefix)
    tokens.extend(["--out", full_out])
    return " ".join(tokens)


def main():
    parser = argparse.ArgumentParser(description="根据输入参数构建一条可运行的 plink2 PCA 命令并打印。")
    parser.add_argument("--plink2", dest="plink2_path", default="./plink2", help="plink2 可执行文件路径（如 ./plink2 或 C:/path/plink2.exe）")
    parser.add_argument("--threads", dest="threads", type=int, default=30, help="线程数（--threads）")
    parser.add_argument("--vcf", dest="vcf", default="*.vcf", help="VCF 文件或通配符（--vcf），如 sample.vcf 或 *.vcf")
    parser.add_argument("--pca", dest="pca_k", type=int, default=10, help="主成分个数（--pca K）")
    parser.add_argument("--out", dest="out_prefix", default="pca10", help="输出前缀（--out）")
    parser.add_argument("--out-dir", dest="out_dir", default=None, help="输出目录（可选，将与前缀拼接为完整输出路径）")

    args = parser.parse_args()

    cmd = build_plink2_command(
        plink2_path=args.plink2_path,
        threads=args.threads,
        vcf_pattern=args.vcf,
        pca_k=args.pca_k,
        out_prefix=args.out_prefix,
        out_dir=args.out_dir,
    )

    print(cmd)


if __name__ == "__main__":
    main()


