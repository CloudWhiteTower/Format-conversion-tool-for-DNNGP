r"""
高性能 TSV/CSV/eigenvec -> PKL 转换脚本（建议在 DNNGP 环境运行，保证 pandas 与模型版本一致）。

特性：
- 自动或手动指定输入格式（auto/eigenvec/csv/tsv）；
- eigenvec 支持含表头（FID IID PC1 ...），可选丢弃 FID/#FID 列；
- 将 PCA/数值列可选下采样为 float32，减少内存/加快序列化；
- 优先使用 pyarrow 引擎读取 CSV/TSV（如可用），否则使用 C 引擎；
- 使用更高的 pickle 协议保存（旧版自动回退）。

示例：
  python tsv2pkl.py \
    --input C:\\Users\\Cloud\\Desktop\\showANDdnngp\\Alignment\\pca250.aligned.eigenvec \
    --output C:\\Users\\Cloud\\Desktop\\showANDdnngp\\Alignment\\pca250.pkl \
    --format eigenvec --drop-fid --downcast-floats

  python tsv2pkl.py \
    --input C:\\Users\\Cloud\\Desktop\\showANDdnngp\\pheno\\maize.hybrid.train_phe.tsv \
    --output C:\\Users\\Cloud\\Desktop\\showANDdnngp\\pheno\\maize.hybrid.train_phe.pkl \
    --format tsv --index-col 0 --downcast-floats
"""

import argparse
import os
import sys
from typing import Optional

import pandas as pd


def choose_engine_for_sep(sep: str) -> Optional[str]:
    """尝试使用更快的读取引擎（pyarrow），否则回退 C 引擎。"""
    try:
        import pyarrow  # noqa: F401
        if sep in {",", "\t"}:
            return "pyarrow"
    except Exception:
        pass
    return "c"


def try_pickle(df: pd.DataFrame, path: str) -> None:
    """以尽量高协议写出 PKL，兼容旧版 pandas。"""
    try:
        import pickle
        protocol = max(4, pickle.HIGHEST_PROTOCOL)
        df.to_pickle(path, protocol=protocol)
    except TypeError:
        df.to_pickle(path)


def load_dataframe(
    input_path: str,
    fmt: str,
    index_col: Optional[str],
    index_col_pos: Optional[int],
    drop_fid: bool,
    downcast_floats: bool,
) -> pd.DataFrame:
    input_path = os.path.abspath(input_path)
    fmt = fmt.lower()

    # 自动推断格式
    if fmt == "auto":
        lower = input_path.lower()
        if lower.endswith(".eigenvec") or "eigenvec" in lower:
            fmt = "eigenvec"
        elif lower.endswith(".csv"):
            fmt = "csv"
        else:
            fmt = "tsv"

    if fmt == "eigenvec":
        # 空白分隔，可能含表头（FID IID PC1 ...）
        df = pd.read_csv(
            input_path,
            delim_whitespace=True,
            header=0,
            engine="c",
        )
        # 兼容旧数据集名 '#FID'
        fid_like_cols = [c for c in df.columns if c.upper() in {"FID", "#FID"}]
        iid_like_cols = [c for c in df.columns if c.upper() in {"IID", "#IID"}]

        # 选择索引列：优先 IID；否则根据用户指定；再否则根据位置
        id_col_name: Optional[str] = None
        if iid_like_cols:
            id_col_name = iid_like_cols[0]
        elif index_col is not None and index_col in df.columns:
            id_col_name = index_col
        elif index_col_pos is not None and 0 <= index_col_pos < len(df.columns):
            id_col_name = df.columns[index_col_pos]
        else:
            # 常见：['FID','IID','PC1',...] 或 ['IID','PC1',...]
            if len(df.columns) >= 2 and df.columns[1].upper() in {"IID", "#IID"}:
                id_col_name = df.columns[1]
            else:
                id_col_name = df.columns[0]

        # 丢弃 FID/#FID（若需要）并确保索引列仍然存在
        if drop_fid and fid_like_cols:
            df = df.drop(columns=fid_like_cols)
            if id_col_name in fid_like_cols or id_col_name not in df.columns:
                # 原先选定的索引列被删除，重新选择可用的 ID 列
                if iid_like_cols:
                    for c in iid_like_cols:
                        if c in df.columns:
                            id_col_name = c
                            break
                if id_col_name not in df.columns:
                    if index_col is not None and index_col in df.columns:
                        id_col_name = index_col
                    elif index_col_pos is not None and 0 <= index_col_pos < len(df.columns):
                        id_col_name = df.columns[index_col_pos]
                    else:
                        id_col_name = df.columns[0]

        # 将以 PC 开头的列尽量转为 float32
        pc_cols = [c for c in df.columns if c.upper().startswith("PC")]
        for c in pc_cols:
            new_col = pd.to_numeric(df[c], errors="coerce")
            df[c] = new_col.astype("float32")

        df = df.set_index(id_col_name)
        return df

    elif fmt in {"csv", "tsv"}:
        sep = "," if fmt == "csv" else "\t"
        engine = choose_engine_for_sep(sep)
        df = pd.read_csv(
            input_path,
            sep=sep,
            header=0,
            engine=engine,
        )

        # 索引列：优先列名、其次位置、最后默认第 0 列
        if index_col is not None and index_col in df.columns:
            id_col_name = index_col
        elif index_col_pos is not None and 0 <= index_col_pos < len(df.columns):
            id_col_name = df.columns[index_col_pos]
        else:
            id_col_name = df.columns[0]

        if downcast_floats:
            float_cols = df.select_dtypes(include=["float64", "float32"]).columns
            if len(float_cols) > 0:
                df[float_cols] = df[float_cols].astype("float32")

        df = df.set_index(id_col_name)
        return df

    else:
        raise ValueError("不支持的输入格式：" + fmt)


def main():
    parser = argparse.ArgumentParser(description="高性能 TSV/CSV/eigenvec -> PKL 转换")
    parser.add_argument("--input", required=True, help="输入文件路径（tsv/csv/eigenvec）")
    parser.add_argument("--output", required=True, help="输出 PKL 路径")
    parser.add_argument("--format", choices=["auto", "eigenvec", "csv", "tsv"], default="auto", help="输入格式（默认 auto 自动检测）")
    parser.add_argument("--index-col", dest="index_col", default=None, help="指定索引列名（优先级高于 --index-col-pos）")
    parser.add_argument("--index-col-pos", dest="index_col_pos", type=int, default=None, help="指定索引列位置（0 起始）")
    parser.add_argument("--drop-fid", action="store_true", help="eigenvec 输入时丢弃 FID/#FID 列")
    parser.add_argument("--downcast-floats", action="store_true", help="将数值列下采样为 float32 以减少体积")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"[Error] 输入文件不存在: {args.input}")
        sys.exit(1)

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"[Info] 读取: {args.input}")
    print(f"[Info] 格式: {args.format}")
    try:
        df = load_dataframe(
            input_path=args.input,
            fmt=args.format,
            index_col=args.index_col,
            index_col_pos=args.index_col_pos,
            drop_fid=args.drop_fid,
            downcast_floats=args.downcast_floats,
        )
    except Exception as e:
        print(f"[Error] 读取失败: {e}")
        sys.exit(2)

    print(f"[Info] 行数: {len(df):,} 列数: {len(df.columns):,}")
    try:
        try_pickle(df, args.output)
    except Exception as e:
        print(f"[Error] 保存失败: {e}")
        sys.exit(3)

    print(f"[Done] 已保存 PKL: {args.output}")


if __name__ == "__main__":
    main()

