import argparse
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from pathlib import Path

DEFAULT_SUPPLEMENT = here("input/gwas_ad/41588_2022_1024_MOESM4_ESM.xlsx")
DEFAULT_GWAS = here("input/gwas_ad/35379992-GCST90027158-MONDO_0004975.h.tsv.gz")
DEFAULT_GENES = here("functional_results/ad_local_ancestry/_h/ad_gene_intervals_grch38.tsv")


def norm_col(name: str) -> str:
    return "".join(ch.lower() for ch in str(name) if ch.isalnum())


def find_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    by_norm = {norm_col(c): c for c in df.columns}
    for cand in candidates:
        if norm_col(cand) in by_norm:
            return by_norm[norm_col(cand)]
    for col in df.columns:
        n = norm_col(col)
        if any(norm_col(cand) in n for cand in candidates):
            return col
    if required:
        raise ValueError(f"Could not find any of columns {candidates} in {list(df.columns)}")
    return None


def clean_chrom(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace("^chr", "", regex=True).str.replace(".0$", "", regex=True)


def read_gwas(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    chrom_col = find_col(df, ["chromosome", "chrom", "hm_chrom", "chr"])
    pos_col = find_col(df, ["base_pair_location", "hm_pos", "position", "pos"])
    p_col = find_col(df, ["p_value", "pvalue", "p", "pval"])
    id_col = find_col(df, ["variant_id", "rsid", "rs_id", "snp", "hm_rsid"], required=False)
    out = pd.DataFrame({
        "chrom": clean_chrom(df[chrom_col]),
        "pos": pd.to_numeric(df[pos_col], errors="coerce"),
        "p_value": pd.to_numeric(df[p_col], errors="coerce"),
    })
    out["variant_id"] = df[id_col].astype(str) if id_col else [f"chr{c}:{p}" for c, p in zip(out["chrom"], out["pos"])]
    out = out.dropna(subset=["chrom", "pos", "p_value"]).copy()
    out["pos"] = out["pos"].astype(np.int64)
    out["source"] = "gwas"
    return out


def read_table5(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "p_value", "source"])
    rows = []
    sheets = pd.read_excel(path, sheet_name=None)
    for sheet_name, df in sheets.items():
        if df.empty:
            continue
        try:
            id_col = find_col(df, ["variant_id", "rsid", "rs_id", "snp", "marker"], required=False)
            chrom_col = find_col(df, ["chromosome", "chrom", "chr"], required=False)
            pos_col = find_col(df, ["base_pair_location", "position", "pos", "bp"], required=False)
            p_col = find_col(df, ["p_value", "pvalue", "p", "pval"], required=False)
        except ValueError:
            continue
        if not (id_col or (chrom_col and pos_col)):
            continue
        tmp = pd.DataFrame()
        tmp["variant_id"] = df[id_col].astype(str) if id_col else pd.NA
        tmp["chrom"] = clean_chrom(df[chrom_col]) if chrom_col else pd.NA
        tmp["pos"] = pd.to_numeric(df[pos_col], errors="coerce") if pos_col else pd.NA
        tmp["p_value"] = pd.to_numeric(df[p_col], errors="coerce") if p_col else pd.NA
        tmp["source"] = f"supplement_table5:{sheet_name}"
        rows.append(tmp)
    if not rows:
        return pd.DataFrame(columns=["variant_id", "chrom", "pos", "p_value", "source"])
    out = pd.concat(rows, ignore_index=True).dropna(how="all", subset=["variant_id", "chrom", "pos"])
    out = out.drop_duplicates(subset=["variant_id", "chrom", "pos"])
    return out


def load_gene_windows(path: str, expand_bp: int) -> pd.DataFrame:
    genes = pd.read_csv(path, sep="\t")
    genes["chrom"] = clean_chrom(genes["chrom"])
    genes["window_start"] = (genes["start"].astype(int) - expand_bp).clip(lower=1)
    genes["window_end"] = genes["end"].astype(int) + expand_bp
    return genes


def annotate_gene_targets(gwas: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    hits = []
    for rec in windows.to_dict("records"):
        mask = (
            (gwas["chrom"].astype(str) == str(rec["chrom"]))
            & (gwas["pos"] >= rec["window_start"])
            & (gwas["pos"] <= rec["window_end"])
        )
        tmp = gwas.loc[mask].copy()
        if tmp.empty:
            continue
        tmp["gene"] = rec["gene"]
        tmp["source"] = tmp["source"] + ";gene_window"
        hits.append(tmp)
    if not hits:
        return pd.DataFrame(columns=[*gwas.columns, "gene"])
    return pd.concat(hits, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description="Prepare Bellenguez AD target SNPs for MSP local ancestry extraction")
    parser.add_argument("--gwas", default=DEFAULT_GWAS)
    parser.add_argument("--supplement", default=DEFAULT_SUPPLEMENT)
    parser.add_argument("--genes", default=DEFAULT_GENES)
    parser.add_argument("--outdir", default=here("functional_results/ad_local_ancestry/_m"))
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--gene-window-bp", type=int, default=100_000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gwas = read_gwas(args.gwas)
    windows = load_gene_windows(args.genes, args.gene_window_bp)
    table5 = read_table5(args.supplement)

    sig = gwas.loc[gwas["p_value"] < args.p_threshold].copy()
    sig["source"] = sig["source"] + ";p_lt_5e-8"
    gene_hits = annotate_gene_targets(gwas, windows)

    targets = pd.concat([sig, table5, gene_hits], ignore_index=True, sort=False)
    targets["chrom"] = clean_chrom(targets["chrom"])
    targets["pos"] = pd.to_numeric(targets["pos"], errors="coerce")
    targets = targets.dropna(subset=["chrom", "pos"]).copy()
    targets["pos"] = targets["pos"].astype(np.int64)
    targets["target_id"] = [f"ad_target_{i + 1}" for i in range(len(targets))]
    targets = targets.drop_duplicates(subset=["chrom", "pos", "variant_id"], keep="first")
    targets = targets.sort_values(["chrom", "pos", "variant_id"]).reset_index(drop=True)

    windows.to_csv(outdir / "ad_gene_windows.tsv", sep="\t", index=False)
    targets.to_csv(outdir / "ad_target_snps.tsv", sep="\t", index=False)
    print(f"Wrote {len(targets)} targets to {outdir / 'ad_target_snps.tsv'}")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
