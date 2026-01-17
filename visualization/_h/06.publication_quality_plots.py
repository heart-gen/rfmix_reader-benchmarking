"""Publication-quality visualization workflow using rfmix-reader helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from rfmix_reader import (
    read_rfmix,
    read_flare,
    read_simu,
    plot_global_ancestry,
    plot_ancestry_by_chromosome,
    generate_tagore_bed,
    plot_local_ancestry_tagore,
)

STRUCTURE_URL = (
    "https://raw.githubusercontent.com/LieberInstitute/aanri_phase1/"
    "main/ancestry_structure/structure.out_ancestry_proportion_raceDemo_compare"
)

ABCA7_REGION = {
    "chrom": "19",
    "start": 1_039_997,
    "end": 1_065_572,
}


@dataclass(frozen=True)
class PlotPaths:
    out_dir: Path

    def for_name(self, name: str) -> str:
        return str(self.out_dir / name)


def _ensure_dataframe(data) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.DataFrame(data)


def _infer_sample_id_column(df: pd.DataFrame) -> str:
    for col in ("sample_id", "sample", "iid", "id"):
        if col in df.columns:
            return col
    return df.index.name or "sample_id"


def _add_sample_id_column(df: pd.DataFrame, sample_col: str | None = None) -> pd.DataFrame:
    working = df.copy()
    sample_col = sample_col or _infer_sample_id_column(working)
    if sample_col not in working.columns:
        working = working.reset_index().rename(columns={"index": sample_col})
    return working


def _infer_chrom_column(df: pd.DataFrame) -> str | None:
    for col in ("chrom", "chromosome", "chr"):
        if col in df.columns:
            return col
    return None


def _infer_ancestry_columns(df: pd.DataFrame, sample_col: str, chrom_col: str | None) -> list[str]:
    exclude = {sample_col}
    if chrom_col:
        exclude.add(chrom_col)
    numeric_cols = [col for col in df.columns if col not in exclude]
    return numeric_cols


def _align_by_sample_id(*frames: pd.DataFrame) -> list[pd.DataFrame]:
    sample_col = _infer_sample_id_column(frames[0])
    sample_sets = [set(frame[sample_col].astype(str)) for frame in frames]
    common_ids = set.intersection(*sample_sets)
    aligned = []
    for frame in frames:
        temp = frame.copy()
        temp[sample_col] = temp[sample_col].astype(str)
        aligned.append(temp[temp[sample_col].isin(common_ids)].sort_values(sample_col))
    return aligned


def tidy_global_ancestry(g_anc: pd.DataFrame) -> pd.DataFrame:
    """Return per-individual global ancestry with consistent sample_id column."""
    g_anc = _ensure_dataframe(g_anc)
    g_anc = _add_sample_id_column(g_anc)
    sample_col = _infer_sample_id_column(g_anc)
    chrom_col = _infer_chrom_column(g_anc)
    ancestry_cols = _infer_ancestry_columns(g_anc, sample_col, chrom_col)
    if chrom_col:
        grouped = (
            g_anc.groupby(sample_col, as_index=False)[ancestry_cols]
            .mean(numeric_only=True)
        )
        return grouped
    return g_anc[[sample_col, *ancestry_cols]].copy()


def make_simulation_global_plots(
    sim_dir: Path,
    rfmix_dir: Path,
    flare_dir: Path,
    out_dir: Path,
    label: str,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = PlotPaths(out_dir)

    loci_sim, g_anc_sim, _ = read_simu(sim_dir)
    _, g_anc_rfmix, _ = read_rfmix(rfmix_dir)
    _, g_anc_flare, _ = read_flare(flare_dir)

    g_anc_sim = _ensure_dataframe(g_anc_sim)
    g_anc_rfmix = _ensure_dataframe(g_anc_rfmix)
    g_anc_flare = _ensure_dataframe(g_anc_flare)

    g_anc_sim, g_anc_rfmix, g_anc_flare = _align_by_sample_id(
        g_anc_sim, g_anc_rfmix, g_anc_flare
    )

    plot_global_ancestry(
        g_anc_sim,
        save_path=paths.for_name(f"global_ancestry_{label}_sim"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plot_global_ancestry(
        g_anc_rfmix,
        save_path=paths.for_name(f"global_ancestry_{label}_rfmix"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plot_global_ancestry(
        g_anc_flare,
        save_path=paths.for_name(f"global_ancestry_{label}_flare"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plot_ancestry_by_chromosome(
        g_anc_sim,
        save_path=paths.for_name(f"by_chrom_{label}_sim"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plot_ancestry_by_chromosome(
        g_anc_rfmix,
        save_path=paths.for_name(f"by_chrom_{label}_rfmix"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plot_ancestry_by_chromosome(
        g_anc_flare,
        save_path=paths.for_name(f"by_chrom_{label}_flare"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )


def _load_structure_ancestry(url: str = STRUCTURE_URL) -> pd.DataFrame:
    structure_df = pd.read_csv(url, sep=None, engine="python")
    return structure_df


def _normalize_ancestry_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    working = df.copy()
    normalized = {}
    for col in working.columns:
        lower = col.lower()
        if any(key in lower for key in ("afr", "african")):
            normalized[col] = f"afr_{prefix}"
        elif any(key in lower for key in ("eur", "europe")):
            normalized[col] = f"eur_{prefix}"
        elif any(key in lower for key in ("amr", "native", "amer")):
            normalized[col] = f"amr_{prefix}"
        elif any(key in lower for key in ("eas", "east", "asian")):
            normalized[col] = f"eas_{prefix}"
    return working.rename(columns=normalized)


def _prepare_structure_df(structure_df: pd.DataFrame) -> pd.DataFrame:
    structure_df = _add_sample_id_column(structure_df)
    sample_col = _infer_sample_id_column(structure_df)
    structure_df = structure_df.rename(columns={sample_col: "sample_id"})
    structure_df = _normalize_ancestry_columns(structure_df, "structure")
    return structure_df


def compare_rfmix_structure(
    g_anc_aanri: pd.DataFrame,
    structure_df: pd.DataFrame,
    out_dir: Path,
    ancestries: Iterable[str] | None = None,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = PlotPaths(out_dir)

    rfmix_df = tidy_global_ancestry(g_anc_aanri)
    rfmix_df = rfmix_df.rename(columns={_infer_sample_id_column(rfmix_df): "sample_id"})
    rfmix_df = _normalize_ancestry_columns(rfmix_df, "rfmix")

    structure_df = _prepare_structure_df(structure_df)

    merged = rfmix_df.merge(structure_df, on="sample_id", how="inner")

    if ancestries is None:
        ancestries = [
            ancestry
            for ancestry in ("afr", "eur", "amr", "eas")
            if f"{ancestry}_rfmix" in merged.columns
            and f"{ancestry}_structure" in merged.columns
        ]

    for ancestry in ancestries:
        rfmix_col = f"{ancestry}_rfmix"
        structure_col = f"{ancestry}_structure"

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.kdeplot(
            data=merged,
            x=rfmix_col,
            ax=ax,
            label="RFMix",
            fill=True,
            alpha=0.4,
        )
        sns.kdeplot(
            data=merged,
            x=structure_col,
            ax=ax,
            label="STRUCTURE",
            fill=True,
            alpha=0.4,
        )
        ax.set_xlabel(f"{ancestry.upper()} global ancestry")
        ax.set_ylabel("Density")
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(
                f"{paths.for_name(f'aanri_{ancestry}_distribution')}" f".{ext}",
                dpi=300,
                bbox_inches="tight",
            )
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5, 5))
        sns.regplot(
            data=merged,
            x=structure_col,
            y=rfmix_col,
            ax=ax,
            scatter_kws={"s": 20, "alpha": 0.6},
            line_kws={"color": "black"},
            ci=95,
        )
        min_val = min(merged[structure_col].min(), merged[rfmix_col].min())
        max_val = max(merged[structure_col].max(), merged[rfmix_col].max())
        ax.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="gray")
        ax.set_xlabel(f"STRUCTURE {ancestry.upper()}")
        ax.set_ylabel(f"RFMix {ancestry.upper()}")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(
                f"{paths.for_name(f'aanri_{ancestry}_scatter')}" f".{ext}",
                dpi=300,
                bbox_inches="tight",
            )
        plt.close(fig)

    return merged


def plot_aanri_rfmix_overview(g_anc_aanri: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = PlotPaths(out_dir)

    plot_global_ancestry(
        g_anc_aanri,
        save_path=paths.for_name("aanri_global_ancestry_rfmix"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plot_ancestry_by_chromosome(
        g_anc_aanri,
        save_path=paths.for_name("aanri_by_chrom_rfmix"),
        save_multi_format=("png", "pdf"),
        dpi=300,
        bbox_inches="tight",
    )


def _tagore_region_filter(bed_df: pd.DataFrame) -> pd.DataFrame:
    bed_df = bed_df.copy()
    chrom_col = _infer_chrom_column(bed_df) or "chrom"
    start_col = "start" if "start" in bed_df.columns else "chromStart"
    end_col = "end" if "end" in bed_df.columns else "chromEnd"
    bed_df[chrom_col] = bed_df[chrom_col].astype(str)
    mask = (
        (bed_df[chrom_col] == ABCA7_REGION["chrom"])
        & (bed_df[start_col] <= ABCA7_REGION["end"])
        & (bed_df[end_col] >= ABCA7_REGION["start"])
    )
    return bed_df.loc[mask]


def _pick_sample_ids_by_iqr(
    g_anc: pd.DataFrame,
    ancestry_col: str,
    chrom: str,
) -> dict[str, str]:
    sample_col = _infer_sample_id_column(g_anc)
    g_anc = g_anc.copy()
    chrom_col = _infer_chrom_column(g_anc)
    if chrom_col:
        g_anc = g_anc[g_anc[chrom_col].astype(str) == chrom]

    stats = g_anc[ancestry_col].describe()
    low_threshold = stats["25%"]
    high_threshold = stats["75%"]

    low = g_anc[g_anc[ancestry_col] <= low_threshold].sort_values(ancestry_col).head(1)
    high = g_anc[g_anc[ancestry_col] >= high_threshold].sort_values(ancestry_col).tail(1)
    mid = g_anc[
        (g_anc[ancestry_col] > low_threshold)
        & (g_anc[ancestry_col] < high_threshold)
    ].sort_values(ancestry_col).head(1)

    return {
        "low": low[sample_col].iloc[0],
        "medium": mid[sample_col].iloc[0],
        "high": high[sample_col].iloc[0],
    }


def _build_tagore_bed(
    loci_df: pd.DataFrame,
    g_anc: pd.DataFrame,
    admix,
    sample_id,
) -> pd.DataFrame:
    return generate_tagore_bed(loci_df, g_anc, admix, sample_id)


def plot_abca7_local_ancestry(
    loci_df: pd.DataFrame,
    g_anc: pd.DataFrame,
    admix,
    out_dir: Path,
    sample_ids: dict[str, str],
    prefix_label: str,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for label, sample_id in sample_ids.items():
        bed_df = _build_tagore_bed(loci_df, g_anc, admix, sample_id)
        bed_locus = _tagore_region_filter(bed_df)
        prefix = out_dir / f"abca7_local_ancestry_{prefix_label}_{label}"
        for ext in ("png", "pdf"):
            plot_local_ancestry_tagore(
                bed_df=bed_locus,
                prefix=str(prefix),
                build="hg38",
                oformat=ext,
                verbose=True,
                force=True,
            )


def maybe_plot_flare_abca7(
    flare_dir: Path,
    out_dir: Path,
    sample_ids: dict[str, str],
) -> None:
    if not flare_dir.exists():
        return
    loci_flare, g_anc_flare, admix_flare = read_flare(flare_dir)
    plot_abca7_local_ancestry(
        loci_flare,
        _add_sample_id_column(_ensure_dataframe(g_anc_flare)),
        admix_flare,
        out_dir,
        sample_ids,
        prefix_label="flare",
    )


def main():
    base_dir = Path("input")
    figures_dir = Path("figures")

    make_simulation_global_plots(
        sim_dir=base_dir / "simulations/two_populations",
        rfmix_dir=base_dir / "simulations/two_populations/_m/rfmix-out",
        flare_dir=base_dir / "simulations/two_populations/_m/flare-out",
        out_dir=figures_dir / "simulations/two_pop",
        label="two_pop",
    )
    make_simulation_global_plots(
        sim_dir=base_dir / "simulations/three_populations",
        rfmix_dir=base_dir / "simulations/three_populations/_m/rfmix-files",
        flare_dir=base_dir / "simulations/three_populations/_m/flare-out",
        out_dir=figures_dir / "simulations/three_pop",
        label="three_pop",
    )

    loci_aanri, g_anc_aanri, admix_aanri = read_rfmix(
        base_dir / "aanri_data/rfmix-version/_m"
    )

    plot_aanri_rfmix_overview(g_anc_aanri, figures_dir / "aanri")

    structure_df = _load_structure_ancestry()
    compare_rfmix_structure(
        g_anc_aanri,
        structure_df,
        figures_dir / "aanri",
    )

    g_anc_aanri = _add_sample_id_column(_ensure_dataframe(g_anc_aanri))
    ancestry_cols = _infer_ancestry_columns(
        g_anc_aanri, _infer_sample_id_column(g_anc_aanri), _infer_chrom_column(g_anc_aanri)
    )
    afr_col = next((col for col in ancestry_cols if "afr" in col.lower()), ancestry_cols[0])
    sample_ids = _pick_sample_ids_by_iqr(g_anc_aanri, afr_col, chrom="19")
    plot_abca7_local_ancestry(
        loci_aanri,
        g_anc_aanri,
        admix_aanri,
        figures_dir / "aanri",
        sample_ids,
        prefix_label="rfmix",
    )
    maybe_plot_flare_abca7(
        base_dir / "aanri_data/flare-version/_m",
        figures_dir / "aanri",
        sample_ids,
    )


if __name__ == "__main__":
    main()
