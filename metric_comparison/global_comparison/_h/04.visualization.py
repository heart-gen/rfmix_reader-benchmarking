def add_identity_line(ax, data_min, data_max):
    ax.plot([data_min, data_max], [data_min, data_max],
            linestyle="--", color="gray", linewidth=1)


def plot_scatter_per_ancestry(df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    comparisons = [
        ("simu", "rfmix", "Simulated vs RFMix"),
        ("simu", "flare", "Simulated vs FLARE"),
        ("flare", "rfmix", "FLARE vs RFMix"),
    ]
    method_labels = {
        "simu": "Simulated",
        "rfmix": "RFMix",
        "flare": "FLARE"
    }
    frames = []
    for left, right, label in comparisons:
        left_df = df[df["method"] == left]
        right_df = df[df["method"] == right]
        merged = left_df.merge(
            right_df, on=["sample_id", "chrom", "ancestry"],
            how="inner", suffixes=("_x", "_y"),
        )
        merged = merged.rename(columns={"g_anc_x": "x", "g_anc_y": "y"})
        merged.loc[:, "comparison"] = label
        merged["x_label"] = method_labels[left]
        merged["y_label"] = method_labels[right]
        frames.append(merged)
    plot_df = pd.concat(frames, ignore_index=True)

    g = sns.FacetGrid(
        plot_df, row="ancestry", col="comparison", hue="ancestry", height=3.4,
        aspect=1.1, sharex=False, sharey=False, margin_titles=True,
    )
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6, s=18)
    g.map_dataframe(sns.regplot, x="x", y="y", scatter=False,
                    ci=95, color="black")
    for (facet_keys, data), ax in zip(g.facet_data(), g.axes.flat):
        data_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        data_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        add_identity_line(ax, data_min, data_max)
        # Set dynamic axis labels
        x_label = data["x_label"].iloc[0]
        y_label = data["y_label"].iloc[0]
        ax.set_xlabel(x_label, fontweight="bold")
        ax.set_ylabel(y_label, fontweight="bold")
    # Save figure
    png_path = output_dir / "global_ancestry_scatter_by_ancestry.png"
    pdf_path = output_dir / "global_ancestry_scatter_by_ancestry.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def plot_scatter_by_comparison(df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    comparisons = [
        ("simu", "rfmix", "Simulated vs RFMix"),
        ("simu", "flare", "Simulated vs FLARE"),
        ("flare", "rfmix", "FLARE vs RFMix"),
    ]
    frames = []
    for left, right, label in comparisons:
        left_df = df[df["method"] == left]
        right_df = df[df["method"] == right]
        merged = left_df.merge(
            right_df, on=["sample_id", "chrom", "ancestry"],
            how="inner", suffixes=("_x", "_y"),
        )
        merged = merged.rename(columns={"g_anc_x": "x", "g_anc_y": "y"})
        merged.loc[:, "comparison"] = label
        frames.append(merged)
    plot_df = pd.concat(frames, ignore_index=True)

    g = sns.FacetGrid(
        plot_df, col="comparison", hue="ancestry", height=4.0,
        aspect=1.1, sharex=False, sharey=False,
    )
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6, s=20)
    g.map_dataframe(sns.regplot, x="x", y="y", scatter=False,
                    ci=95, color="black")
    for ax in g.axes.flat:
        data_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        data_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        add_identity_line(ax, data_min, data_max)
        ax.set_xlabel("", fontweight="bold")
        ax.set_ylabel("", fontweight="bold")
    g.set_axis_labels("Global Ancestry (x)", "Global Ancestry (y)")
    g.add_legend(title="Ancestry", fontsize=10, title_fontsize=11)

    png_path = output_dir / "global_ancestry_scatter_by_comparison.png"
    pdf_path = output_dir / "global_ancestry_scatter_by_comparison.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def plot_boxplots(summary_df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    plot_df = summary_df.copy()
    plot_df = plot_df.rename(
        columns={"comparison": "Comparison", "ancestry": "Ancestry"}
    )
    plot_df = plot_df.melt(
        id_vars=["Comparison", "Ancestry"], value_vars=["corr", "bias"],
        var_name="Metric", value_name="Value",
    )

    g = sns.catplot(
        data=plot_df, x="Comparison", y="Value", hue="Ancestry",
        col="Metric", kind="box", sharey=False, height=4.0,
        aspect=1.1,
    )

    for ax in g.axes.flat:
        sns.stripplot(
            data=plot_df[plot_df["Metric"] == ax.get_title().split(" = ")[-1]],
            x="Comparison", y="Value", hue="Ancestry", dodge=True,
            alpha=0.6, size=4, ax=ax,
        )
        ax.set_xlabel("", fontweight="bold")
        ax.set_ylabel("", fontweight="bold")
        ax.legend_.remove()
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment('right')

    handles, labels = g.axes.flat[0].get_legend_handles_labels()
    g.fig.legend(handles, labels, title="Ancestry", loc="upper right")

    png_path = output_dir / "global_ancestry_metrics_boxplot.png"
    pdf_path = output_dir / "global_ancestry_metrics_boxplot.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")

