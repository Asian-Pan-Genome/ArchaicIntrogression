#!/usr/bin/env python3
"""Draw TP1-TP4 and DP1-DP4 with one consistent visual specification."""

from __future__ import annotations

import argparse
from pathlib import Path


# Display-only values chosen to keep arrows visually separated. They are not
# fitted pulse-time estimates from the fastsimcoal2 analysis.
ILLUSTRATIVE_PULSE_TIMES = {
    1: (1200,),
    2: (950, 1500),
    3: (850, 1250, 1650),
    4: (800, 1100, 1400, 1700),
}

# Generations before present; approximately 283, 307, 363, and 420 kya
# with 29 years per generation.
DONOR_SPLITS = (9759, 10599, 12517, 14483)
ILLUSTRATIVE_TOTAL_ADMIXTURE = 0.006
MODEL_CHOICES = tuple(f"{family}{k}" for family in ("TP", "DP") for k in range(1, 5))

# Fixed across every panel. DENI and D1 are kept distinct because they denote
# different modeled populations, even though D1 is the closest DP analogue of
# the shared TP donor.
POPULATION_COLOURS = {
    "AFR": "#4E79A7",
    "EAS": "#E15759",
    "DENI": "#F28E2B",
    "D1": "#EDC948",
    "D2": "#B07AA1",
    "D3": "#59A14F",
    "D4": "#76B7B2",
    "DENS": "#9C755F",
}


def configure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.weight": "normal",
            "axes.titleweight": "normal",
            "axes.labelweight": "normal",
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
        }
    )
    import matplotlib.pyplot as plt

    return plt


def add_fixed_modern_history(demography) -> None:
    demography.add_population(name="AFR", initial_size=43477.1976)
    demography.add_population(
        name="EAS", initial_size=94703.382, growth_rate=0.0022539867
    )
    demography.add_population_parameters_change(
        time=1765, population="EAS", initial_size=726.3523, growth_rate=0
    )
    demography.add_mass_migration(time=2151, source="EAS", dest="AFR", proportion=1)
    demography.add_population_parameters_change(
        time=3640, population="AFR", initial_size=19308.0474
    )


def build_model(family: str, k: int):
    import msprime

    demography = msprime.Demography()
    add_fixed_modern_history(demography)

    if family == "TP":
        demography.add_population(
            name="DENI", initial_size=2534.7038, description="shared ghost donor"
        )
        demography.add_population(
            name="DENS", initial_size=2534.7038, description="Altai Denisovan"
        )
        for time in ILLUSTRATIVE_PULSE_TIMES[k]:
            demography.add_mass_migration(
                time=time,
                source="EAS",
                dest="DENI",
                proportion=ILLUSTRATIVE_TOTAL_ADMIXTURE / k,
            )
        demography.add_mass_migration(
            time=DONOR_SPLITS[1], source="DENS", dest="DENI", proportion=1
        )
        demography.add_mass_migration(
            time=20395, source="DENI", dest="AFR", proportion=1
        )
    else:
        for donor_index in range(1, k + 1):
            demography.add_population(
                name=f"D{donor_index}",
                initial_size=2534.7038,
                description=f"ghost donor D{donor_index}",
            )
        demography.add_population(
            name="DENS", initial_size=2534.7038, description="Altai Denisovan"
        )
        for donor_index, time in enumerate(ILLUSTRATIVE_PULSE_TIMES[k], start=1):
            demography.add_mass_migration(
                time=time,
                source="EAS",
                dest=f"D{donor_index}",
                proportion=ILLUSTRATIVE_TOTAL_ADMIXTURE / k,
            )
        demography.add_mass_migration(
            time=DONOR_SPLITS[0], source="DENS", dest="D1", proportion=1
        )
        for donor_index in range(2, k + 1):
            demography.add_mass_migration(
                time=DONOR_SPLITS[donor_index - 1],
                source=f"D{donor_index}",
                dest="D1",
                proportion=1,
            )
        demography.add_mass_migration(
            time=20395, source="D1", dest="AFR", proportion=1
        )

    demography.sort_events()
    demography.debug()
    return demography


def panel_title(family: str, k: int) -> str:
    pulse_word = "pulse" if k == 1 else "pulses"
    if family == "TP":
        return f"TP{k}: {k} {pulse_word}; shared donor (307 ka split)"
    split_text = "/".join(str(x) for x in (283, 307, 363, 420)[:k])
    donor_word = "donor" if k == 1 else "donors"
    return f"DP{k}: {k} {pulse_word}; {k} divergent {donor_word} ({split_text} ka)"


def draw_panel(ax, family: str, k: int, *, combined: bool) -> None:
    import demesdraw

    graph = build_model(family, k).to_demes()
    active_colours = {
        deme.name: POPULATION_COLOURS[deme.name] for deme in graph.demes
    }
    demesdraw.tubes(
        graph,
        ax=ax,
        log_time=True,
        colours=active_colours,
        seed=20260908,
    )
    ax.set_ylim(1, 30000)
    ax.set_yticks([1, 10, 100, 1000, 10000])
    ax.set_yticklabels(
        [r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"]
    )
    ax.set_xlabel("Population", fontsize=7.5 if combined else 11)
    ax.set_ylabel("Time ago (generations)", fontsize=7.5 if combined else 11)
    ax.set_title(
        panel_title(family, k),
        fontsize=8.1 if combined else 12,
        fontweight="normal",
        pad=4 if combined else 7,
    )
    ax.tick_params(
        axis="x",
        which="major",
        labelsize=7.0 if combined else 10,
        length=2.4 if combined else 3.5,
        pad=1.5 if combined else 3,
    )
    ax.tick_params(
        axis="y",
        which="major",
        # Mathtext exponents are rendered at about 70% of the base size.
        # A 10-point base therefore keeps the exponent near 7 points.
        labelsize=10,
        length=2.4 if combined else 3.5,
        pad=1.5 if combined else 3,
    )
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("normal")


def save_individuals(output_dir: Path, plt, selected_model: str | None = None) -> None:
    individual_dir = output_dir / "individual"
    individual_dir.mkdir(parents=True, exist_ok=True)
    for k in range(1, 5):
        for family in ("TP", "DP"):
            model = f"{family}{k}"
            if selected_model is not None and model != selected_model:
                continue
            fig, ax = plt.subplots(figsize=(7.1, 4.8))
            draw_panel(ax, family, k, combined=False)
            fig.tight_layout()
            stem = individual_dir / f"{family}{k}_model"
            fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
            fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
            plt.close(fig)


def save_combined(output_dir: Path, plt, layout: str) -> None:
    combined_dir = output_dir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)

    if layout == "paired":
        panel_order = [item for k in range(1, 5) for item in (("TP", k), ("DP", k))]
        nrows, ncols = 4, 2
        figsize = (8.27, 11.69)  # A4 portrait
        layout_tag = "2col_4row"
    else:
        panel_order = [("TP", k) for k in range(1, 5)] + [("DP", k) for k in range(1, 5)]
        nrows, ncols = 2, 4
        figsize = (11.69, 8.27)  # A4 landscape
        layout_tag = "4col_2row"

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    for ax, (family, k) in zip(axes.flat, panel_order):
        draw_panel(ax, family, k, combined=True)
    fig.suptitle(
        "Denisovan introgression models: shared versus divergent donors",
        fontsize=10.5,
        fontweight="normal",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.0, w_pad=0.8)

    stem = combined_dir / f"eight_model_schematics_{layout_tag}_A4"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--model",
        choices=MODEL_CHOICES,
        help="Generate one individual model only (for example, TP1 or DP3).",
    )
    parser.add_argument(
        "--layout",
        choices=("paired", "families"),
        default="paired",
        help="paired = TPk/DPk in each row; families = TP row followed by DP row",
    )
    args = parser.parse_args()

    plt = configure_matplotlib()
    save_individuals(args.output_dir, plt, selected_model=args.model)
    if args.model is None:
        save_combined(args.output_dir, plt, args.layout)
    print(f"Wrote figures to {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
