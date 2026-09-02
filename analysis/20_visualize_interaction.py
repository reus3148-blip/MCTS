"""Figures for the two-dimensional sensitivity grids (fig29).

A one-at-a-time tornado says which assumption matters most. It cannot show
whether an assumption's effect *depends* on another, and that is exactly what a
surface plot makes obvious: parallel contours mean the axes act independently,
tilted or twisted ones mean they interact.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REPORT_DIR = ROOT / "reports" / "interaction-sensitivity-v0.8"
TABLE_DIR = REPORT_DIR / "tables"
FIGURE_DIR = REPORT_DIR / "figures"
PUBLIC_DIR = ROOT / "public" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PUBLIC_DIR.mkdir(parents=True, exist_ok=True)

matplotlib.rcParams.update({
    "font.family": "Malgun Gothic",
    "axes.unicode_minus": False,
    "figure.facecolor": "#fafaf9",
    "axes.facecolor": "#fafaf9",
    "axes.edgecolor": "#292524",
    "axes.labelcolor": "#292524",
    "axes.titlecolor": "#1c1917",
    "xtick.color": "#57534e",
    "ytick.color": "#57534e",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
})

INK = "#1c1917"

AXIS_LABELS = {
    "reward.recurrence_free_year": "무재발 1년의 가치",
    "reward.acute_toxicity_penalty": "급성 독성 페널티",
    "toxicity.chemo.intensified": "강화항암 독성확률",
    "discount_rate_annual": "연 할인율",
}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig29_interaction_surfaces() -> None:
    detail = pd.read_csv(TABLE_DIR / "interaction_grid.csv")
    metrics = json.loads((REPORT_DIR / "metrics.json").read_text(encoding="utf-8"))
    summaries = {row["grid"]: row for row in metrics["grids"]}
    grids = list(summaries)

    # One diverging scale across all three panels, centred on zero, so a cell's
    # colour means the same thing everywhere and sign changes are visible.
    limit = float(detail["utility_gap"].abs().max())
    cmap = matplotlib.colormaps["RdYlGn"]
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

    fig, axes = plt.subplots(1, len(grids), figsize=(14.5, 4.6))
    for ax, name in zip(np.atleast_1d(axes), grids):
        block = detail[detail["grid"] == name]
        x_values = sorted(block["x_value"].unique())
        y_values = sorted(block["y_value"].unique())
        surface = np.array([
            [float(block[(block.x_value == x) & (block.y_value == y)]
                   ["utility_gap"].iloc[0]) for y in y_values]
            for x in x_values
        ])
        ax.imshow(surface, cmap=cmap, norm=norm, origin="lower", aspect="auto")
        for i in range(len(x_values)):
            for j in range(len(y_values)):
                ax.text(j, i, f"{surface[i, j]:+.3f}", ha="center", va="center",
                        fontsize=9.5, color=INK)
        summary = summaries[name]
        ax.set_xticks(range(len(y_values)))
        ax.set_xticklabels([f"{v:g}" for v in y_values])
        ax.set_yticks(range(len(x_values)))
        ax.set_yticklabels([f"{v:g}" for v in x_values])
        ax.set_xlabel(AXIS_LABELS.get(summary["y_parameter"], summary["y_parameter"]))
        ax.set_ylabel(AXIS_LABELS.get(summary["x_parameter"], summary["x_parameter"]))
        did = summary["corner_difference_in_differences"]
        se = summary["corner_did_standard_error"]
        z = summary["corner_did_z"]
        # Say whether the interaction clears seed noise, not just how big it is.
        verdict = "잡음 초과" if abs(z) >= 2 else "잡음 수준"
        ax.set_title(
            f"{summary['title']}\n{summary['kind']}\n"
            f"상호작용 {did:+.4f} ± {se:.4f}  (z={z:+.2f}, {verdict})",
            fontsize=10.5, pad=10)

    # Panel titles run to three lines, so the suptitle needs real clearance and
    # the colour bar has to sit below them rather than beside.
    fig.subplots_adjust(top=0.70, bottom=0.22)
    mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    bar = fig.colorbar(mappable, ax=np.atleast_1d(axes).tolist(),
                       orientation="horizontal", fraction=0.05, pad=0.22,
                       aspect=55)
    bar.set_label("효용 격차 (MCTS - NCCN)")
    fig.suptitle("가정은 따로 놀지 않는다 — 한 번에 하나씩 흔들면 못 보는 것",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig29_interaction_surfaces.png")


def main() -> None:
    fig29_interaction_surfaces()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
