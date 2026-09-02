"""Figure for the utility-interaction confirmation (fig31).

v0.8 reported that the two reward weights reverse each other's direction. v0.9
re-ran the same grid with four times the seeds and the reversal did not survive.
This figure puts the two grids side by side, because the correction is only
convincing if you can see which cells moved and by how much.
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

REPORT_DIR = ROOT / "reports" / "utility-interaction-v0.9"
TABLE_DIR = REPORT_DIR / "tables"
FIGURE_DIR = REPORT_DIR / "figures"
PRIOR_DIR = ROOT / "reports" / "interaction-sensitivity-v0.8"
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
GRAY = "#78716c"
GREEN = "#15803d"
RED = "#dc2626"

METRICS = json.loads((REPORT_DIR / "metrics.json").read_text(encoding="utf-8"))


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig31_confirmation() -> None:
    new = pd.read_csv(TABLE_DIR / "confirmation_grid.csv")
    old = pd.read_csv(PRIOR_DIR / "tables" / "interaction_grid.csv")
    old = old[old["grid"] == "judgement_x_judgement"]

    rewards = sorted(new["recurrence_free_year"].unique())
    penalties = sorted(new["acute_toxicity_penalty"].unique())
    old_surface = np.array([[float(old[(old.x_value == r) & (old.y_value == p)]
                                   ["utility_gap"].iloc[0]) for p in penalties]
                            for r in rewards])
    new_surface = np.array([[float(new[(new.recurrence_free_year == r)
                                       & (new.acute_toxicity_penalty == p)]
                                   ["utility_gap"].iloc[0]) for p in penalties]
                            for r in rewards])

    limit = float(max(np.abs(old_surface).max(), np.abs(new_surface).max()))
    cmap = matplotlib.colormaps["RdYlGn"]
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

    fig = plt.figure(figsize=(13.5, 4.8))
    grid = fig.add_gridspec(1, 3, width_ratios=(1, 1, 1.15), wspace=0.32)

    for column, (surface, label, seeds) in enumerate((
        (old_surface, "v0.8", 3), (new_surface, "v0.9", 12),
    )):
        ax = fig.add_subplot(grid[0, column])
        ax.imshow(surface, cmap=cmap, norm=norm, origin="lower", aspect="auto")
        for i in range(len(rewards)):
            for j in range(len(penalties)):
                ax.text(j, i, f"{surface[i, j]:+.3f}", ha="center", va="center",
                        fontsize=10, color=INK)
        ax.set_xticks(range(len(penalties)))
        ax.set_xticklabels([f"{p:g}" for p in penalties])
        ax.set_yticks(range(len(rewards)))
        ax.set_yticklabels([f"{r:g}" for r in rewards])
        ax.set_xlabel("급성 독성 페널티")
        if column == 0:
            ax.set_ylabel("무재발 1년의 가치")
        ax.set_title(f"{label} — 시드 {seeds}개", fontsize=12, pad=10)

    ax = fig.add_subplot(grid[0, 2])
    old_did = METRICS["v0_8_reference"]["corner_difference_in_differences"]
    old_se = METRICS["v0_8_reference"]["standard_error"]
    new_did = METRICS["corner_difference_in_differences"]
    new_se = METRICS["paired_standard_error"]
    for position, (value, se, color, label) in enumerate((
        (old_did, old_se, GRAY, f"v0.8 (3시드)  z={METRICS['v0_8_reference']['z']:+.2f}"),
        (new_did, new_se, GREEN, f"v0.9 (12시드)  z={METRICS['paired_z']:+.2f}"),
    )):
        ax.errorbar([value], [position], xerr=[[1.96 * se], [1.96 * se]], fmt="o",
                    color=color, markersize=9, capsize=5, linewidth=2.0)
        ax.text(value, position + 0.22, f"{value:+.4f}", ha="center",
                fontsize=10, color=INK)
        ax.text(value, position - 0.30, label, ha="center", fontsize=9, color=color)
    ax.axvline(0, color=INK, linewidth=1.2)
    ax.set_yticks([])
    ax.set_ylim(-0.65, 1.6)
    ax.set_xlabel("상호작용 (차이-의-차이), 95% CI")
    ax.set_title("시드를 4배로 늘리자 효과가 절반이 됐다",
                 fontsize=12, pad=10, color=RED)

    fig.suptitle("v0.9 — 우리가 한 턴 전에 찾았다고 한 '방향 반전'은 살아남지 못했다",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig31_interaction_confirmation.png")


def main() -> None:
    fig31_confirmation()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
