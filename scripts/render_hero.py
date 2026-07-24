"""
render_hero.py
===============

Renders a single high-quality static image of the channel at its baseline
inlet/outlet angles, with wall thickness, for use as the README banner.

Usage:
    python scripts/render_hero.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from channel_model import channel_shell, BASELINE_INLET_DEG, BASELINE_OUTLET_DEG

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "media", "hero_channel.png")


def main():
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    n_theta, n_w = 220, 18
    xt, yt, zt = channel_shell(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG, +1.0,
                                n_theta=n_theta, n_w=n_w)
    xb, yb, zb = channel_shell(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG, -1.0,
                                n_theta=n_theta, n_w=n_w)

    ax.plot_surface(xt, yt, zt, color="#e8823c", alpha=0.75, linewidth=0, antialiased=True)
    ax.plot_surface(xb, yb, zb, color="#c23b32", alpha=0.75, linewidth=0, antialiased=True)

    ax.set_box_aspect((1, 1, 0.35))
    ax.set_axis_off()
    ax.view_init(elev=34, azim=-58)
    ax.set_title(
        "Single-channel pump skeleton — log-spiral ribbon, finite wall thickness\n"
        f"inlet = {BASELINE_INLET_DEG:.0f}°   outlet = {BASELINE_OUTLET_DEG:.0f}°",
        fontsize=12,
    )

    fig.tight_layout()
    fig.savefig(OUT, dpi=160, facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
