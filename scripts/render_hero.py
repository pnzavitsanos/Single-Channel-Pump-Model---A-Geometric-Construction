"""
render_hero.py
===============

Banner image: the blade with wall thickness, standing between the two
parallel planes z = 0 and z = H.

Usage:
    python scripts/render_hero.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, draw_planes, setup_3d, trim_png
from channel_model import (DEFAULT_PARAMS as P, BASELINE_INLET_DEG, BASELINE_OUTLET_DEG,
                           channel_shell, frame)

OUT = os.path.join(MEDIA, "hero_channel.png")
RADIUS = 0.78


def main():
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    fig.patch.set_facecolor("white")

    draw_planes(ax, RADIUS, P.height)
    for sign, color in ((-1.0, "#256abf"), (+1.0, "#6da7ec")):
        x, y, z = channel_shell(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG, sign,
                                n_theta=260, n_w=20, offset="horizontal")
        ax.plot_surface(x, y, z, color=color, alpha=0.95, linewidth=0, antialiased=True, shade=True)

    th = np.linspace(np.deg2rad(BASELINE_INLET_DEG), np.deg2rad(BASELINE_OUTLET_DEG), 300)
    r = frame(th)[0]
    ax.plot(r[0], r[1], r[2], color=INK, linewidth=1.2)

    setup_3d(ax, RADIUS, P.height, elev=30, azim=-62, zoom=1.18)
    ax.set_title("Log-spiral blade between two parallel planes\n"
                 r"$X(\theta,w)=r(\theta)+\frac{w}{\cos\beta(\theta)}\,F(\theta)$,"
                 r"   $0\leq w\leq H$", fontsize=12, color=INK)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.88)
    fig.savefig(OUT, dpi=150, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
