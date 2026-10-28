"""
render_hero.py
===============

media/hero_impeller.png: (a) the single-blade impeller (hub, blade, and a
cut-away shroud) and (b) the blade alone with its straight rulings.

Usage:
    python scripts/render_hero.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import (MEDIA, INK, BLADE, draw_meshes, draw_impeller,
                        draw_planes, setup_3d, trim_png)
from impeller_geometry import NOTEBOOK, impeller_meshes, surface, grid

OUT = os.path.join(MEDIA, "hero_impeller.png")


def main(d=NOTEBOOK):
    parts = impeller_meshes(d)
    R = np.linalg.norm(parts["hub"][..., :2], axis=-1).max()

    fig = plt.figure(figsize=(13, 5.6))
    fig.patch.set_facecolor("white")

    ax = fig.add_subplot(1, 2, 1, projection="3d", computed_zorder=False)
    draw_impeller(ax, d)
    setup_3d(ax, R, d.height + 0.05, elev=32, azim=-60, zoom=1.25)
    ax.set_zlim(-0.03, d.height + 0.05)
    ax.set_title("(a) Impeller: hub, blade, and shroud (front half cut away)", color=INK, fontsize=13)

    ax = fig.add_subplot(1, 2, 2, projection="3d", computed_zorder=False)
    draw_planes(ax, R, d.height)
    draw_meshes(ax, [(parts["blade"], BLADE, 1.0)])
    th, wv, TH, W = grid(d, 25, 2)
    X = surface(TH, W, d)
    for i in range(len(th)):
        ax.plot(X[0, i], X[1, i], X[2, i], color=INK, linewidth=0.7)
    setup_3d(ax, R, d.height, elev=32, azim=-60, zoom=1.25)
    ax.set_title("(b) Blade: straight rulings from hub plane to shroud plane", color=INK, fontsize=13)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.92, wspace=0.02)
    fig.savefig(OUT, dpi=140, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
