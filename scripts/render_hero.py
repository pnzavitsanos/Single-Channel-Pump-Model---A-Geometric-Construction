"""
render_hero.py
===============

media/hero.png: (a) the pump: bottom plate, blade, and part of the top plate
(the rest is cut away to show the blade), and (b) the blade alone between the
two planes, with its straight rulings.

Usage:
    python scripts/render_hero.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, BLADE, draw_meshes, draw_pump, draw_planes, setup_3d, trim_png
from blade_geometry import NOTEBOOK, blade_mesh, surface, grid

OUT = os.path.join(MEDIA, "hero.png")


def main(b=NOTEBOOK):
    fig = plt.figure(figsize=(13, 5.6))
    fig.patch.set_facecolor("white")

    ax = fig.add_subplot(1, 2, 1, projection="3d", computed_zorder=False)
    R = draw_pump(ax, b)
    setup_3d(ax, R, b.height + 0.05, elev=32, azim=-60, zoom=1.25)
    ax.set_zlim(-0.03, b.height + 0.05)
    ax.set_title("(a) Blade between two plates (top plate partly cut away)", color=INK, fontsize=13)

    ax = fig.add_subplot(1, 2, 2, projection="3d", computed_zorder=False)
    draw_planes(ax, R, b.height)
    draw_meshes(ax, [(blade_mesh(b), BLADE, 1.0)])
    th, wv, TH, W = grid(b, 25, 2)
    X = surface(TH, W, b)
    for i in range(len(th)):
        ax.plot(X[0, i], X[1, i], X[2, i], color=INK, linewidth=0.7)
    setup_3d(ax, R, b.height, elev=32, azim=-60, zoom=1.25)
    ax.set_title("(b) Straight rulings from the plane z = 0 to the plane z = H", color=INK, fontsize=13)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.92, wspace=0.02)
    fig.savefig(OUT, dpi=140, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
