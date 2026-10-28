"""
render_fibre.py
================

media/fibre_and_spiral.png: the two ingredients of the construction and how
they combine.

  (a) a circle carrying a fibre that turns at half speed: the Moebius strip;
  (b) a spiral is a circle whose radius depends on the angle;
  (c) a spiral carrying a fibre that turns slowly, stretched to two planes:
      the blade.

Usage:
    python scripts/render_fibre.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

from plot_utils import MEDIA, INK, MUTED, BLADE, ACCENT, BLUES, draw_meshes, draw_planes, setup_3d, trim_png
from blade_geometry import (NOTEBOOK, Circle, LogSpiral, ArchimedeanSpiral, base_frame, fibre,
                            mobius_strip, surface, grid, ruling_angle)

OUT = os.path.join(MEDIA, "fibre_and_spiral.png")
LIGHT = LightSource(azdeg=315, altdeg=35)
FIBRE = "#3987e5"


def arrow3d(ax, p, v, color, lw=2.2):
    q = p + v
    ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], color=color, linewidth=lw, zorder=20)
    ax.scatter([q[0]], [q[1]], [q[2]], color=color, s=12, zorder=21, depthshade=False)


def grid_tris(X):
    """Triangles of a (3, n, m) grid surface."""
    A, B, C, D = X[:, :-1, :-1], X[:, 1:, :-1], X[:, 1:, 1:], X[:, :-1, 1:]
    t1 = np.stack([A, B, C], axis=-1).reshape(3, -1, 3)
    t2 = np.stack([A, C, D], axis=-1).reshape(3, -1, 3)
    return np.concatenate([t1, t2], axis=1).transpose(1, 2, 0)


def mobius_panel(ax):
    half = 0.45
    th = np.linspace(0, 2 * np.pi, 121)
    v = np.linspace(-half, half, 7)
    TH, V = np.meshgrid(th, v, indexing="ij")
    draw_meshes(ax, [(grid_tris(mobius_strip(TH, V)), BLUES[1], 1.0)])
    for t in np.linspace(0, 2 * np.pi, 25)[:-1]:
        seg = mobius_strip(np.array([t, t]), np.array([-half, half]))
        ax.plot(seg[0], seg[1], seg[2], color=INK, linewidth=0.6, zorder=10)
    C = base_frame(th, Circle(1.0))[0]
    ax.plot(C[0], C[1], C[2], color=MUTED, linewidth=1.0, ls="--", zorder=9)
    for t in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2):
        P = base_frame(t, Circle(1.0))[0]
        arrow3d(ax, P, 0.62 * fibre(t, np.pi / 2 - t / 2, Circle(1.0)), ACCENT)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_zlim(-0.65, 0.65)
    ax.set_box_aspect((2.8, 2.8, 1.3), zoom=1.3)
    ax.set_axis_off()
    ax.view_init(elev=34, azim=-62)


def spiral_panel(ax):
    th = np.linspace(0, 2 * np.pi * 0.95, 400)
    curves = [(Circle(0.45), MUTED, "--", r"circle  $\rho = $ const"),
              (ArchimedeanSpiral(0.2, 0.0771), BLUES[2], "-", r"Archimedean  $\rho = a + b\theta$"),
              (LogSpiral(0.2, 0.2), BLADE, "-", r"logarithmic  $\rho = a\,e^{b\theta}$")]
    for sp, col, ls, lab in curves:
        P = base_frame(th, sp)[0]
        ax.plot(P[0], P[1], color=col, ls=ls, linewidth=2.6 if ls == "-" else 1.4, label=lab)
    for t in np.radians([0, 60, 120, 180, 240, 300]):
        P = base_frame(t, LogSpiral(0.2, 0.2))[0]
        ax.annotate("", xy=(P[0], P[1]), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=1.2, shrinkA=0, shrinkB=0))
    P = base_frame(np.radians(240), LogSpiral(0.2, 0.2))[0]
    ax.text(0.5 * P[0] - 0.02, 0.5 * P[1], r"$\rho(\theta)$", color=ACCENT, fontsize=13, ha="right")
    ax.plot(0, 0, "o", color=INK, ms=3)
    ax.set_aspect("equal")
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.72, 0.62)
    ax.axis("off")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), fontsize=10, frameon=False, ncol=1)


def blade_panel(ax):
    b = NOTEBOOK
    draw_planes(ax, 0.78, b.height)
    th, wv, TH, W = grid(b, 160, 12)
    X = surface(TH, W, b)
    ax.plot_surface(X[0], X[1], X[2], color=BLUES[1], alpha=0.85, linewidth=0, shade=True,
                    lightsource=LIGHT)
    for i in range(0, len(th), 8):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], X[2, i, [0, -1]], color=INK, linewidth=0.6, zorder=10)
    ax.plot(X[0, :, 0], X[1, :, 0], X[2, :, 0], color=MUTED, linewidth=1.0, ls="--", zorder=9)
    for t in (0.0, b.wrap / 3, 2 * b.wrap / 3, b.wrap):
        P = base_frame(t, b.spiral)[0]
        beta = float(ruling_angle(t, b))
        arrow3d(ax, P, b.height / np.cos(beta) * fibre(t, beta, b.spiral), ACCENT)
    setup_3d(ax, 0.78, b.height, elev=28, azim=-60, zoom=1.2)


def main():
    fig = plt.figure(figsize=(16, 6.4))
    fig.patch.set_facecolor("white")

    mobius_panel(fig.add_axes([0.0, 0.22, 0.34, 0.66], projection="3d", computed_zorder=False))
    spiral_panel(fig.add_axes([0.37, 0.27, 0.26, 0.6]))
    blade_panel(fig.add_axes([0.66, 0.22, 0.34, 0.66], projection="3d", computed_zorder=False))

    heads = [(0.17, "(a) Circle + fibre turning at half speed"),
             (0.5, "(b) A spiral is a circle whose radius\ndepends on the angle"),
             (0.83, "(c) Spiral + slowly turning fibre,\nstretched between two planes")]
    for x, t in heads:
        fig.text(x, 0.93, t, ha="center", va="center", color=INK, fontsize=13)
    formulas = [
        (0.17, r"$X = C(\theta) + w\,F(\theta),\quad C(\theta) = (\cos\theta,\ \sin\theta,\ 0)$" "\n"
               r"$F = \cos\frac{\theta}{2}\,R + \sin\frac{\theta}{2}\,Z$", "Möbius strip"),
        (0.5, r"$P(\theta) = \rho(\theta)\,(\cos\theta,\ \sin\theta,\ 0)$", None),
        (0.83, r"$X = P(\theta) + \frac{w}{\cos\beta}\,F(\theta)$" "\n"
               r"$F = -\sin\beta\,N + \cos\beta\,Z,\quad \beta = \frac{\pi}{4} - \frac{\theta}{8}$",
         "pump blade"),
    ]
    for x, f, name in formulas:
        fig.text(x, 0.12 if name else 0.07, f, ha="center", va="center", color=INK, fontsize=14,
                 linespacing=1.6)
        if name:
            fig.text(x, 0.02, name, ha="center", color=ACCENT, fontsize=13)
    fig.savefig(OUT, dpi=130, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
