"""
render_geometry.py
===================

Three views of the mathematics behind the blade (media/geometry.png):

  (a) the rulings: straight segments in the normal planes of the spiral, each
      running from z = 0 to z = H, coloured by Gaussian curvature;
  (b) horizontal slices: at height c the blade is the offset curve
      r - c tan(beta) N of the base spiral;
  (c) the Gaussian curvature K(theta, w) in closed form: negative everywhere,
      strongest at the inlet foot, and zero only if the lean is constant.

Usage:
    python scripts/render_geometry.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from plot_utils import MEDIA, INK, MUTED, BLUES, CURVATURE_CMAP, draw_planes, setup_3d, trim_png
from channel_model import (DEFAULT_PARAMS as P, BASELINE_INLET_DEG, BASELINE_OUTLET_DEG,
                           frame, surface, gaussian_curvature, level_curve, notebook_lean)

OUT = os.path.join(MEDIA, "geometry.png")
RADIUS = 0.78


def main():
    th = np.linspace(np.deg2rad(BASELINE_INLET_DEG), np.deg2rad(BASELINE_OUTLET_DEG), 240)
    wv = np.linspace(0, P.height, 40)
    TH, W = np.meshgrid(th, wv, indexing="ij")
    X = surface(TH, W)
    K = gaussian_curvature(TH, W)
    norm = Normalize(vmin=0, vmax=np.abs(K).max())

    fig = plt.figure(figsize=(16, 5.4))
    fig.patch.set_facecolor("white")

    # (a) rulings, coloured by |K|
    ax = fig.add_subplot(1, 3, 1, projection="3d", computed_zorder=False)
    draw_planes(ax, RADIUS, P.height)
    Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
    ax.plot_surface(X[0], X[1], X[2], facecolors=CURVATURE_CMAP(norm(np.abs(Kc))),
                    linewidth=0, antialiased=False, shade=False, alpha=0.9)
    for i in range(0, len(th), 12):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], X[2, i, [0, -1]], color=INK, linewidth=0.8)
    setup_3d(ax, RADIUS, P.height, elev=28, azim=-62, zoom=1.15)
    ax.set_title("(a) Rulings run from plane $z=0$ to plane $z=H$", color=INK, fontsize=12)

    # (b) horizontal slices are offset curves of the spiral
    ax = fig.add_subplot(1, 3, 2)
    heights = np.linspace(0, P.height, 5)
    for c, col in zip(heights, BLUES[2::1]):
        p = level_curve(c, th)
        label = "$c=0$ (base spiral)" if c == 0 else f"$c={c:.3f}$"
        ax.plot(p[0], p[1], color=col if c > 0 else INK, linewidth=2 if c > 0 else 1.6, label=label)
    for i in range(0, len(th), 12):
        seg = X[:, i, [0, -1]]
        ax.plot(seg[0], seg[1], color=MUTED, linewidth=0.7)
    ax.set_aspect("equal")
    ax.set_title(r"(b) Slice at height $c$:  $r(\theta)-c\,\tan\beta(\theta)\,N(\theta)$",
                 color=INK, fontsize=12)
    ax.set_xlim(-0.62, 1.02)
    ax.legend(loc="upper right", fontsize=9, frameon=False, title="slice height", title_fontsize=9)
    ax.tick_params(colors=MUTED, labelsize=9)
    for sp in ax.spines.values():
        sp.set_color("#d0d7de")

    # (c) Gaussian curvature map
    ax = fig.add_subplot(1, 3, 3)
    im = ax.pcolormesh(np.rad2deg(TH), W, np.abs(K), cmap=CURVATURE_CMAP, norm=norm, shading="gouraud")
    cs = ax.contour(np.rad2deg(TH), W, K, levels=[-0.3, -0.2, -0.1, -0.05, -0.02],
                    colors=INK, linewidths=0.6)
    ax.clabel(cs, fmt="%.2f", fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label(r"$|K|$", color=INK)
    ax.set_xlabel(r"$\theta$ (degrees along the spiral)", color=INK)
    ax.set_ylabel(r"height $w$", color=INK)
    ax.set_title(r"(c) Gaussian curvature $K(\theta,w)$, negative everywhere", color=INK, fontsize=12)
    ax.tick_params(colors=MUTED, labelsize=9)

    fig.subplots_adjust(left=0.0, right=0.95, bottom=0.12, top=0.86, wspace=0.18)
    fig.savefig(OUT, dpi=140, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
