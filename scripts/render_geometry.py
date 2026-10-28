"""
render_geometry.py
===================

media/geometry.png, three views of the blade's differential geometry:

  (a) the rulings: straight segments in the normal planes of the base curve,
      each running from the plane z = 0 to the plane z = H, coloured by |K|;
  (b) horizontal slices: at height c the blade is the offset curve
      P - c tan(beta) N of the base curve;
  (c) the Gaussian curvature K(theta, w) in closed form.

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
from blade_geometry import NOTEBOOK, grid, surface, gaussian_curvature, level_curve, spiral_angle, slice_spiral_angle

OUT = os.path.join(MEDIA, "geometry.png")
RADIUS = 0.78


def main(b=NOTEBOOK):
    th, wv, TH, W = grid(b, 240, 40)
    X = surface(TH, W, b)
    K = gaussian_curvature(TH, W, b)
    norm = Normalize(vmin=0, vmax=np.abs(K).max())

    fig = plt.figure(figsize=(16, 5.4))
    fig.patch.set_facecolor("white")

    ax = fig.add_subplot(1, 3, 1, projection="3d", computed_zorder=False)
    draw_planes(ax, RADIUS, b.height)
    Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
    ax.plot_surface(X[0], X[1], X[2], facecolors=CURVATURE_CMAP(norm(np.abs(Kc))),
                    linewidth=0, antialiased=False, shade=False, alpha=0.95)
    for i in range(0, len(th), 12):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], X[2, i, [0, -1]], color=INK, linewidth=0.8)
    setup_3d(ax, RADIUS, b.height, elev=30, azim=-90, zoom=1.15)
    ax.set_title("(a) Rulings from z = 0 to z = H, coloured by |K|", color=INK, fontsize=12)

    ax = fig.add_subplot(1, 3, 2)
    for c, col in zip(np.linspace(0, b.height, 5), BLUES[2:]):
        p = level_curve(c, th, b)
        label = "$c=0$ (base curve)" if c == 0 else (f"$c={c:.3f}$" + (" (top edge)" if c == b.height else ""))
        ax.plot(p[0], p[1], color=INK if c == 0 else col, linewidth=1.6 if c == 0 else 2, label=label)
    for i in range(0, len(th), 12):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], color=MUTED, linewidth=0.7)
    ends = np.array([0.0, b.wrap])
    p_in, p_out = np.degrees(spiral_angle(ends, b))
    t_in, t_out = np.degrees(slice_spiral_angle(b.height, ends, b))
    ax.text(0.02, 0.03, f"spiral angle at the top: inlet {t_in:.1f}°, outlet {t_out:.1f}°\n"
            f"(bottom: {p_in:.1f}° and {p_out:.1f}°)", transform=ax.transAxes, fontsize=9, color=INK)
    ax.set_aspect("equal")
    ax.set_xlim(-0.62, 1.22)
    ax.set_ylim(-0.8, 0.52)
    ax.legend(loc="upper right", fontsize=9, frameon=False, title="slice height", title_fontsize=9)
    ax.set_title(r"(b) Slice at height $c$:  $P(\theta)-c\,\tan\beta(\theta)\,N(\theta)$",
                 color=INK, fontsize=12)
    ax.tick_params(colors=MUTED, labelsize=9)
    for sp in ax.spines.values():
        sp.set_color("#d0d7de")

    ax = fig.add_subplot(1, 3, 3)
    im = ax.pcolormesh(np.rad2deg(TH), W, np.abs(K), cmap=CURVATURE_CMAP, norm=norm, shading="gouraud")
    cs = ax.contour(np.rad2deg(TH), W, K, levels=[-0.3, -0.2, -0.1, -0.05, -0.02],
                    colors=INK, linewidths=0.6)
    ax.clabel(cs, fmt="%.2f", fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label(r"$|K|$", color=INK)
    ax.set_xlabel(r"$\theta$ (degrees along the base curve)", color=INK)
    ax.set_ylabel(r"height $w$", color=INK)
    ax.set_title(r"(c) Gaussian curvature $K(\theta,w)$, negative everywhere", color=INK, fontsize=12)
    ax.tick_params(colors=MUTED, labelsize=9)

    fig.subplots_adjust(left=0.0, right=0.95, bottom=0.12, top=0.86, wspace=0.18)
    fig.savefig(OUT, dpi=140, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
