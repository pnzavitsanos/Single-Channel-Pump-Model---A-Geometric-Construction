"""
render_geometry.py
===================

media/geometry.png, three views of the blade's differential geometry:

  (a) the rulings: straight segments in the normal planes of the camber line,
      each running from the hub plane to the shroud plane, coloured by |K|;
  (b) horizontal slices: at height c the blade is the offset curve
      P - c tan(lambda) N of the camber line;
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
from impeller_geometry import NOTEBOOK, grid, surface, gaussian_curvature, level_curve, slice_blade_angle

OUT = os.path.join(MEDIA, "geometry.png")
RADIUS = 0.78


def main(d=NOTEBOOK):
    th, wv, TH, W = grid(d, 240, 40)
    X = surface(TH, W, d)
    K = gaussian_curvature(TH, W, d)
    norm = Normalize(vmin=0, vmax=np.abs(K).max())

    fig = plt.figure(figsize=(16, 5.4))
    fig.patch.set_facecolor("white")

    ax = fig.add_subplot(1, 3, 1, projection="3d", computed_zorder=False)
    draw_planes(ax, RADIUS, d.height)
    Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
    ax.plot_surface(X[0], X[1], X[2], facecolors=CURVATURE_CMAP(norm(np.abs(Kc))),
                    linewidth=0, antialiased=False, shade=False, alpha=0.95)
    for i in range(0, len(th), 12):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], X[2, i, [0, -1]], color=INK, linewidth=0.8)
    setup_3d(ax, RADIUS, d.height, elev=30, azim=-90, zoom=1.15)
    ax.set_title("(a) Rulings from hub plane to shroud plane", color=INK, fontsize=12)

    ax = fig.add_subplot(1, 3, 2)
    for c, col in zip(np.linspace(0, d.height, 5), BLUES[2:]):
        p = level_curve(c, th, d)
        label = "$c=0$ (hub: camber line)" if c == 0 else (f"$c={c:.3f}$" + (" (shroud)" if c == d.height else ""))
        ax.plot(p[0], p[1], color=INK if c == 0 else col, linewidth=1.6 if c == 0 else 2, label=label)
    for i in range(0, len(th), 12):
        ax.plot(X[0, i, [0, -1]], X[1, i, [0, -1]], color=MUTED, linewidth=0.7)
    b_in, b_out = np.degrees(slice_blade_angle(d.height, np.array([0.0, d.wrap]), d))
    ax.text(0.02, 0.03, f"blade angle at the shroud: inlet {b_in:.1f}°, outlet {b_out:.1f}°\n"
            f"(hub: {d.beta1_deg:.1f}° and {d.beta2_deg:.1f}°)",
            transform=ax.transAxes, fontsize=9, color=INK)
    ax.set_aspect("equal")
    ax.set_xlim(-0.62, 1.22)
    ax.set_ylim(-0.8, 0.52)
    ax.legend(loc="upper right", fontsize=9, frameon=False, title="slice height", title_fontsize=9)
    ax.set_title(r"(b) Slice at height $c$:  $P(\theta)-c\,\tan\lambda(\theta)\,N(\theta)$",
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
    ax.set_xlabel(r"$\theta$ (degrees along the blade)", color=INK)
    ax.set_ylabel(r"height $w$", color=INK)
    ax.set_title(r"(c) Gaussian curvature $K(\theta,w)$, negative everywhere", color=INK, fontsize=12)
    ax.tick_params(colors=MUTED, labelsize=9)

    fig.subplots_adjust(left=0.0, right=0.95, bottom=0.12, top=0.86, wspace=0.18)
    fig.savefig(OUT, dpi=140, facecolor="white")
    trim_png(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
