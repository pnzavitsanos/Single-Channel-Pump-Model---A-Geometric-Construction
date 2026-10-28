"""
make_lean_animation.py
=======================

Renders media/lean_sweep.gif: what the lean of the rulings does to the
curvature of the blade.

The camber line (r1, r2, beta1, beta2) and the height H stay fixed. The inlet
lean is 45 deg and the outlet lean swings 45 -> -45 -> 45 deg:

  * outlet lean = inlet lean: constant lean, K = 0 everywhere, the blade is
    developable and can be rolled from flat sheet;
  * otherwise: the lean changes along the blade and K < 0 (saddle-shaped),
    with K = -(dlambda/dsigma)^2 along the hub.

Left: the blade coloured by |K| (same colour scale in every frame).
Right: the lean law, with the band |lambda| >= 90 deg in which a ruling could
never reach the shroud.

Usage:
    python scripts/make_lean_animation.py
"""
import os
import shutil
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, PowerNorm

from plot_utils import MEDIA, INK, MUTED, BLUES, CURVATURE_CMAP, draw_planes, setup_3d, save_gif
from impeller_geometry import NOTEBOOK, grid, surface, gaussian_curvature, lean

OUT = os.path.join(MEDIA, "lean_sweep.gif")
N_FRAMES = 36
FPS = 8
RADIUS = 0.8
LEAN_IN = 45.0
K_NORM = PowerNorm(gamma=0.5, vmin=0.0, vmax=1.7)
LIGHT = LightSource(azdeg=315, altdeg=35)


def lean_out_schedule(t):
    return LEAN_IN - 45.0 * (1.0 - np.cos(2 * np.pi * t))


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_lean_")
    fig = plt.figure(figsize=(11, 5))
    fig.patch.set_facecolor("white")

    paths = []
    for i in range(N_FRAMES):
        lean_out = lean_out_schedule(i / N_FRAMES)
        d = NOTEBOOK.with_(lean1_deg=LEAN_IN, lean2_deg=lean_out)
        fig.clf()

        ax3d = fig.add_subplot(1, 2, 1, projection="3d", computed_zorder=False)
        draw_planes(ax3d, RADIUS, d.height)
        th, wv, TH, W = grid(d, 90, 14)
        X = surface(TH, W, d)
        K = gaussian_curvature(TH, W, d)
        Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
        ax3d.plot_surface(X[0], X[1], X[2], facecolors=CURVATURE_CMAP(0.15 + 0.85 * K_NORM(np.abs(Kc))),
                          linewidth=0, antialiased=False, shade=True, lightsource=LIGHT)
        for edge in (np.s_[:, 0], np.s_[:, -1], np.s_[0, :], np.s_[-1, :]):
            ax3d.plot(X[0][edge], X[1][edge], X[2][edge], color=INK, linewidth=0.9)
        setup_3d(ax3d, RADIUS, d.height, elev=30, azim=-90, zoom=1.25)
        kmax = np.abs(K).max()
        state = "constant lean: developable, $K=0$" if kmax < 1e-12 else f"max $|K|$ = {kmax:.2f}"
        ax3d.set_title(f"blade coloured by $|K|$\n{state}", fontsize=11, color=INK)

        ax = fig.add_subplot(1, 2, 2)
        th_deg = np.degrees(np.linspace(0, d.wrap, 200))
        lam_deg = np.degrees(lean(np.radians(th_deg), d))
        for lo, hi in ((90, 120), (-120, -90)):
            ax.axhspan(lo, hi, color="#d0d7de", alpha=0.6, linewidth=0)
        ax.text(d.wrap_deg / 2, 100, "ruling cannot reach the shroud", ha="center", va="center",
                fontsize=9, color=MUTED)
        ax.axhline(0, color="#d0d7de", linewidth=1)
        ax.plot(th_deg, lam_deg, color=BLUES[4], linewidth=2.5)
        ax.text(0.02, 0.06, rf"$d\lambda/d\theta$ = {d.lean_rate:+.3f}", transform=ax.transAxes,
                fontsize=10, color=INK)
        ax.set_xlim(0, d.wrap_deg)
        ax.set_ylim(-120, 120)
        ax.set_xlabel(r"$\theta$ (degrees along the blade)", color=INK)
        ax.set_ylabel(r"lean $\lambda(\theta)$ (degrees)", color=INK)
        ax.set_title(rf"inlet lean $\lambda_1$ = {LEAN_IN:.0f}°    outlet lean $\lambda_2$ = {lean_out:+.1f}°",
                     fontsize=11, color=INK)
        ax.tick_params(colors=MUTED, labelsize=9)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

        fig.subplots_adjust(left=0.0, right=0.97, top=0.86, bottom=0.12, wspace=0.12)
        p = os.path.join(tmpdir, f"f{i:03d}.png")
        fig.savefig(p, dpi=85, facecolor="white")
        paths.append(p)

    save_gif(paths, OUT, FPS, hold_first=8)  # pause on the developable case
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
