"""
make_control_function_animation.py
====================================

Renders media/control_function_sweep.gif: what the lean law does to the
curvature of the blade.

The lean runs linearly from beta_in = 45 deg at the inlet to beta_out at the
outlet, with the 1/cos(beta) stretch kept, so every frame still fits exactly
between the planes z = 0 and z = H. beta_out swings 45 -> -45 -> 45 deg:

  * beta_out = 45 deg: constant lean, beta' = 0, so K = 0 everywhere. The
    blade is developable and could be rolled from a flat sheet.
  * any other beta_out: beta' != 0 and K < 0 (saddle-shaped), with
    K = -beta'^2 / s^2 along the bottom edge.

Left: the blade coloured by |K| (fixed colour scale for all frames).
Right: the lean law beta(theta), with the forbidden band |beta| >= 90 deg in
which a ruling could never reach the top plane.

Usage:
    python scripts/make_control_function_animation.py
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
from channel_model import (DEFAULT_PARAMS as P, BASELINE_OUTLET_DEG,
                           control_function_surface, control_function_values)

OUT = os.path.join(MEDIA, "control_function_sweep.gif")
N_FRAMES = 36
FPS = 8
N_THETA, N_W = 90, 14
RADIUS = 0.8
THETA_MAX_DEG = BASELINE_OUTLET_DEG
BETA_IN_DEG = 45.0
K_NORM = PowerNorm(gamma=0.5, vmin=0.0, vmax=1.7)
LIGHT = LightSource(azdeg=315, altdeg=35)


def beta_out_schedule(t):
    """t in [0, 1): 45 -> -45 -> 45 deg, smooth."""
    return BETA_IN_DEG - 45.0 * (1.0 - np.cos(2 * np.pi * t))


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_cf_")
    fig = plt.figure(figsize=(11, 5))
    fig.patch.set_facecolor("white")

    paths = []
    for i in range(N_FRAMES):
        beta_out = beta_out_schedule(i / N_FRAMES)
        fig.clf()

        ax3d = fig.add_subplot(1, 2, 1, projection="3d", computed_zorder=False)
        draw_planes(ax3d, RADIUS, P.height)
        x, y, z, K, th, wv = control_function_surface(THETA_MAX_DEG, BETA_IN_DEG, beta_out,
                                                      n_theta=N_THETA, n_w=N_W)
        Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
        ax3d.plot_surface(x, y, z, facecolors=CURVATURE_CMAP(0.15 + 0.85 * K_NORM(np.abs(Kc))),
                          linewidth=0, antialiased=False, shade=True, lightsource=LIGHT)
        for edge in (np.s_[:, 0], np.s_[:, -1], np.s_[0, :], np.s_[-1, :]):
            ax3d.plot(x[edge], y[edge], z[edge], color=INK, linewidth=0.9)
        setup_3d(ax3d, RADIUS, P.height, elev=30, azim=-62, zoom=1.25)
        kmax = np.abs(K).max()
        state = "constant lean: developable, $K=0$" if kmax < 1e-12 else f"max $|K|$ = {kmax:.2f}"
        ax3d.set_title(f"blade coloured by $|K|$\n{state}", fontsize=11, color=INK)

        ax = fig.add_subplot(1, 2, 2)
        th_deg, beta_deg = control_function_values(THETA_MAX_DEG, BETA_IN_DEG, beta_out, n_theta=200)
        for lo, hi in ((90, 120), (-120, -90)):
            ax.axhspan(lo, hi, color="#d0d7de", alpha=0.6, linewidth=0)
        ax.text(THETA_MAX_DEG / 2, 100, "ruling cannot reach $z=H$", ha="center", va="center",
                fontsize=9, color=MUTED)
        ax.axhline(0, color="#d0d7de", linewidth=1)
        ax.plot(th_deg, beta_deg, color=BLUES[4], linewidth=2.5)
        slope = (beta_out - BETA_IN_DEG) / THETA_MAX_DEG
        ax.text(0.02, 0.06, rf"$\beta'$ = {slope:+.3f} (deg/deg)", transform=ax.transAxes,
                fontsize=10, color=INK)
        ax.set_xlim(0, THETA_MAX_DEG)
        ax.set_ylim(-120, 120)
        ax.set_xlabel(r"$\theta$ (degrees along the spiral)", color=INK)
        ax.set_ylabel(r"lean $\beta(\theta)$ (degrees)", color=INK)
        ax.set_title(rf"$\beta_{{in}}$ = {BETA_IN_DEG:.0f}°    $\beta_{{out}}$ = {beta_out:+.1f}°",
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
