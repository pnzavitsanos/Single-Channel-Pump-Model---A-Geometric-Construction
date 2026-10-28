"""
make_angles_animation.py
=========================

Renders media/inlet_outlet_angles.gif: the inlet angle and the outlet angle of
the rulings, and what they do to the curvature of the blade.

The base curve (the notebook's logarithmic spiral), its angular extent and the
height stay fixed. First the inlet angle is swept, then the outlet angle. When
the two are equal the ruling angle is constant, K = 0 and the blade could be
rolled from flat sheet.

Left: the blade between the two planes, coloured by |K|, with the inlet and
outlet rulings and their angles from the vertical.
Right: the ruling angle along the blade, from inlet to outlet.

Usage:
    python scripts/make_angles_animation.py [--mp4 path]
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, PowerNorm

from plot_utils import (MEDIA, INK, MUTED, BLUES, ACCENT, CURVATURE_CMAP, draw_planes, setup_3d,
                        save_gif, save_mp4)
from blade_geometry import NOTEBOOK, grid, surface, gaussian_curvature, ruling_angle, frame

OUT = os.path.join(MEDIA, "inlet_outlet_angles.gif")
FPS = 12
N_SWEEP = 40
HOLD = 8
RADIUS = 1.0
K_NORM = PowerNorm(gamma=0.4, vmin=0.0, vmax=2.0, clip=True)
LIGHT = LightSource(azdeg=315, altdeg=35)


def ease(t):
    return 0.5 * (1 - np.cos(2 * np.pi * t))


def schedule():
    b_in, b_out = NOTEBOOK.inlet_angle_deg, NOTEBOOK.outlet_angle_deg
    frames = [("in", b_in - 70 * ease(i / N_SWEEP), b_out) for i in range(N_SWEEP)]
    frames += [("out", b_in, b_out + (b_in - b_out) * ease(i / N_SWEEP)) for i in range(N_SWEEP)]
    holds = {0: HOLD, N_SWEEP: HOLD, N_SWEEP + N_SWEEP // 2: HOLD}
    return frames, holds


def ruling_with_angle(ax, b, theta, label, color):
    """Draw one ruling, the vertical through its foot, and the angle between them."""
    P, T, N, _, _ = frame(theta, b)
    beta = float(ruling_angle(theta, b))
    top = surface(theta, b.height, b)
    ax.plot([P[0], top[0]], [P[1], top[1]], [P[2], top[2]], color=color, linewidth=3, zorder=10)
    ax.plot([P[0], P[0]], [P[1], P[1]], [0, b.height], color=MUTED, linewidth=1.2, ls="--", zorder=10)
    phi = np.linspace(0, beta, 20)
    r = 0.6 * b.height
    arc = P[:, None] + r * (np.cos(phi) * np.array([0, 0, 1.0])[:, None] - np.sin(phi) * N[:, None])
    ax.plot(arc[0], arc[1], arc[2], color=color, linewidth=1.8, zorder=10)
    tip = P + 1.15 * b.height * (np.cos(beta / 2) * np.array([0, 0, 1.0]) - np.sin(beta / 2) * N)
    ax.text(tip[0], tip[1], tip[2] + 0.02, label, color=color, fontsize=14, zorder=11)


def draw_frame(fig, which, a_in, a_out):
    b = NOTEBOOK.with_(inlet_angle_deg=a_in, outlet_angle_deg=a_out)
    fig.clf()

    ax = fig.add_axes([-0.02, 0.0, 0.64, 0.92], projection="3d", computed_zorder=False)
    draw_planes(ax, RADIUS, b.height)
    th, wv, TH, W = grid(b, 110, 14)
    X = surface(TH, W, b)
    K = gaussian_curvature(TH, W, b)
    Kc = 0.25 * (K[:-1, :-1] + K[1:, :-1] + K[:-1, 1:] + K[1:, 1:])
    ax.plot_surface(X[0], X[1], X[2], facecolors=CURVATURE_CMAP(0.12 + 0.88 * K_NORM(np.abs(Kc))),
                    linewidth=0, antialiased=False, shade=True, lightsource=LIGHT)
    for edge in (np.s_[:, 0], np.s_[:, -1]):
        ax.plot(X[0][edge], X[1][edge], X[2][edge], color=INK, linewidth=0.8)
    ruling_with_angle(ax, b, 0.0, r"$\beta_{in}$", ACCENT if which == "in" else INK)
    ruling_with_angle(ax, b, b.wrap, r"$\beta_{out}$", ACCENT if which == "out" else INK)
    setup_3d(ax, RADIUS, b.height, elev=24, azim=-90, zoom=1.32)

    ax2 = fig.add_axes([0.64, 0.14, 0.33, 0.5])
    th_deg = np.degrees(np.linspace(0, b.wrap, 100))
    beta_deg = np.degrees(ruling_angle(np.radians(th_deg), b))
    for lo, hi in ((90, 110), (-110, -90)):
        ax2.axhspan(lo, hi, color="#d0d7de", alpha=0.6, linewidth=0)
    ax2.axhline(0, color="#d0d7de", linewidth=1)
    ax2.plot(th_deg, beta_deg, color=BLUES[4], linewidth=2.5)
    for x, y, txt, col, ha in ((0, a_in, "  inlet", ACCENT if which == "in" else INK, "left"),
                               (b.wrap_deg, a_out, "outlet  ", ACCENT if which == "out" else INK, "right")):
        ax2.plot(x, y, "o", color=col, ms=9, clip_on=False, zorder=5)
        ax2.text(x, y + 10, txt, color=col, fontsize=11, ha=ha)
    ax2.set_xlim(0, b.wrap_deg)
    ax2.set_ylim(-110, 110)
    ax2.set_xlabel(r"$\theta$ along the base curve (degrees)", color=INK)
    ax2.set_ylabel(r"ruling angle $\beta$ (degrees)", color=INK)
    ax2.text(b.wrap_deg / 2, 100, "a ruling at 90° never reaches the top plane", ha="center",
             va="center", fontsize=8.5, color=MUTED)
    ax2.tick_params(colors=MUTED, labelsize=9)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)

    kmax = np.abs(K).max()
    fig.text(0.64, 0.86, rf"inlet angle    $\beta_{{in}}$ = {a_in:5.1f}°",
             color=ACCENT if which == "in" else INK, fontsize=14)
    fig.text(0.64, 0.80, rf"outlet angle  $\beta_{{out}}$ = {a_out:5.1f}°",
             color=ACCENT if which == "out" else INK, fontsize=14)
    state = ("constant angle: K = 0, can be rolled from flat sheet" if kmax < 1e-12
             else f"varying angle: saddle-shaped, max |K| = {kmax:.2f}")
    fig.text(0.64, 0.74, state, color=INK, fontsize=11.5)

    title = r"Varying the inlet angle $\beta_{in}$" if which == "in" else r"Varying the outlet angle $\beta_{out}$"
    fig.text(0.5, 0.95, title, ha="center", color=ACCENT, fontsize=16)
    fig.text(0.29, 0.03, "blade coloured by Gaussian curvature |K| (light = flat)", ha="center",
             color=MUTED, fontsize=10.5)


def main(mp4=None):
    frames, holds = schedule()
    tmpdir = tempfile.mkdtemp(prefix="pump_angles_")
    fig = plt.figure(figsize=(13, 5.8))
    fig.patch.set_facecolor("white")
    paths = []
    for i, f in enumerate(frames):
        draw_frame(fig, *f)
        p = os.path.join(tmpdir, f"f{i:03d}.png")
        fig.savefig(p, dpi=80, facecolor="white")
        paths.append(p)
    save_gif(paths, OUT, FPS, holds=holds)
    if mp4:
        save_mp4(paths, mp4, FPS)
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--preview":
        fig = plt.figure(figsize=(13, 5.8))
        for k in map(int, args[1:]):
            draw_frame(fig, *schedule()[0][k])
            fig.savefig(f"preview_{k}.png", dpi=80, facecolor="white")
    else:
        main(args[1] if len(args) > 1 and args[0] == "--mp4" else None)
