"""
make_fibre_animation.py
========================

Renders media/rotating_fibre.gif: a fibre carried twice around a circle while
it turns about the tangent at half the speed of the point,
F = cos(theta/2) R + sin(theta/2) Z. The fibre sweeps out the Moebius strip.

After one turn (theta = 360 deg) the fibre has rotated by 180 deg and points
the opposite way; only after two turns (720 deg) is it back where it started.

Left: the strip, the moving point and its fibre. Right: the plane across the
circle at that point (directions R and Z), with the fibre as a clock hand.

Usage:
    python scripts/make_fibre_animation.py [--mp4 path]
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, MUTED, BLUES, ACCENT, draw_meshes, save_gif, save_mp4
from blade_geometry import Circle, base_frame, fibre, mobius_strip
from render_fibre import grid_tris, arrow3d

OUT = os.path.join(MEDIA, "rotating_fibre.gif")
FPS = 12
N_FRAMES = 72
HOLD = 14
HALF = 0.45


def strip_tris():
    th = np.linspace(0, 2 * np.pi, 121)
    v = np.linspace(-HALF, HALF, 7)
    TH, V = np.meshgrid(th, v, indexing="ij")
    return grid_tris(mobius_strip(TH, V))


def state_text(deg):
    if abs(deg - 360) < 1e-6:
        return "one turn: the fibre has turned 180° and points the other way"
    if abs(deg - 720) < 1e-6:
        return "two turns: the fibre is back where it started"
    return ""


def draw_frame(fig, theta, tris):
    fig.clf()
    deg = np.degrees(theta)
    circle = Circle(1.0)

    ax = fig.add_axes([0.0, 0.02, 0.58, 0.86], projection="3d", computed_zorder=False)
    draw_meshes(ax, [(tris, BLUES[1], 1.0)])
    for t in np.linspace(0, 2 * np.pi, 25)[:-1]:
        rul = mobius_strip(np.array([t, t]), np.array([-HALF, HALF]))
        ax.plot(rul[0], rul[1], rul[2], color=INK, linewidth=0.4, alpha=0.6, zorder=8)
    C = base_frame(np.linspace(0, 2 * np.pi, 200), circle)[0]
    ax.plot(C[0], C[1], C[2], color=MUTED, linewidth=1.0, ls="--", zorder=9)
    P = base_frame(theta, circle)[0]
    R = np.array([np.cos(theta), np.sin(theta), 0.0])
    Z = np.array([0.0, 0.0, 1.0])
    F = fibre(theta, np.pi / 2 - theta / 2, circle)
    seg = mobius_strip(np.array([theta, theta]), np.array([-HALF, HALF]))
    ax.plot(seg[0], seg[1], seg[2], color=INK, linewidth=2.0, zorder=15)
    for vec, name in ((R, "R"), (Z, "Z")):
        q = P + 0.35 * vec
        ax.plot([P[0], q[0]], [P[1], q[1]], [P[2], q[2]], color=MUTED, linewidth=1.4, zorder=16)
        ax.text(q[0], q[1], q[2], f" {name}", color=MUTED, fontsize=11, zorder=17)
    arrow3d(ax, P, 0.62 * F, ACCENT, lw=3)
    ax.text(*(P + 0.7 * F), " F", color=ACCENT, fontsize=13, zorder=22)
    ax.scatter([P[0]], [P[1]], [P[2]], color=INK, s=18, zorder=23, depthshade=False)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_zlim(-0.65, 0.65)
    ax.set_box_aspect((2.8, 2.8, 1.3), zoom=1.3)
    ax.set_axis_off()
    ax.view_init(elev=34, azim=-62)

    ax2 = fig.add_axes([0.62, 0.2, 0.3, 0.6])
    t = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(np.cos(t), np.sin(t), color="#d0d7de", linewidth=1)
    ax2.annotate("", xy=(1.25, 0), xytext=(-1.25, 0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
    ax2.annotate("", xy=(0, 1.25), xytext=(0, -1.25), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
    ax2.text(1.3, -0.04, "R", color=MUTED, fontsize=12, va="center")
    ax2.text(0.05, 1.3, "Z", color=MUTED, fontsize=12)
    phi = theta / 2
    arc = np.linspace(0, phi, 100)
    r_arc = 0.35 + 0.08 * arc / np.pi
    ax2.plot(r_arc * np.cos(arc), r_arc * np.sin(arc), color=ACCENT, linewidth=1.5)
    ax2.annotate("", xy=(np.cos(phi), np.sin(phi)), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=3, mutation_scale=18))
    ax2.text(1.12 * np.cos(phi + 0.2), 1.12 * np.sin(phi + 0.2), "F", color=ACCENT, fontsize=13, ha="center", va="center")
    ax2.set_xlim(-1.45, 1.45)
    ax2.set_ylim(-1.45, 1.45)
    ax2.set_aspect("equal")
    ax2.axis("off")

    fig.text(0.77, 0.86, rf"point:  $\theta$ = {deg:5.0f}°", color=INK, fontsize=14, ha="center")
    fig.text(0.77, 0.80, rf"fibre:  $\theta/2$ = {deg / 2:5.0f}°", color=ACCENT, fontsize=14, ha="center")
    fig.text(0.77, 0.12, state_text(deg), color=INK, fontsize=11.5, ha="center")
    fig.text(0.5, 0.94, r"The rotating fibre:  $F = \cos\frac{\theta}{2}\,R + \sin\frac{\theta}{2}\,Z$",
             ha="center", color=INK, fontsize=16)
    fig.text(0.29, 0.03, "carried around the circle, the fibre sweeps out a Möbius strip", ha="center",
             color=MUTED, fontsize=10.5)


def main(mp4=None):
    tris = strip_tris()
    thetas = np.linspace(0, 4 * np.pi, N_FRAMES + 1)
    tmpdir = tempfile.mkdtemp(prefix="pump_fibre_")
    fig = plt.figure(figsize=(12, 5.6))
    fig.patch.set_facecolor("white")
    paths = []
    for i, th in enumerate(thetas):
        draw_frame(fig, th, tris)
        p = os.path.join(tmpdir, f"f{i:03d}.png")
        fig.savefig(p, dpi=80, facecolor="white")
        paths.append(p)
    save_gif(paths, OUT, FPS, holds={0: HOLD // 2, N_FRAMES // 2: HOLD, N_FRAMES: HOLD})
    if mp4:
        save_mp4(paths, mp4, FPS)
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--preview":
        fig = plt.figure(figsize=(12, 5.6))
        tris = strip_tris()
        for k in map(float, args[1:]):
            draw_frame(fig, np.radians(k), tris)
            fig.savefig(f"preview_f{int(k)}.png", dpi=80, facecolor="white")
    else:
        main(args[1] if len(args) > 1 and args[0] == "--mp4" else None)
