"""
make_angle_animation.py
========================

Renders media/blade_angles.gif: what the inlet and outlet blade angles do to
the impeller.

Left: plan view of the camber line on the hub, with the blade angle at the
inlet (beta1) and at the outlet (beta2), each measured from the
circumferential direction. Right: the impeller (hub and blade).

The radii r1, r2, the height H and the lean law stay fixed. First the outlet
angle beta2 is swept, then the inlet angle beta1. The wrap angle is not a free
parameter: it follows from r1, r2, beta1 and beta2, so the blade gets
shorter as the angles open up.

Usage:
    python scripts/make_angle_animation.py
"""
import os
import shutil
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, MUTED, BLADE, ACCENT, draw_impeller, setup_3d, save_gif
from impeller_geometry import NOTEBOOK, frame, blade_angle

OUT = os.path.join(MEDIA, "blade_angles.gif")
FPS = 12
N_SWEEP = 36          # frames per sweep (there and back)
HOLD = 6              # frames held before each sweep
LIM = 0.78            # plan-view half width
R_VIEW = 0.76         # 3D half width


def ease(t):
    """0 -> 1 -> 0 with zero speed at the ends."""
    return 0.5 * (1 - np.cos(2 * np.pi * t))


def schedule():
    b0 = NOTEBOOK.beta1_deg
    frames = [("beta2", b0, b0)] * HOLD
    frames += [("beta2", b0, b0 + (30 - b0) * ease(i / N_SWEEP)) for i in range(N_SWEEP)]
    frames += [("beta1", b0, b0)] * HOLD
    frames += [("beta1", b0 + (25 - b0) * ease(i / N_SWEEP), b0) for i in range(N_SWEEP)]
    return frames


def angle_marker(ax, theta, beta, label, color, d):
    """Circumferential direction, blade tangent, and the arc between them."""
    P = frame(theta, d)[0]
    a_c = theta + np.pi / 2            # circumferential direction
    a_t = a_c - beta                   # blade tangent
    L, rr = 0.3, 0.2
    for a, col, lw, ls in ((a_c, MUTED, 1.5, "--"), (a_t, color, 2.4, "-")):
        ax.plot([P[0], P[0] + L * np.cos(a)], [P[1], P[1] + L * np.sin(a)], color=col, lw=lw, ls=ls)
    arc = np.linspace(a_t, a_c, 20)
    ax.plot(P[0] + rr * np.cos(arc), P[1] + rr * np.sin(arc), color=color, lw=1.8)
    am = 0.5 * (a_t + a_c)
    side = am - np.pi / 2              # outside of the angle
    ax.text(P[0] + 0.25 * np.cos(am) + 0.1 * np.cos(side), P[1] + 0.25 * np.sin(am) + 0.1 * np.sin(side),
            label, color=color, fontsize=15, ha="center", va="center")


def draw_frame(fig, which, b1, b2):
    d = NOTEBOOK.with_(beta1_deg=b1, beta2_deg=b2)
    fig.clf()

    ax = fig.add_axes([0.08, 0.08, 0.42, 0.8])
    t = np.linspace(0, 2 * np.pi, 200)
    for r, name in ((d.r1, "$r_1$"), (d.r2, "$r_2$")):
        ax.plot(r * np.cos(t), r * np.sin(t), color="#d0d7de", lw=1, ls=":")
        ax.text((r + 0.05) * np.cos(2.36), (r + 0.05) * np.sin(2.36), name, color=MUTED, fontsize=11, ha="center", va="center")
    th = np.linspace(0, d.wrap, 300)
    P = frame(th, d)[0]
    ax.plot(P[0], P[1], color=BLADE, lw=3.2, solid_capstyle="round")
    ax.plot(0, 0, marker="+", color=MUTED, ms=8)
    c1 = ACCENT if which == "beta1" else INK
    c2 = ACCENT if which == "beta2" else INK
    angle_marker(ax, 0.0, float(blade_angle(0.0, d)), r"$\beta_1$", c1, d)
    angle_marker(ax, d.wrap, float(blade_angle(d.wrap, d)), r"$\beta_2$", c2, d)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Camber line on the hub (plan view)", color=INK, fontsize=13)

    fig.text(0.015, 0.84, rf"inlet angle    $\beta_1$ = {b1:4.1f}°", color=c1, fontsize=14)
    fig.text(0.015, 0.78, rf"outlet angle  $\beta_2$ = {b2:4.1f}°", color=c2, fontsize=14)
    fig.text(0.015, 0.72, rf"wrap angle    $\Phi$ = {d.wrap_deg:5.1f}°", color=INK, fontsize=14)

    ax3 = fig.add_axes([0.5, 0.06, 0.5, 0.84], projection="3d", computed_zorder=False)
    draw_impeller(ax3, d, shroud_sector=None)
    setup_3d(ax3, R_VIEW, d.height, elev=38, azim=-90, zoom=1.3)
    ax3.set_title("Impeller (hub and blade)", color=INK, fontsize=13)

    varying = "outlet angle $\\beta_2$" if which == "beta2" else "inlet angle $\\beta_1$"
    fig.text(0.5, 0.955, f"Varying the {varying}", ha="center", color=ACCENT, fontsize=16)
    fig.text(0.5, 0.025, f"fixed: $r_1$ = {d.r1:.2f},  $r_2$ = {d.r2:.2f},  $H$ = {d.height:.2f},  "
             f"lean {d.lean1_deg:.0f}° → {d.lean2_deg:.1f}°.   The wrap angle follows from "
             r"$r_1, r_2, \beta_1, \beta_2$.", ha="center", color=MUTED, fontsize=11.5)


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_angles_")
    fig = plt.figure(figsize=(13, 5.8))
    fig.patch.set_facecolor("white")
    paths = []
    for i, (which, b1, b2) in enumerate(schedule()):
        draw_frame(fig, which, b1, b2)
        p = os.path.join(tmpdir, f"f{i:03d}.png")
        fig.savefig(p, dpi=80, facecolor="white")
        paths.append(p)
    save_gif(paths, OUT, FPS)
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:   # preview single frames: python make_angle_animation.py 20 50
        fig = plt.figure(figsize=(13, 5.8))
        for k in map(int, sys.argv[1:]):
            draw_frame(fig, *schedule()[k])
            fig.savefig(f"preview_{k}.png", dpi=80, facecolor="white")
    else:
        main()
