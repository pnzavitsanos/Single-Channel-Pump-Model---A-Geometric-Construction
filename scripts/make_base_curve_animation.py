"""
make_base_curve_animation.py
=============================

Renders media/base_curves.gif: the same construction on different spirals.

The inlet and outlet angles, the angular extent and the height stay fixed.
First the opening rate b of the logarithmic spiral rho = a e^{b theta} is
varied (tighter or looser winding), then the base curve is moved continuously
from the logarithmic spiral to an Archimedean spiral rho = a + b theta and
back. The blade still runs from plane to plane at every step, because that
property does not depend on the base curve.

Left: the base curve (solid) and the top edge of the blade (dashed), seen from
above. Right: the bottom plate and the blade.

Usage:
    python scripts/make_base_curve_animation.py [--mp4 path]
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, MUTED, BLADE, ACCENT, draw_pump, setup_3d, save_gif, save_mp4
from blade_geometry import (NOTEBOOK, LogSpiral, ArchimedeanSpiral, Blend, frame, level_curve,
                            spiral_angle, curvature)

OUT = os.path.join(MEDIA, "base_curves.gif")
FPS = 12
N_SWEEP = 40
HOLD = 8
LIM = 0.95
R_VIEW = 0.95
ARCH = ArchimedeanSpiral(0.2, 0.0771)


def schedule():
    frames = [("log", LogSpiral(0.2, 0.2 + 0.05 * np.sin(2 * np.pi * i / N_SWEEP))) for i in range(N_SWEEP)]
    frames += [("blend", Blend(LogSpiral(0.2, 0.2), ARCH, 0.5 * (1 - np.cos(2 * np.pi * i / N_SWEEP))))
               for i in range(N_SWEEP)]
    holds = {0: HOLD, N_SWEEP: HOLD, N_SWEEP + N_SWEEP // 2: HOLD}
    return frames, holds


def describe(kind, sp):
    if kind == "log":
        return ("Logarithmic spiral: how fast it opens", rf"$\rho = 0.2\,e^{{b\theta}}$,   b = {sp.b:.3f}")
    if sp.t < 0.02:
        name = r"logarithmic  $\rho = 0.2\,e^{0.2\theta}$"
    elif sp.t > 0.98:
        name = r"Archimedean  $\rho = 0.2 + 0.077\,\theta$"
    else:
        name = f"{100 * (1 - sp.t):.0f}% logarithmic + {100 * sp.t:.0f}% Archimedean"
    return ("From a logarithmic to an Archimedean spiral", name)


def draw_frame(fig, kind, sp):
    b = NOTEBOOK.with_(spiral=sp)
    fig.clf()
    title, name = describe(kind, sp)

    ax = fig.add_axes([0.05, 0.06, 0.4, 0.74])
    th = np.linspace(0, b.wrap, 300)
    P = frame(th, b)[0]
    top = level_curve(b.height, th, b)
    ax.plot(top[0], top[1], color=BLADE, lw=1.4, ls="--", label="top edge (z = H)")
    ax.plot(P[0], P[1], color=BLADE, lw=3.2, solid_capstyle="round", label="base curve (z = 0)")
    ax.plot(P[0, 0], P[1, 0], "o", color=INK, ms=6)
    ax.text(P[0, 0] + 0.04, P[1, 0], "inlet", color=INK, fontsize=10, va="center")
    ax.plot(P[0, -1], P[1, -1], "o", color=INK, ms=6)
    ax.text(P[0, -1] + 0.04, P[1, -1], "outlet", color=INK, fontsize=10, va="center")
    ax.plot(0, 0, marker="+", color=MUTED, ms=8)
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=9, frameon=False)

    ends = np.array([0.0, b.wrap])
    psi = np.degrees(spiral_angle(ends, b))
    kap = curvature(ends, b)
    fig.text(0.02, 0.86, name, color=ACCENT, fontsize=13.5)
    fig.text(0.02, 0.81, f"angle to the circle  inlet {psi[0]:4.1f}°   outlet {psi[1]:4.1f}°      "
             f"curvature  inlet {kap[0]:.2f}   outlet {kap[1]:.2f}", color=INK, fontsize=10.5)

    ax3 = fig.add_axes([0.47, 0.02, 0.53, 0.84], projection="3d", computed_zorder=False)
    draw_pump(ax3, b, top_sector=None, n_theta=120)
    setup_3d(ax3, R_VIEW, b.height, elev=34, azim=-90, zoom=1.35)

    fig.text(0.5, 0.95, title, ha="center", color=ACCENT, fontsize=16)
    fig.text(0.73, 0.04, f"fixed: inlet angle {b.inlet_angle_deg:.0f}°, outlet angle "
             f"{b.outlet_angle_deg:.1f}°, extent {b.wrap_deg:.0f}°, H = {b.height}",
             ha="center", color=MUTED, fontsize=10.5)


def main(mp4=None):
    frames, holds = schedule()
    tmpdir = tempfile.mkdtemp(prefix="pump_curves_")
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
            fig.savefig(f"preview_c{k}.png", dpi=80, facecolor="white")
    else:
        main(args[1] if len(args) > 1 and args[0] == "--mp4" else None)
