"""
make_ruling_animation.py
=========================

Renders media/ruling_sweep.gif: the blade traced out by a single straight
ruling as its foot runs along the logarithmic spiral.

At each angle theta the ruling lies in the normal plane span{N, Z} of the
spiral, leans from the vertical by beta(theta) = pi/4 - theta/8 (a rotation
about the tangent T, as in the Moebius strip), and is stretched by
1/cos(beta) so that it always ends on the top plane z = H.

Usage:
    python scripts/make_ruling_animation.py
"""
import os
import shutil
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_utils import MEDIA, INK, MUTED, draw_planes, setup_3d, save_gif
from channel_model import (DEFAULT_PARAMS as P, BASELINE_INLET_DEG, BASELINE_OUTLET_DEG,
                           frame, surface, notebook_lean)

OUT = os.path.join(MEDIA, "ruling_sweep.gif")
RADIUS = 0.78
N_FRAMES = 44
FPS = 10


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_ruling_")
    fig = plt.figure(figsize=(6.4, 4.6))
    fig.patch.set_facecolor("white")
    th_end = np.linspace(np.deg2rad(BASELINE_INLET_DEG), np.deg2rad(BASELINE_OUTLET_DEG),
                         N_FRAMES + 1)[1:]
    wv = np.linspace(0, P.height, 8)

    paths = []
    for i, th_k in enumerate(th_end):
        fig.clf()
        ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
        draw_planes(ax, RADIUS, P.height)

        th = np.linspace(np.deg2rad(BASELINE_INLET_DEG), th_k, max(int(80 * th_k / th_end[-1]), 2))
        TH, W = np.meshgrid(th, wv, indexing="ij")
        X = surface(TH, W)
        ax.plot_surface(X[0], X[1], X[2], color="#6da7ec", alpha=0.85, linewidth=0, shade=True)

        r, T, N, _ = frame(th_k)
        beta, _ = notebook_lean(th_k)
        top = surface(th_k, P.height)
        ax.plot([r[0], top[0]], [r[1], top[1]], [r[2], top[2]], color=INK, linewidth=2.6)
        for vec, name in ((T, "T"), (N, "N"), (np.array([0.0, 0.0, 1.0]), "Z")):
            tip = r + 0.13 * vec
            ax.plot([r[0], tip[0]], [r[1], tip[1]], [r[2], tip[2]], color=MUTED, linewidth=1.4)
            ax.text(tip[0], tip[1], tip[2], f" {name}", color=MUTED, fontsize=9)

        setup_3d(ax, RADIUS, P.height, elev=30, azim=-62, zoom=1.45)
        ax.set_title("A straight ruling swept along the spiral\n"
                     rf"$\theta$ = {np.rad2deg(th_k):5.1f}°    lean $\beta$ = {np.rad2deg(beta):4.1f}°"
                     r"    length $H/\cos\beta$", fontsize=11, color=INK)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=0.88)
        p = os.path.join(tmpdir, f"f{i:03d}.png")
        fig.savefig(p, dpi=90, facecolor="white")
        paths.append(p)

    save_gif(paths, OUT, FPS, hold_last=12)  # hold on the finished blade
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
