"""
make_control_function_animation.py
====================================

Renders ``media/control_function_sweep.gif``.

Unlike ``make_animation.py`` (which fixes the twist law and animates the
inlet/outlet *theta* domain), this animation fixes the theta domain at the
baseline sweep ``[0, 342 deg]`` (1.9*pi) and instead animates the **control
function itself**:

    beta(theta) = beta_max * theta / theta_max

by sweeping ``beta_max`` from 0 deg up to 180 deg and back down, smoothly and
periodically. At beta_max = 0 the channel is untwisted (a flat offset band);
at beta_max = 180 deg it carries a full Mobius-style half-twist from one end
to the other, applied over the partial spiral sweep instead of a full turn.

Two panels are drawn side by side every frame:

  * left:  the 3D channel surface for the current beta_max
  * right: the control function beta(theta) itself, i.e. the straight line
           from (0, 0) to (theta_max, beta_max) in a (theta, beta) plot with
           fixed axes -- so you watch the *law* rotate into place, which is
           the thing this animation is meant to show changing.

Usage:
    python scripts/make_control_function_animation.py
"""
import os
import sys
import shutil
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from channel_model import control_function_surface, control_function_values, BASELINE_OUTLET_DEG

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "media", "control_function_sweep.gif")

N_FRAMES = 30
FPS = 10
N_THETA, N_W = 60, 9

THETA_MAX_DEG = BASELINE_OUTLET_DEG  # 342 deg = 1.9*pi, held fixed throughout


def beta_max_schedule(frame_t: float) -> float:
    """frame_t in [0, 1); smooth 0 -> 180 -> 0 deg ping-pong via a raised cosine."""
    return 90.0 * (1.0 - np.cos(2 * np.pi * frame_t))


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_cf_anim_")
    fig = plt.figure(figsize=(11, 5.5))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.04, right=0.98, top=0.88, bottom=0.12, wspace=0.25)

    frame_paths = []
    for i in range(N_FRAMES):
        beta_max = beta_max_schedule(i / N_FRAMES)

        # --- left: 3D surface ---
        ax3d.clear()
        x, y, z, th, wv = control_function_surface(THETA_MAX_DEG, beta_max,
                                                     n_theta=N_THETA, n_w=N_W)
        theta_norm = (th - th.min()) / (th.max() - th.min() + 1e-9)
        facecolors = plt.cm.plasma(np.tile(theta_norm[:, None], (1, N_W)))
        ax3d.plot_surface(x, y, z, facecolors=facecolors, linewidth=0,
                           antialiased=False, shade=False)
        ax3d.set_box_aspect((1, 1, 0.4))
        ax3d.set_axis_off()
        ax3d.view_init(elev=32, azim=-58)
        ax3d.set_xlim(-0.6, 0.7)
        ax3d.set_ylim(-0.7, 0.5)
        ax3d.set_zlim(-0.15, 0.15)
        ax3d.set_title(f"channel surface\n$\\beta_{{max}}$ = {beta_max:5.1f}°", fontsize=11)

        # --- right: the control function itself ---
        ax2d.clear()
        th_deg, beta_deg = control_function_values(THETA_MAX_DEG, beta_max, n_theta=200)
        ax2d.plot(th_deg, beta_deg, color="#c23b32", linewidth=2.5)
        ax2d.fill_between(th_deg, 0, beta_deg, color="#e8823c", alpha=0.15)
        ax2d.axhline(90, color="gray", linestyle=":", linewidth=1)
        ax2d.set_xlim(0, THETA_MAX_DEG)
        ax2d.set_ylim(0, 180)
        ax2d.set_xlabel(r"$\theta$ (deg)  —  fixed inlet$\to$outlet sweep")
        ax2d.set_ylabel(r"control function $\beta(\theta)$ (deg)")
        ax2d.set_title(r"$\beta(\theta) = \beta_{max}\cdot\theta/\theta_{max}$", fontsize=11)
        ax2d.spines["top"].set_visible(False)
        ax2d.spines["right"].set_visible(False)

        p = os.path.join(tmpdir, f"frame_{i:03d}.png")
        fig.savefig(p, dpi=100, facecolor="white")
        frame_paths.append(p)

    images = [Image.open(p).convert("P", palette=Image.ADAPTIVE, colors=200) for p in frame_paths]
    images[0].save(
        OUT, save_all=True, append_images=images[1:],
        duration=int(1000 / FPS), loop=0, optimize=False,
    )
    shutil.rmtree(tmpdir, ignore_errors=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
