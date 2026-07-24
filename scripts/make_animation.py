"""
make_animation.py
==================

Renders ``media/inlet_outlet_sweep.gif``: an animation of the log-spiral
channel (see ``src/channel_model.py``) as its inlet and outlet angles are
perturbed by small, independent sinusoidal deviations around the baseline
values used in ``notebooks/05_log_spiral_v2_final.nb``.

This is meant to demonstrate that the channel is a *parametric* design: the
inlet angle (where the spiral sweep begins) and the outlet angle (where it
ends) can be dialed in independently without touching the underlying spiral
or twist law.

Usage:
    python scripts/make_animation.py
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
from channel_model import channel_surface, BASELINE_INLET_DEG, BASELINE_OUTLET_DEG

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "media", "inlet_outlet_sweep.gif")

N_FRAMES = 48
FPS = 12
N_THETA, N_W = 70, 10

# Small-deviation amplitudes (degrees). Independent, 90 deg out of phase so
# the viewer can see the inlet and outlet move on their own schedules rather
# than always in lockstep.
INLET_AMP_DEG = 12.0
OUTLET_AMP_DEG = 15.0


def inlet_outlet(frame_t: float):
    """frame_t in [0, 1); returns (inlet_deg, outlet_deg)."""
    inlet = BASELINE_INLET_DEG + INLET_AMP_DEG * np.sin(2 * np.pi * frame_t)
    outlet = BASELINE_OUTLET_DEG + OUTLET_AMP_DEG * np.sin(2 * np.pi * frame_t + np.pi / 2)
    return inlet, outlet


def main():
    tmpdir = tempfile.mkdtemp(prefix="pump_anim_")
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("white")

    frame_paths = []
    for i in range(N_FRAMES):
        ax.clear()
        t = i / N_FRAMES
        inlet, outlet = inlet_outlet(t)
        x, y, z, th, wv = channel_surface(inlet, outlet, n_theta=N_THETA, n_w=N_W)

        theta_norm = (th - th.min()) / (th.max() - th.min() + 1e-9)
        facecolors = plt.cm.plasma(np.tile(theta_norm[:, None], (1, N_W)))

        ax.plot_surface(x, y, z, facecolors=facecolors, linewidth=0,
                         antialiased=False, shade=False)
        ax.set_box_aspect((1, 1, 0.4))
        ax.set_axis_off()
        ax.view_init(elev=34, azim=-58)
        ax.set_xlim(-0.6, 0.7)
        ax.set_ylim(-0.6, 0.5)
        ax.set_zlim(-0.15, 0.15)
        ax.set_title(
            "Single-channel pump skeleton — inlet/outlet angle sensitivity\n"
            f"inlet = {inlet:6.1f}°    outlet = {outlet:6.1f}°   "
            f"(baseline {BASELINE_INLET_DEG:.0f}° / {BASELINE_OUTLET_DEG:.0f}°)",
            fontsize=11,
        )
        p = os.path.join(tmpdir, f"frame_{i:03d}.png")
        fig.savefig(p, dpi=110, facecolor="white")
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
