"""Shared styling for the figures and animations."""
import os
import sys

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

MEDIA = os.path.join(os.path.dirname(__file__), "..", "media")

INK = "#1f2328"
MUTED = "#6e7781"
PLANE = "#8c959f"

# One-hue blue ramp: light = flat (K = 0), dark = strongly curved.
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CURVATURE_CMAP = LinearSegmentedColormap.from_list("curvature", BLUES)


def draw_planes(ax, radius, height, n=96):
    """The two parallel planes z = 0 and z = H, drawn as translucent discs."""
    t = np.linspace(0, 2 * np.pi, n)
    rr = np.linspace(0, radius, 2)
    T, R = np.meshgrid(t, rr)
    for z in (0.0, height):
        ax.plot_surface(R * np.cos(T), R * np.sin(T), np.full_like(T, z),
                        color=PLANE, alpha=0.08, linewidth=0, shade=False)
        ax.plot(radius * np.cos(t), radius * np.sin(t), np.full_like(t, z),
                color=PLANE, linewidth=0.8)


def setup_3d(ax, radius, height, elev=24, azim=-62, zoom=1.0):
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(0, height)
    ax.set_box_aspect((2 * radius, 2 * radius, height), zoom=zoom)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)


def save_gif(frame_paths, out, fps, hold_last=0, hold_first=0):
    """Write an animated GIF with one shared palette (consistent colours on every frame)."""
    from PIL import Image
    frames = [Image.open(p).convert("RGB") for p in frame_paths]
    sample = Image.new("RGB", (frames[0].width, frames[0].height * 4))
    for k, idx in enumerate(np.linspace(0, len(frames) - 1, 4).astype(int)):
        sample.paste(frames[idx], (0, frames[0].height * k))
    palette = sample.quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    images = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    durations = [int(1000 / fps)] * len(images)
    durations[0] += hold_first * int(1000 / fps)
    durations[-1] += hold_last * int(1000 / fps)
    images[0].save(out, save_all=True, append_images=images[1:], duration=durations,
                   loop=0, optimize=False, disposal=1)


def trim_png(path, pad=12):
    """Crop the white margin that 3D axes leave around a saved PNG."""
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if box:
        l, t, r, b = box
        im.crop((max(l - pad, 0), max(t - pad, 0), min(r + pad, im.width), min(b + pad, im.height))).save(path)
