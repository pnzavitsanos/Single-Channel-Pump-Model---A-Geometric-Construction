"""Shared styling and drawing helpers for the figures and animations."""
import os
import sys

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

MEDIA = os.path.join(os.path.dirname(__file__), "..", "media")

INK = "#1f2328"
MUTED = "#6e7781"
PLANE = "#8c959f"
PLATE_BOTTOM = "#8c959f"
PLATE_TOP = "#c9d1d9"
BLADE = "#3987e5"
ACCENT = "#c2410c"

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


def setup_3d(ax, radius, height, elev=30, azim=-90, zoom=1.0):
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(0, height)
    ax.set_box_aspect((2 * radius, 2 * radius, height), zoom=zoom)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)


def save_gif(frame_paths, out, fps, hold_last=0, hold_first=0, holds=None):
    """Animated GIF with one shared palette. holds: {frame index: extra frames to hold}."""
    from PIL import Image
    frames = [Image.open(p).convert("RGB") for p in frame_paths]
    sample = Image.new("RGB", (frames[0].width, frames[0].height * 4))
    for k, idx in enumerate(np.linspace(0, len(frames) - 1, 4).astype(int)):
        sample.paste(frames[idx], (0, frames[0].height * k))
    palette = sample.quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    images = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    step = int(1000 / fps)
    durations = [step] * len(images)
    durations[0] += hold_first * step
    durations[-1] += hold_last * step
    for i, extra in (holds or {}).items():
        durations[i] += extra * step
    images[0].save(out, save_all=True, append_images=images[1:], duration=durations,
                   loop=0, optimize=False, disposal=1)


def save_mp4(frame_paths, out, fps):
    """Same frames as an MP4 (needs ffmpeg on the PATH)."""
    import shutil
    import subprocess
    if not shutil.which("ffmpeg"):
        return
    folder = os.path.dirname(frame_paths[0])
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", os.path.join(folder, "f%03d.png"), "-pix_fmt", "yuv420p",
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-c:v", "libx264", "-crf", "18",
                    "-movflags", "+faststart", out], check=True)


def trim_png(path, pad=12):
    """Crop the white margin that 3D axes leave around a saved PNG."""
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if box:
        l, t, r, b = box
        im.crop((max(l - pad, 0), max(t - pad, 0), min(r + pad, im.width), min(b + pad, im.height))).save(path)


# ---------------------------------------------------------------------------
# Shaded triangle meshes
# ---------------------------------------------------------------------------

def shaded(tris, color, alpha=1.0, light=(-0.4, -0.6, 0.8), ambient=0.45):
    """Per-triangle Lambert shading (two-sided)."""
    from matplotlib.colors import to_rgb
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
    L = np.asarray(light, float)
    L /= np.linalg.norm(L)
    k = ambient + (1 - ambient) * np.abs(n @ L)
    rgb = np.clip(np.outer(k, to_rgb(color)), 0, 1)
    return np.column_stack([rgb, np.full(len(tris), alpha)])


def draw_meshes(ax, parts, zorder=1):
    """parts: list of (tris, color, alpha), drawn as one depth-sorted collection."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    polys = np.concatenate([p[0] for p in parts])
    colors = np.concatenate([shaded(*p) for p in parts])
    edges = colors.copy()
    edges[edges[:, 3] < 1, 3] = 0.0  # hide seams on opaque parts only
    pc = Poly3DCollection(polys, facecolors=colors, edgecolors=edges, linewidths=0.35, zorder=zorder)
    ax.add_collection3d(pc)
    return pc


def draw_pump(ax, b, top_sector=(np.deg2rad(20), np.deg2rad(200)), plate=0.03, n_theta=140):
    """Bottom plate, blade and a sector of the top plate, layered bottom to top.
    Seen from above, each layer is behind the next, so this order avoids
    depth-sorting errors between parts."""
    from blade_geometry import blade_mesh, disc_mesh
    blade = blade_mesh(b, n_theta=n_theta, n_w=12)
    R = np.linalg.norm(blade[..., :2], axis=-1).max() + 0.02
    layers = [[(disc_mesh(0.04, R, -plate, 0.0), PLATE_BOTTOM, 1.0)], [(blade, BLADE, 1.0)]]
    if top_sector is not None:
        eye = float(b.spiral(0.0)[0])
        layers.append([(disc_mesh(eye, R, b.height, b.height + plate, sector=top_sector), PLATE_TOP, 1.0)])
    for z, layer in enumerate(layers):
        draw_meshes(ax, layer, zorder=z + 1)
    return R
