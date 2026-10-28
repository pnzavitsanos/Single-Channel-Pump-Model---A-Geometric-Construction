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


# ---------------------------------------------------------------------------
# Shaded triangle meshes in one depth-sorted collection
# ---------------------------------------------------------------------------

HUB = "#8c959f"
BLADE = "#3987e5"
SHROUD = "#c9d1d9"
ACCENT = "#c2410c"


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


def draw_impeller(ax, d, shroud_sector=(np.deg2rad(20), np.deg2rad(200)), blade_color=None,
                  hub_thickness=0.03, shroud_thickness=0.02):
    """Hub, blade and a sector of the shroud, layered bottom to top. Seen from
    above, the hub is always behind the blade and the blade behind the shroud,
    so drawing them in that order avoids depth-sorting errors."""
    from impeller_geometry import blade_mesh, disc_mesh
    blade = blade_mesh(d, n_theta=140, n_w=12)
    R = np.linalg.norm(blade[..., :2], axis=-1).max() + 0.02
    hub = disc_mesh(0.04, R, -hub_thickness, 0.0)
    layers = [[(hub, HUB, 1.0)], [(blade, blade_color or BLADE, 1.0)]]
    if shroud_sector is not None:
        shroud = disc_mesh(d.r1, R, d.height, d.height + shroud_thickness, sector=shroud_sector)
        layers.append([(shroud, SHROUD, 1.0)])
    for z, layer in enumerate(layers):
        draw_meshes(ax, layer, zorder=z + 1)
    return R


def cutaway(tris, angle_from, angle_to):
    """Keep only the triangles whose centroid lies in the given polar-angle range."""
    c = tris.mean(axis=1)
    a = np.mod(np.arctan2(c[:, 1], c[:, 0]) - angle_from, 2 * np.pi)
    return tris[a <= np.mod(angle_to - angle_from, 2 * np.pi)]
