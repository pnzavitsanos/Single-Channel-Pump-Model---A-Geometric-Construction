"""
blade_geometry.py
==================

The blade of a single-channel pump as a ruled surface over a spiral.

Construction
------------
1.  Base curve: any spiral rho(theta) in the plane z = 0, for 0 <= theta <= Phi.
    Logarithmic, Archimedean and arbitrary formulas are provided. Everything
    below needs only rho, rho' and rho'':
        P = rho (cos theta, sin theta, 0),   s = |P'| = sqrt(rho^2 + rho'^2)
        T = P'/s,   N = Z x T   (N points to the inside of the turn)
        spiral angle  psi = arctan(rho'/rho)   (between the curve and the circle)
        turning rate  alpha' = 1 - psi',   curvature kappa = alpha'/s
                    = (rho^2 + 2 rho'^2 - rho rho'') / (rho^2 + rho'^2)^(3/2)

2.  Rulings: at every point a straight segment is attached in the normal plane
    span{N, Z}, turned about the tangent T by the angle beta(theta), which runs
    linearly from the inlet angle beta_in to the outlet angle beta_out, and
    stretched by 1/cos(beta):
        X(theta, w) = P + (w / cos beta) (-sin beta N + cos beta Z)
                    = P - w tan(beta) N + w Z,        0 <= w <= H.
    The height of every point is exactly w, so every ruling runs from the
    plane z = 0 to the plane z = H.

3.  Gaussian curvature, in closed form for every base curve:
        K = -(s beta')^2 / [ (s + alpha' w tan beta)^2 + (w beta' / cos beta)^2 ]^2,
    so the blade is developable exactly when beta_in = beta_out.

Parameters: the base spiral (2 shape numbers for the log and Archimedean
families, plus the angular extent Phi), the inlet and outlet angles, the
height H, and the wall thickness h.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field, replace

import numpy as np


# ---------------------------------------------------------------------------
# Base curves: anything that returns rho, rho', rho''
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LogSpiral:
    """rho = a e^{b theta}: a sets the start radius, b how fast the spiral opens."""
    a: float = 0.2
    b: float = 0.2

    def __call__(self, th):
        r = self.a * np.exp(self.b * th)
        return r, self.b * r, self.b ** 2 * r

    def formula(self):
        return rf"$\rho = {self.a:g}\,e^{{{self.b:.3g}\theta}}$"


@dataclass(frozen=True)
class ArchimedeanSpiral:
    """rho = a + b theta: the radius grows by the same amount every turn."""
    a: float = 0.2
    b: float = 0.0771

    def __call__(self, th):
        th = np.asarray(th, float)
        return self.a + self.b * th, np.full_like(th, self.b), np.zeros_like(th)

    def formula(self):
        return rf"$\rho = {self.a:g} + {self.b:.3g}\,\theta$"


@dataclass(frozen=True)
class CustomSpiral:
    """Any spiral given as a formula in theta, e.g. "0.2 + 0.05*theta + 0.01*theta**2"."""
    expression: str = "0.2*exp(0.2*theta)"
    _f: tuple = field(default=None, compare=False, repr=False)

    def __post_init__(self):
        import sympy as sp
        t = sp.Symbol("theta", real=True)
        e = sp.sympify(self.expression, locals={"theta": t})
        fs = tuple(sp.lambdify(t, sp.diff(e, t, k), "numpy") for k in range(3))
        object.__setattr__(self, "_f", fs)

    def __call__(self, th):
        th = np.asarray(th, float)
        return tuple(np.broadcast_to(f(th), th.shape).astype(float) for f in self._f)

    def formula(self):
        body = self.expression.replace("**", "^").replace("*", r"\,").replace("theta", r"\theta")
        return rf"$\rho = {body}$"


@dataclass(frozen=True)
class Blend:
    """(1 - t) * first + t * second: a continuous path between two spirals."""
    first: object
    second: object
    t: float

    def __call__(self, th):
        return tuple((1 - self.t) * u + self.t * v for u, v in zip(self.first(th), self.second(th)))

    def formula(self):
        return f"{1 - self.t:.2f} × log + {self.t:.2f} × Archimedean"


# ---------------------------------------------------------------------------
# The blade
# ---------------------------------------------------------------------------

_NB_WRAP_DEG = 1.9 * 180.0


@dataclass(frozen=True)
class Blade:
    """Defaults reproduce the final Mathematica notebook."""

    spiral: object = LogSpiral(0.2, 0.2)
    wrap_deg: float = _NB_WRAP_DEG          # angular extent Phi of the base curve
    inlet_angle_deg: float = 45.0           # beta_in, ruling angle at the inlet
    outlet_angle_deg: float = 45.0 - _NB_WRAP_DEG / 8   # beta_out (2.25 deg)
    height: float = 0.3                     # H, distance between the two planes
    thickness: float = 0.05                 # h

    def __post_init__(self):
        for name in ("inlet_angle_deg", "outlet_angle_deg"):
            if not -90 < getattr(self, name) < 90:
                raise ValueError(f"{name} must lie strictly between -90 and 90 degrees: "
                                 "a horizontal ruling can never reach the top plane")
        if self.wrap_deg <= 0 or self.height <= 0 or self.thickness < 0:
            raise ValueError("wrap and height must be positive, thickness non-negative")

    def with_(self, **changes) -> "Blade":
        return replace(self, **changes)

    @property
    def wrap(self) -> float:
        return float(np.radians(self.wrap_deg))

    @property
    def angle_rate(self) -> float:
        """beta' (rad per rad)."""
        return float(np.radians(self.outlet_angle_deg - self.inlet_angle_deg) / self.wrap)


NOTEBOOK = Blade()


def ruling_angle(theta, b: Blade = NOTEBOOK):
    """beta(theta), linear from the inlet angle to the outlet angle."""
    return np.radians(b.inlet_angle_deg) + b.angle_rate * np.asarray(theta, float)


def frame(theta, b: Blade = NOTEBOOK):
    """Point P, tangent T, normal N (shape (3, ...)), speed s and turning rate alpha'."""
    theta = np.asarray(theta, float)
    r, r1, r2 = b.spiral(theta)
    c, sn = np.cos(theta), np.sin(theta)
    zero = np.zeros_like(theta)
    P = np.stack([r * c, r * sn, zero])
    dP = np.stack([r1 * c - r * sn, r1 * sn + r * c, zero])
    s = np.sqrt(r * r + r1 * r1)
    T = dP / s
    N = np.stack([-T[1], T[0], zero])
    dalpha = (r * r + 2 * r1 * r1 - r * r2) / (r * r + r1 * r1)
    return P, T, N, s, dalpha


def spiral_angle(theta, b: Blade = NOTEBOOK):
    """psi = arctan(rho'/rho): angle between the base curve and the circle through the point."""
    r, r1, _ = b.spiral(np.asarray(theta, float))
    return np.arctan2(r1, r)


def curvature(theta, b: Blade = NOTEBOOK):
    """kappa = (rho^2 + 2 rho'^2 - rho rho'') / (rho^2 + rho'^2)^(3/2)."""
    _, _, _, s, dalpha = frame(theta, b)
    return dalpha / s


def _Z(like):
    Z = np.zeros_like(like)
    Z[2] = 1.0
    return Z


def surface(theta, w, b: Blade = NOTEBOOK):
    """X(theta, w) = P - w tan(beta) N + w Z."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    P, _, N, _, _ = frame(theta, b)
    return P - w * np.tan(ruling_angle(theta, b)) * N + w * _Z(P)


def normal(theta, w, b: Blade = NOTEBOOK):
    """X_theta x X_w = B T - A N - A tan(beta) Z,  A = s + alpha' w tan(beta),  B = -w beta' / cos^2(beta)."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, T, N, s, dalpha = frame(theta, b)
    beta = ruling_angle(theta, b)
    A = s + dalpha * w * np.tan(beta)
    B = -w * b.angle_rate / np.cos(beta) ** 2
    return B * T - A * N - A * np.tan(beta) * _Z(T)


def gaussian_curvature(theta, w, b: Blade = NOTEBOOK):
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, _, _, s, dalpha = frame(theta, b)
    beta, m = ruling_angle(theta, b), b.angle_rate
    denom = (s + dalpha * w * np.tan(beta)) ** 2 + (w * m / np.cos(beta)) ** 2
    return -(s * m) ** 2 / denom ** 2


def level_curve(c: float, theta, b: Blade = NOTEBOOK):
    """Slice at height c: the offset curve P - c tan(beta) N."""
    theta = np.asarray(theta, float)
    return surface(theta, np.full_like(theta, c), b)


def slice_spiral_angle(c: float, theta, b: Blade = NOTEBOOK):
    """Spiral angle (rad) of the slice at height c, against the circle through each of its points."""
    theta = np.asarray(theta, float)
    Q = level_curve(c, theta, b)
    _, T, N, s, dalpha = frame(theta, b)
    beta = ruling_angle(theta, b)
    tangent = (s + dalpha * c * np.tan(beta)) * T - c * b.angle_rate / np.cos(beta) ** 2 * N
    e_r = Q[:2] / np.linalg.norm(Q[:2], axis=0)
    e_phi = np.stack([-e_r[1], e_r[0]])
    return np.arctan2((tangent[:2] * e_r).sum(0), (tangent[:2] * e_phi).sum(0))


# ---------------------------------------------------------------------------
# Grids, wall faces, solids and STL
# ---------------------------------------------------------------------------

def grid(b: Blade = NOTEBOOK, n_theta: int = 160, n_w: int = 24):
    th = np.linspace(0.0, b.wrap, n_theta)
    wv = np.linspace(0.0, b.height, n_w)
    return th, wv, *np.meshgrid(th, wv, indexing="ij")


def blade_face(sign: float, b: Blade = NOTEBOOK, n_theta: int = 160, n_w: int = 24,
               offset: str = "horizontal"):
    """One wall face, offset by sign * h/2.

    offset="horizontal": inside each horizontal slice, so the face stays in 0 <= z <= H
    offset="normal":     along the unit surface normal (the notebook's CAD export)
    """
    _, _, TH, W = grid(b, n_theta, n_w)
    n = normal(TH, W, b)
    if offset == "horizontal":
        n[2] = 0.0
    elif offset != "normal":
        raise ValueError("offset must be 'horizontal' or 'normal'")
    return surface(TH, W, b) + sign * b.thickness / 2 * n / np.linalg.norm(n, axis=0)


def _tube(loops, closed=True):
    tris = []
    n_l, n_p, _ = loops.shape
    for i in range(n_l - 1):
        p, q = loops[i], loops[i + 1]
        for k in range(n_p if closed else n_p - 1):
            k2 = (k + 1) % n_p
            tris.append((p[k], q[k], q[k2]))
            tris.append((p[k], q[k2], p[k2]))
    return np.array(tris)


def _edges_consistent(tris):
    seen = set()
    for t in np.round(tris, 12):
        for i, j in ((0, 1), (1, 2), (2, 0)):
            key = (tuple(t[i]), tuple(t[j]))
            if key in seen:
                return False
            seen.add(key)
    return True


def _orient_outward(tris):
    vol = np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6
    return tris if vol > 0 else tris[:, ::-1]


def _with_caps(side, first, second, n_half):
    """Close a tube whose cross-section is two rows of n_half points."""
    caps = []
    for L, flip in ((first, True), (second, False)):
        lo, hi = L[:n_half], L[n_half:][::-1]
        for j in range(n_half - 1):
            t1, t2 = (lo[j], lo[j + 1], hi[j + 1]), (lo[j], hi[j + 1], hi[j])
            caps += [t1[::-1], t2[::-1]] if flip else [t1, t2]
    tris = np.concatenate([side, np.array(caps)])
    if not _edges_consistent(tris):
        tris = np.concatenate([side, np.array(caps)[:, ::-1]])
    return _orient_outward(tris)


def blade_mesh(b: Blade = NOTEBOOK, n_theta: int = 160, n_w: int = 16):
    """Closed triangle mesh of the blade with its wall thickness."""
    A = np.moveaxis(blade_face(-1, b, n_theta, n_w), 0, -1)
    B = np.moveaxis(blade_face(+1, b, n_theta, n_w), 0, -1)
    loops = np.concatenate([A, B[:, ::-1]], axis=1)
    return _with_caps(_tube(loops), loops[0], loops[-1], n_w)


def disc_mesh(r_in: float, r_out: float, z0: float, z1: float, n: int = 120, n_r: int = 14,
              sector: tuple[float, float] | None = None):
    """Closed annular disc; with sector=(a0, a1) only that angular piece, with flat ends."""
    rr = np.linspace(r_in, r_out, n_r)
    profile = np.array([(r, z0) for r in rr] + [(r, z1) for r in rr[::-1]])

    def loop(a):
        return np.array([(r * np.cos(a), r * np.sin(a), z) for r, z in profile])

    if sector is None:
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        return _orient_outward(_tube(np.array([loop(a) for a in t] + [loop(0.0)])))
    t = np.linspace(sector[0], sector[1], max(int(n * (sector[1] - sector[0]) / (2 * np.pi)), 2))
    loops = np.array([loop(a) for a in t])
    return _with_caps(_tube(loops), loops[0], loops[-1], n_r)


def pump_meshes(b: Blade = NOTEBOOK, plate: float = 0.03, bore: float = 0.04, eye: float | None = None):
    """Blade plus the two parallel plates (bottom with a shaft bore, top with an inlet eye)."""
    blade = blade_mesh(b)
    r_max = np.linalg.norm(blade[..., :2], axis=-1).max() + 0.02
    eye = float(b.spiral(0.0)[0]) if eye is None else eye
    return {"blade": blade,
            "bottom_plate": disc_mesh(bore, r_max, -plate, 0.0),
            "top_plate": disc_mesh(eye, r_max, b.height, b.height + plate)}


def write_stl(tris, path: str, name: str = "part"):
    """Binary STL."""
    tris = np.asarray(tris, dtype=np.float32)
    nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        for n, t in zip(nrm, tris):
            f.write(struct.pack("<12fH", *n, *t.ravel(), 0))


# ---------------------------------------------------------------------------
# Summary and an independent symbolic check
# ---------------------------------------------------------------------------

def summary(b: Blade = NOTEBOOK):
    th, wv, TH, W = grid(b)
    K = gaussian_curvature(TH, W, b)
    ends = np.array([0.0, b.wrap])
    rho = b.spiral(ends)[0]
    psi = np.degrees(spiral_angle(ends, b))
    psi_top = np.degrees(slice_spiral_angle(b.height, ends, b))
    kap = curvature(ends, b)
    rows = [
        ("base curve", b.spiral.formula().strip("$").replace(r"\,", " ").replace(r"\theta", "θ")
                       .replace(r"\rho", "ρ")),
        ("angular extent Phi", f"{b.wrap_deg:.1f} deg"),
        ("inlet angle beta_in", f"{b.inlet_angle_deg:.2f} deg"),
        ("outlet angle beta_out", f"{b.outlet_angle_deg:.2f} deg"),
        ("height H", f"{b.height:.4f}"),
        ("thickness h", f"{b.thickness:.4f}"),
        ("-> radius inlet / outlet", f"{rho[0]:.4f} / {rho[1]:.4f}"),
        ("-> spiral angle inlet / outlet (bottom)", f"{psi[0]:.2f} / {psi[1]:.2f} deg"),
        ("-> spiral angle inlet / outlet (top)", f"{psi_top[0]:.2f} / {psi_top[1]:.2f} deg"),
        ("-> curvature of base curve inlet / outlet", f"{kap[0]:.3f} / {kap[1]:.3f}"),
        ("-> Gaussian curvature range", f"{K.min() + 0.0:.4f} to {K.max() + 0.0:.4f}"),
        ("-> developable", "yes" if np.abs(K).max() < 1e-12 else "no"),
    ]
    for k, v in rows:
        print(f"  {k:42s} {v}")


def symbolic_check():
    """Rebuild the construction from the definitions with SymPy, for three base
    curves, and compare with the closed forms used above."""
    import sympy as sp

    th, v = sp.symbols("theta v", real=True)
    cases = [
        ("log", "0.2*exp(0.2*theta)", LogSpiral(0.2, 0.2)),
        ("Archimedean", "0.2 + 0.0771*theta", ArchimedeanSpiral(0.2, 0.0771)),
        ("custom", "0.2 + 0.05*theta + 0.01*theta**2", CustomSpiral("0.2 + 0.05*theta + 0.01*theta**2")),
    ]
    for name, expr, spiral in cases:
        b = Blade(spiral=spiral, inlet_angle_deg=40.0, outlet_angle_deg=-10.0)
        rho = sp.sympify(expr, locals={"theta": th})
        P = sp.Matrix([rho * sp.cos(th), rho * sp.sin(th), 0])
        Pp = P.diff(th)
        s = sp.sqrt(Pp.dot(Pp))
        T = Pp / s
        Tp = T.diff(th)
        N = Tp / sp.sqrt(Tp.dot(Tp))
        beta = sp.Float(np.radians(b.inlet_angle_deg)) + sp.Float(b.angle_rate) * th
        F = -sp.sin(beta) * N + sp.cos(beta) * sp.Matrix([0, 0, 1])
        X = P + v * F
        Xt, Xv = X.diff(th), X.diff(v)
        n = Xt.cross(Xv)
        K = (Xt.diff(th).dot(n) * Xv.diff(v).dot(n) - Xt.diff(v).dot(n) ** 2) / n.dot(n) ** 2
        kappa = sp.sqrt(Tp.dot(Tp)) / s
        f = sp.lambdify((th, v), [K, kappa, (P + v / sp.cos(beta) * F)[2]], "numpy")
        err = np.zeros(3)
        for t0 in np.linspace(0.1, b.wrap - 0.1, 6):
            for v0 in (0.05, 0.2, 0.35):
                K0, k0, z0 = f(t0, v0)
                w0 = v0 * np.cos(float(ruling_angle(t0, b)))
                err = np.maximum(err, [abs(K0 - gaussian_curvature(t0, w0, b)),
                                       abs(k0 - curvature(t0, b)), abs(z0 - v0)])
        print(f"  {name:12s} max error  K {err[0]:.1e}   curvature {err[1]:.1e}   height {err[2]:.1e}")


if __name__ == "__main__":
    print("Notebook design:")
    summary(NOTEBOOK)
    print("\nSame inlet/outlet angles on an Archimedean spiral:")
    summary(NOTEBOOK.with_(spiral=ArchimedeanSpiral(0.2, 0.0771)))
    print("\nClosed forms vs symbolic derivation:")
    symbolic_check()
