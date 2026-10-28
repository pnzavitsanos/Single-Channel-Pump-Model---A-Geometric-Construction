"""
impeller_geometry.py
=====================

Parametric single-blade impeller, built as a ruled surface over a plane curve.

Design parameters (8)
---------------------
Camber line in the hub plane (4):
    r1, r2        inlet and outlet radius
    beta1, beta2  inlet and outlet blade angle, measured from the
                  circumferential direction
Blade surface (3):
    H             blade height = distance between hub (z = 0) and shroud (z = H)
    lean1, lean2  inlet and outlet lean of the straight rulings from the axial
                  direction
Wall (1):
    h             blade thickness

The wrap angle Phi is not free: it follows from r1, r2, beta1, beta2.

Construction
------------
1.  Camber line.  The blade angle beta(theta) runs linearly from beta1 to
    beta2. A polar curve rho(theta) has blade angle tan(beta) = rho'/rho, so
        rho(theta) = r1 * (cos beta1 / cos beta(theta))^(1/k),  k = dbeta/dtheta,
    and rho = r1 e^{theta tan beta1} (a logarithmic spiral) when k = 0.
    The tangent direction is alpha = theta + pi/2 - beta, the speed is
    s = rho / cos(beta), and the curvature is
        kappa = alpha' / s = (1 - beta') cos(beta) / rho.

2.  Rulings.  At every point a straight segment is attached in the normal
    plane span{N, Z}, turned about the tangent T by the lean lambda(theta)
    (linear from lean1 to lean2) and stretched by 1/cos(lambda):
        X(theta, w) = P(theta) - w tan(lambda) N + w Z,   0 <= w <= H.
    The height of every point is exactly w, so every ruling runs from the hub
    plane to the shroud plane.

3.  Gaussian curvature (closed form):
        K = -(s lambda')^2 / [ (s + alpha' w tan lambda)^2 + (w lambda' / cos lambda)^2 ]^2.
    K = 0 everywhere (developable) exactly when the lean is constant, for any
    choice of blade angles.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, replace

import numpy as np

_NB_B = 0.2                      # notebook spiral growth rate b
_NB_PHI = 1.9 * np.pi            # notebook wrap angle


@dataclass(frozen=True)
class BladeDesign:
    """The eight design parameters. Defaults reproduce the Mathematica notebook."""

    r1: float = 0.2
    r2: float = 0.2 * float(np.exp(_NB_B * _NB_PHI))
    beta1_deg: float = float(np.degrees(np.arctan(_NB_B)))
    beta2_deg: float = float(np.degrees(np.arctan(_NB_B)))
    height: float = 0.3
    lean1_deg: float = 45.0
    lean2_deg: float = float(45.0 - np.degrees(_NB_PHI / 8))
    thickness: float = 0.05

    def __post_init__(self):
        if not 0 < self.r1 < self.r2:
            raise ValueError("need 0 < r1 < r2")
        for name in ("beta1_deg", "beta2_deg"):
            if not 0 < getattr(self, name) < 90:
                raise ValueError(f"{name} must lie strictly between 0 and 90 degrees")
        for name in ("lean1_deg", "lean2_deg"):
            if not -90 < getattr(self, name) < 90:
                raise ValueError(f"{name} must lie strictly between -90 and 90 degrees: "
                                 "a horizontal ruling can never reach the shroud")
        if self.height <= 0 or self.thickness < 0:
            raise ValueError("height must be positive and thickness non-negative")
        if self.blade_angle_rate >= 1:
            raise ValueError("blade angle grows too fast: the camber line would have an inflection")

    def with_(self, **changes) -> "BladeDesign":
        return replace(self, **changes)

    # --- derived quantities ---------------------------------------------
    @property
    def beta1(self):
        return np.deg2rad(self.beta1_deg)

    @property
    def beta2(self):
        return np.deg2rad(self.beta2_deg)

    @property
    def wrap(self) -> float:
        """Wrap angle Phi (rad), fixed by r1, r2, beta1, beta2."""
        L = np.log(self.r2 / self.r1)
        d = self.beta2 - self.beta1
        if abs(d) < 1e-9:
            return L / np.tan(0.5 * (self.beta1 + self.beta2))
        return d * L / np.log(np.cos(self.beta1) / np.cos(self.beta2))

    @property
    def wrap_deg(self) -> float:
        return float(np.degrees(self.wrap))

    @property
    def blade_angle_rate(self) -> float:
        """k = dbeta/dtheta (rad per rad)."""
        L = np.log(self.r2 / self.r1)
        d = self.beta2 - self.beta1
        if abs(d) < 1e-9:
            return 0.0
        return float(np.log(np.cos(self.beta1) / np.cos(self.beta2)) / L)

    @property
    def lean_rate(self) -> float:
        """m = dlambda/dtheta (rad per rad)."""
        return float(np.deg2rad(self.lean2_deg - self.lean1_deg) / self.wrap)


NOTEBOOK = BladeDesign()


# ---------------------------------------------------------------------------
# Laws along the blade
# ---------------------------------------------------------------------------

def blade_angle(theta, d: BladeDesign = NOTEBOOK):
    theta = np.asarray(theta, float)
    return d.beta1 + d.blade_angle_rate * theta


def lean(theta, d: BladeDesign = NOTEBOOK):
    theta = np.asarray(theta, float)
    return np.deg2rad(d.lean1_deg) + d.lean_rate * theta


def radius(theta, d: BladeDesign = NOTEBOOK):
    """rho(theta) = r1 exp( integral of tan(beta) ) in closed form."""
    theta = np.asarray(theta, float)
    k = d.blade_angle_rate
    if abs(k) < 1e-12:
        return d.r1 * np.exp(theta * np.tan(d.beta1))
    return d.r1 * (np.cos(d.beta1) / np.cos(blade_angle(theta, d))) ** (1.0 / k)


def frame(theta, d: BladeDesign = NOTEBOOK):
    """Point P, tangent T, normal N (shape (3, ...)), speed s and turning rate alpha'."""
    theta = np.asarray(theta, float)
    rho = radius(theta, d)
    beta = blade_angle(theta, d)
    alpha = theta + np.pi / 2 - beta
    zero = np.zeros_like(theta)
    P = np.stack([rho * np.cos(theta), rho * np.sin(theta), zero])
    T = np.stack([np.cos(alpha), np.sin(alpha), zero])
    N = np.stack([-np.sin(alpha), np.cos(alpha), zero])
    s = rho / np.cos(beta)
    dalpha = np.full_like(theta, 1.0 - d.blade_angle_rate)
    return P, T, N, s, dalpha


def camber_curvature(theta, d: BladeDesign = NOTEBOOK):
    """kappa = (1 - beta') cos(beta) / rho."""
    _, _, _, s, dalpha = frame(theta, d)
    return dalpha / s


def _Z(like):
    Z = np.zeros_like(like)
    Z[2] = 1.0
    return Z


def surface(theta, w, d: BladeDesign = NOTEBOOK):
    """X(theta, w) = P - w tan(lambda) N + w Z."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    P, _, N, _, _ = frame(theta, d)
    return P - w * np.tan(lean(theta, d)) * N + w * _Z(P)


def normal(theta, w, d: BladeDesign = NOTEBOOK):
    """X_theta x X_w = B T - A N - A tan(lambda) Z."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, T, N, s, dalpha = frame(theta, d)
    lam, m = lean(theta, d), d.lean_rate
    A = s + dalpha * w * np.tan(lam)
    B = -w * m / np.cos(lam) ** 2
    return B * T - A * N - A * np.tan(lam) * _Z(T)


def gaussian_curvature(theta, w, d: BladeDesign = NOTEBOOK):
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, _, _, s, dalpha = frame(theta, d)
    lam, m = lean(theta, d), d.lean_rate
    denom = (s + dalpha * w * np.tan(lam)) ** 2 + (w * m / np.cos(lam)) ** 2
    return -(s * m) ** 2 / denom ** 2


def level_curve(c: float, theta, d: BladeDesign = NOTEBOOK):
    """Slice at height c: the offset curve P - c tan(lambda) N."""
    theta = np.asarray(theta, float)
    return surface(theta, np.full_like(theta, c), d)


def slice_blade_angle(c: float, theta, d: BladeDesign = NOTEBOOK):
    """Blade angle (rad) of the slice at height c, from the circumferential direction."""
    theta = np.asarray(theta, float)
    Q = level_curve(c, theta, d)
    _, T, N, s, dalpha = frame(theta, d)
    lam = lean(theta, d)
    tangent = (s + dalpha * c * np.tan(lam)) * T - c * d.lean_rate / np.cos(lam) ** 2 * N
    e_r = Q[:2] / np.linalg.norm(Q[:2], axis=0)
    e_phi = np.stack([-e_r[1], e_r[0]])
    return np.arctan2((tangent[:2] * e_r).sum(0), (tangent[:2] * e_phi).sum(0))


# ---------------------------------------------------------------------------
# Grids, wall faces, and the impeller solid
# ---------------------------------------------------------------------------

def grid(d: BladeDesign = NOTEBOOK, n_theta: int = 160, n_w: int = 24):
    th = np.linspace(0.0, d.wrap, n_theta)
    wv = np.linspace(0.0, d.height, n_w)
    return th, wv, *np.meshgrid(th, wv, indexing="ij")


def blade_face(sign: float, d: BladeDesign = NOTEBOOK, n_theta: int = 160, n_w: int = 24,
               offset: str = "horizontal"):
    """One wall face, offset by sign * h/2.

    offset="horizontal": inside each horizontal slice, so the face stays in 0 <= z <= H
    offset="normal":     along the unit surface normal (the notebook's CAD export)
    """
    _, _, TH, W = grid(d, n_theta, n_w)
    n = normal(TH, W, d)
    if offset == "horizontal":
        n[2] = 0.0
    elif offset != "normal":
        raise ValueError("offset must be 'horizontal' or 'normal'")
    return surface(TH, W, d) + sign * d.thickness / 2 * n / np.linalg.norm(n, axis=0)


def _tube(loops):
    """Triangles joining consecutive closed loops, shape (n_loops, n_pts, 3)."""
    tris = []
    n_l, n_p, _ = loops.shape
    for i in range(n_l - 1):
        a, b = loops[i], loops[i + 1]
        for k in range(n_p):
            k2 = (k + 1) % n_p
            tris.append((a[k], b[k], b[k2]))
            tris.append((a[k], b[k2], a[k2]))
    return np.array(tris)


def _orient_outward(tris):
    vol = np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6
    return tris if vol > 0 else tris[:, ::-1]


def blade_mesh(d: BladeDesign = NOTEBOOK, n_theta: int = 160, n_w: int = 16):
    """Closed triangle mesh of the blade with thickness (horizontal offset)."""
    A = np.moveaxis(blade_face(-1, d, n_theta, n_w), 0, -1)   # (n_theta, n_w, 3)
    B = np.moveaxis(blade_face(+1, d, n_theta, n_w), 0, -1)
    loops = np.concatenate([A, B[:, ::-1]], axis=1)            # closed cross-sections
    side = _tube(loops)
    caps = []
    for i, flip in ((0, True), (-1, False)):
        for j in range(n_w - 1):
            t1 = (A[i, j], A[i, j + 1], B[i, j + 1])
            t2 = (A[i, j], B[i, j + 1], B[i, j])
            caps += [t1[::-1], t2[::-1]] if flip else [t1, t2]
    tris = np.concatenate([side, np.array(caps)])
    # make the caps agree with the side walls, then point everything outward
    if not _edges_consistent(tris):
        caps = np.array(caps)[:, ::-1]
        tris = np.concatenate([side, caps])
    return _orient_outward(tris)


def _edges_consistent(tris):
    seen = set()
    for t in tris:
        for e in ((0, 1), (1, 2), (2, 0)):
            key = (tuple(np.round(t[e[0]], 12)), tuple(np.round(t[e[1]], 12)))
            if key in seen:
                return False
            seen.add(key)
    return True


def disc_mesh(r_in: float, r_out: float, z0: float, z1: float, n: int = 120, n_r: int = 14,
              sector: tuple[float, float] | None = None):
    """Closed annular disc (a revolved rectangle), split into rings for clean rendering.
    With sector=(a0, a1) only that angular piece is built, closed by flat end faces."""
    rr = np.linspace(r_in, r_out, n_r)
    profile = np.array([(r, z0) for r in rr] + [(r, z1) for r in rr[::-1]])

    def loop(a):
        return np.array([(r * np.cos(a), r * np.sin(a), z) for r, z in profile])

    if sector is None:
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        loops = np.array([loop(a) for a in t] + [loop(0.0)])
        return _orient_outward(_tube(loops))
    t = np.linspace(sector[0], sector[1], max(int(n * (sector[1] - sector[0]) / (2 * np.pi)), 2))
    loops = np.array([loop(a) for a in t])
    side = _tube(loops)
    caps = []
    for L, flip in ((loops[0], True), (loops[-1], False)):
        bot, top = L[:n_r], L[n_r:][::-1]
        for j in range(n_r - 1):
            t1, t2 = (bot[j], bot[j + 1], top[j + 1]), (bot[j], top[j + 1], top[j])
            caps += [t1[::-1], t2[::-1]] if flip else [t1, t2]
    tris = np.concatenate([side, np.array(caps)])
    if not _edges_consistent(tris):
        tris = np.concatenate([side, np.array(caps)[:, ::-1]])
    return _orient_outward(tris)


def impeller_meshes(d: BladeDesign = NOTEBOOK, hub_thickness: float = 0.03,
                    shroud_thickness: float = 0.02, bore: float = 0.04, eye: float | None = None):
    """Blade, hub disc and shroud disc as closed triangle meshes."""
    blade = blade_mesh(d)
    r_max = np.linalg.norm(blade[..., :2], axis=-1).max() + 0.02
    eye = d.r1 if eye is None else eye
    hub = disc_mesh(bore, r_max, -hub_thickness, 0.0)
    shroud = disc_mesh(eye, r_max, d.height, d.height + shroud_thickness)
    return {"blade": blade, "hub": hub, "shroud": shroud}


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
# Independent check against a symbolic derivation (needs SymPy)
# ---------------------------------------------------------------------------

def symbolic_check(d: BladeDesign | None = None):
    """Rebuild the construction symbolically and compare with the closed forms."""
    import sympy as sp

    d = d or NOTEBOOK.with_(beta2_deg=25.0, lean2_deg=-10.0)
    th, v = sp.symbols("theta v", real=True)
    k, m = sp.Float(d.blade_angle_rate), sp.Float(d.lean_rate)
    beta = d.beta1 + k * th
    lam = sp.Float(np.deg2rad(d.lean1_deg)) + m * th
    rho = d.r1 * sp.exp(sp.integrate(sp.tan(d.beta1 + k * sp.Symbol("u")), (sp.Symbol("u"), 0, th)))
    P = sp.Matrix([rho * sp.cos(th), rho * sp.sin(th), 0])
    Pp = P.diff(th)
    s = sp.sqrt(Pp.dot(Pp))
    T = Pp / s
    Tp = T.diff(th)
    N = Tp / sp.sqrt(Tp.dot(Tp))
    F = -sp.sin(lam) * N + sp.cos(lam) * sp.Matrix([0, 0, 1])
    X = P + v * F
    Xt, Xv = X.diff(th), X.diff(v)
    n = Xt.cross(Xv)
    K_expr = (Xt.diff(th).dot(n) * Xv.diff(v).dot(n) - Xt.diff(v).dot(n) ** 2) / n.dot(n) ** 2
    kappa_expr = sp.sqrt(Tp.dot(Tp)) / s
    radial = (Pp.dot(P) / sp.sqrt(P.dot(P)))
    circ = (Pp[1] * P[0] - Pp[0] * P[1]) / sp.sqrt(P.dot(P))
    f = sp.lambdify((th, v), [K_expr, kappa_expr, sp.atan2(radial, circ), rho, (P + v * F)[2]], "numpy")

    errs = {"K": 0.0, "kappa": 0.0, "blade angle": 0.0, "radius": 0.0, "height": 0.0}
    for t0 in np.linspace(0.05, d.wrap - 0.05, 5):
        for v0 in (0.05, 0.15, 0.3):
            K0, kap0, b0, rho0, z0 = f(t0, v0)
            w0 = v0 * np.cos(float(lean(t0, d)))
            errs["K"] = max(errs["K"], abs(K0 - gaussian_curvature(t0, w0, d)))
            errs["kappa"] = max(errs["kappa"], abs(kap0 - camber_curvature(t0, d)))
            errs["blade angle"] = max(errs["blade angle"], abs(b0 - blade_angle(t0, d)))
            errs["radius"] = max(errs["radius"], abs(rho0 - radius(t0, d)))
            errs["height"] = max(errs["height"], abs(z0 - w0))  # height = w
    print(f"test design: beta {d.beta1_deg:.2f} -> {d.beta2_deg:.2f} deg, "
          f"lean {d.lean1_deg:.2f} -> {d.lean2_deg:.2f} deg, wrap {d.wrap_deg:.1f} deg")
    for key, val in errs.items():
        print(f"  max error, {key:12s}: {val:.1e}")
    print(f"  radius at the outlet: {radius(d.wrap, d):.6f} (r2 = {d.r2:.6f})")


def summary(d: BladeDesign = NOTEBOOK):
    th, wv, TH, W = grid(d)
    K = gaussian_curvature(TH, W, d)
    rows = [
        ("inlet radius r1", f"{d.r1:.4f}"),
        ("outlet radius r2", f"{d.r2:.4f}"),
        ("inlet blade angle beta1", f"{d.beta1_deg:.2f} deg"),
        ("outlet blade angle beta2", f"{d.beta2_deg:.2f} deg"),
        ("blade height H", f"{d.height:.4f}"),
        ("inlet lean lambda1", f"{d.lean1_deg:.2f} deg"),
        ("outlet lean lambda2", f"{d.lean2_deg:.2f} deg"),
        ("thickness h", f"{d.thickness:.4f}"),
        ("-> wrap angle", f"{d.wrap_deg:.1f} deg"),
        ("-> blade angle at shroud, inlet/outlet",
         f"{np.degrees(slice_blade_angle(d.height, 0.0, d)):.2f} / "
         f"{np.degrees(slice_blade_angle(d.height, d.wrap, d)):.2f} deg"),
        ("-> Gaussian curvature range", f"{K.min():.4f} to {K.max():.4f}"),
        ("-> developable", "yes" if np.abs(K).max() < 1e-12 else "no"),
    ]
    for k_, v_ in rows:
        print(f"  {k_:40s} {v_}")


if __name__ == "__main__":
    print("Notebook design:")
    summary(NOTEBOOK)
    print()
    symbolic_check()
