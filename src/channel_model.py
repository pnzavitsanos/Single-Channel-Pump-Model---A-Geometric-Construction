"""
channel_model.py
=================

Closed-form model of the blade surface developed in
``notebooks/05_log_spiral_v2_final.nb``.

The construction
----------------
1.  Base curve: a logarithmic spiral in the plane z = 0,
        r(theta) = a e^{b theta} (cos theta, sin theta, 0).
    Its tangent turns at exactly unit rate in theta (kappa * s = 1), so the
    Frenet frame is explicit:
        T(theta) = ( cos(theta + phi),  sin(theta + phi), 0)
        N(theta) = (-sin(theta + phi),  cos(theta + phi), 0),   cot(phi) = b,
    with speed s(theta) = a sqrt(1 + b^2) e^{b theta}, which is also the
    radius of curvature.

2.  Rulings: at each point of the spiral a straight segment is attached in
    the normal plane span{N, Z}, leaning from the vertical by beta(theta):
        F(theta) = -sin(beta) N + cos(beta) Z.
    Because N x Z = T, the ruling is rotated about the tangent: the Moebius
    mechanism, transported along a spiral.

3.  Two parallel planes: each ruling is scaled by 1 / cos(beta), so
        X(theta, w) = r + (w / cos beta) F = r - w tan(beta) N + w Z.
    The height of every point is exactly w, so for w in [0, H] each ruling
    runs from the plane z = 0 to the plane z = H, whatever its lean. The
    slice at height c is the offset curve r - c tan(beta) N.

4.  Curvature: the surface is ruled and
        K = -(s beta')^2 / [ (s + w tan beta)^2 + (w beta' / cos beta)^2 ]^2,
    so it is developable (K = 0, can be rolled from flat sheet) exactly when
    the lean is constant, and saddle-shaped (K < 0) wherever the lean changes.

5.  Wall thickness: the notebook offsets the surface by +/- h/2 along its
    unit normal. Because the normal is tilted, that moves the edges slightly
    off the two planes; the "horizontal" option offsets each slice inside its
    own plane instead, so both faces stay exactly between the planes.

The lean law of the notebook is beta(theta) = pi/4 - theta/8. The general law
used for exploration is linear between an inlet and an outlet lean; both must
stay strictly between -90 and 90 degrees, since a horizontal ruling can never
reach the top plane.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpiralParams:
    """Geometric parameters of the blade (defaults are the notebook values)."""

    a: float = 0.2                # spiral scale
    b: float = 0.2                # spiral growth rate
    height: float = 0.3           # distance H between the two parallel planes
    wall_thickness: float = 0.05  # h


DEFAULT_PARAMS = SpiralParams()

# Angular extent of the blade in the notebook's final render and CAD export.
BASELINE_INLET_DEG = 0.0
BASELINE_OUTLET_DEG = 1.9 * 180.0  # 1.9*pi rad


# ---------------------------------------------------------------------------
# Lean laws: each maps theta (rad) to (beta, dbeta/dtheta) in radians
# ---------------------------------------------------------------------------

def notebook_lean(theta):
    """beta(theta) = pi/4 - theta/8, the law used in the final notebook."""
    theta = np.asarray(theta, dtype=float)
    return np.pi / 4 - theta / 8, np.full_like(theta, -1.0 / 8)


def linear_lean(theta_max_deg: float, beta_in_deg: float, beta_out_deg: float):
    """Lean law that runs linearly from beta_in at theta = 0 to beta_out at theta_max."""
    for name, val in (("beta_in_deg", beta_in_deg), ("beta_out_deg", beta_out_deg)):
        if not -90.0 < val < 90.0:
            raise ValueError(f"{name} = {val}: the lean must stay strictly between -90 and 90 "
                             "degrees, otherwise the ruling cannot reach the top plane.")
    b_in = np.deg2rad(beta_in_deg)
    slope = (np.deg2rad(beta_out_deg) - b_in) / np.deg2rad(theta_max_deg)

    def law(theta):
        theta = np.asarray(theta, dtype=float)
        return b_in + slope * theta, np.full_like(theta, slope)

    return law


# ---------------------------------------------------------------------------
# Closed-form geometry
# ---------------------------------------------------------------------------

def frame(theta, params: SpiralParams = DEFAULT_PARAMS):
    """Base point r, tangent T, normal N (each shape (3, ...)) and speed s."""
    theta = np.asarray(theta, dtype=float)
    a, b = params.a, params.b
    phi = np.arctan2(1.0, b)  # cot(phi) = b
    rho = a * np.exp(b * theta)
    zero = np.zeros_like(theta)
    r = np.stack([rho * np.cos(theta), rho * np.sin(theta), zero])
    T = np.stack([np.cos(theta + phi), np.sin(theta + phi), zero])
    N = np.stack([-np.sin(theta + phi), np.cos(theta + phi), zero])
    s = a * np.sqrt(1 + b * b) * np.exp(b * theta)
    return r, T, N, s


def _Z(like):
    Z = np.zeros_like(like)
    Z[2] = 1.0
    return Z


def surface(theta, w, lean=notebook_lean, params: SpiralParams = DEFAULT_PARAMS):
    """X(theta, w) = r - w tan(beta) N + w Z on broadcastable arrays."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    r, _, N, _ = frame(theta, params)
    beta, _ = lean(theta)
    return r - w * np.tan(beta) * N + w * _Z(r)


def normal(theta, w, lean=notebook_lean, params: SpiralParams = DEFAULT_PARAMS):
    """Un-normalized normal X_theta x X_w = B T - A N - A tan(beta) Z,
    with A = s + w tan(beta) and B = -w beta' / cos^2(beta)."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, T, N, s = frame(theta, params)
    beta, dbeta = lean(theta)
    A = s + w * np.tan(beta)
    B = -w * dbeta / np.cos(beta) ** 2
    return B * T - A * N - A * np.tan(beta) * _Z(T)


def unit_normal(theta, w, lean=notebook_lean, params: SpiralParams = DEFAULT_PARAMS):
    n = normal(theta, w, lean, params)
    return n / np.linalg.norm(n, axis=0)


def gaussian_curvature(theta, w, lean=notebook_lean, params: SpiralParams = DEFAULT_PARAMS):
    """K = -(s beta')^2 / [(s + w tan beta)^2 + (w beta' / cos beta)^2]^2."""
    theta, w = np.broadcast_arrays(np.asarray(theta, float), np.asarray(w, float))
    _, _, _, s = frame(theta, params)
    beta, dbeta = lean(theta)
    denom = (s + w * np.tan(beta)) ** 2 + (w * dbeta / np.cos(beta)) ** 2
    return -(s * dbeta) ** 2 / denom ** 2


def level_curve(height_c: float, theta, lean=notebook_lean,
                params: SpiralParams = DEFAULT_PARAMS):
    """Slice of the blade at z = c: the offset curve r - c tan(beta) N."""
    theta = np.asarray(theta, float)
    return surface(theta, np.full_like(theta, height_c), lean, params)


# ---------------------------------------------------------------------------
# Grid helpers used by the scripts
# ---------------------------------------------------------------------------

def _grid(theta_min_deg, theta_max_deg, params, n_theta, n_w):
    th = np.linspace(np.deg2rad(theta_min_deg), np.deg2rad(theta_max_deg), n_theta)
    wv = np.linspace(0.0, params.height, n_w)
    TH, W = np.meshgrid(th, wv, indexing="ij")
    return th, wv, TH, W


def channel_surface(theta_min_deg: float, theta_max_deg: float,
                    params: SpiralParams = DEFAULT_PARAMS, lean=notebook_lean,
                    n_theta: int = 140, n_w: int = 24):
    """Mid-surface between the planes z = 0 and z = H, sampled on a grid."""
    th, wv, TH, W = _grid(theta_min_deg, theta_max_deg, params, n_theta, n_w)
    X = surface(TH, W, lean, params)
    return X[0], X[1], X[2], th, wv


def channel_shell(theta_min_deg: float, theta_max_deg: float, sign: float,
                  params: SpiralParams = DEFAULT_PARAMS, lean=notebook_lean,
                  n_theta: int = 140, n_w: int = 24, offset: str = "normal"):
    """One wall face, offset by sign * h/2.

    offset="normal"     : along the unit surface normal (the notebook's export)
    offset="horizontal" : inside each horizontal slice, so z stays in [0, H]
    """
    _, _, TH, W = _grid(theta_min_deg, theta_max_deg, params, n_theta, n_w)
    n = normal(TH, W, lean, params)
    if offset == "horizontal":
        n[2] = 0.0
    elif offset != "normal":
        raise ValueError("offset must be 'normal' or 'horizontal'")
    X = surface(TH, W, lean, params) + sign * params.wall_thickness / 2 * n / np.linalg.norm(n, axis=0)
    return X[0], X[1], X[2]


def control_function_surface(theta_max_deg: float, beta_in_deg: float, beta_out_deg: float,
                             params: SpiralParams = DEFAULT_PARAMS,
                             n_theta: int = 140, n_w: int = 24):
    """Surface and Gaussian curvature for a linear lean law beta_in -> beta_out."""
    lean = linear_lean(theta_max_deg, beta_in_deg, beta_out_deg)
    th, wv, TH, W = _grid(0.0, theta_max_deg, params, n_theta, n_w)
    X = surface(TH, W, lean, params)
    K = gaussian_curvature(TH, W, lean, params)
    return X[0], X[1], X[2], K, th, wv


def control_function_values(theta_max_deg: float, beta_in_deg: float, beta_out_deg: float,
                            n_theta: int = 140):
    """(theta_deg, beta_deg) samples of the linear lean law."""
    th = np.linspace(0.0, np.deg2rad(theta_max_deg), n_theta)
    beta, _ = linear_lean(theta_max_deg, beta_in_deg, beta_out_deg)(th)
    return np.rad2deg(th), np.rad2deg(beta)


# ---------------------------------------------------------------------------
# Independent check against a symbolic derivation (needs SymPy)
# ---------------------------------------------------------------------------

def symbolic_check(params: SpiralParams = DEFAULT_PARAMS):
    """Rebuild the notebook construction symbolically and compare with the closed forms."""
    import sympy as sp

    th, v = sp.symbols("theta v", real=True)
    a, b = sp.nsimplify(params.a), sp.nsimplify(params.b)
    r = sp.Matrix([a * sp.exp(b * th) * sp.cos(th), a * sp.exp(b * th) * sp.sin(th), 0])
    rp = r.diff(th)
    s = sp.sqrt(sp.simplify(rp.dot(rp)))
    T = rp / s
    Tp = T.diff(th)
    kappa_s = sp.simplify(sp.sqrt(sp.simplify(Tp.dot(Tp))))
    N = sp.simplify(Tp / kappa_s)
    beta = sp.pi / 4 - th / 8
    F = -sp.sin(beta) * N + sp.cos(beta) * sp.Matrix([0, 0, 1])

    print("kappa * s                 =", kappa_s)
    print("height of r + (v/cos b) F =", sp.simplify((r + v / sp.cos(beta) * F)[2]))
    print("det(r', F, F') / s        =", sp.simplify(sp.Matrix.hstack(rp, F, F.diff(th)).det() / s))

    # Gaussian curvature from the fundamental forms of r + v F (v = w / cos beta)
    X = r + v * F
    Xt, Xv = X.diff(th), X.diff(v)
    n = Xt.cross(Xv)
    K_expr = (Xt.diff(th).dot(n) * Xv.diff(v).dot(n) - Xt.diff(v).dot(n) ** 2) / n.dot(n) ** 2
    K_num = sp.lambdify((th, v), K_expr, "numpy")
    worst = 0.0
    for t0, v0 in [(0.3, 0.1), (2.0, 0.2), (5.0, 0.35)]:
        w0 = v0 * np.cos(np.pi / 4 - t0 / 8)
        worst = max(worst, abs(K_num(t0, v0) - gaussian_curvature(t0, w0, params=params)))
    print("closed-form K, max error  =", f"{worst:.1e}")


if __name__ == "__main__":
    x, y, z, th, wv = channel_surface(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG)
    TH, W = np.meshgrid(th, wv, indexing="ij")
    K = gaussian_curvature(TH, W)
    print("surface grid              :", x.shape)
    print("height range              :", z.min(), "to", z.max())
    print("Gaussian curvature range  :", K.min(), "to", K.max())
    for sign in (+1, -1):
        for mode in ("normal", "horizontal"):
            zs = channel_shell(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG, sign, offset=mode)[2]
            print(f"face {sign:+d}, {mode:10s} offset: z in [{zs.min():.4f}, {zs.max():.4f}]")
    symbolic_check()
