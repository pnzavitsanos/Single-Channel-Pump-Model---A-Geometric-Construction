"""
channel_model.py
=================

Python re-implementation of the log-spiral channel geometry developed in
``notebooks/05_log_spiral_v2_final.nb`` ("Log try 2"), the final iteration of
the modeling progression documented in this repository (see ``README.md``).

The construction, in order:

1.  A logarithmic spiral base curve ``r(theta) = a * exp(b*theta) * (cos theta, sin theta)``.
2.  Its Frenet data: speed, unit tangent ``T(theta)``, and unit normal ``N(theta)``
    (obtained symbolically, exactly as in the notebook, to avoid any numerical
    drift in the derivative chain).
3.  A slowly-varying "fiber" direction ``F(theta)`` that blends the in-plane
    normal ``N(theta)`` with the vertical axis ``Z`` through an angle
    ``beta(theta) = pi/4 - theta/8``. This plays the role of the classic
    Mobius strip's half-twist, but at a much slower, partial rate.
4.  A ribbon/channel surface ``X(theta, w)`` that sweeps ``F(theta)`` along the
    spiral, offset by the channel's half-width coordinate ``w``.
5.  A finite wall thickness obtained by offsetting ``X`` along its own surface
    normal by +/- h/2, giving two shells (``X_thick``) that can be rendered as
    a solid channel.

Everything is derived symbolically once with SymPy (mirroring the original
Mathematica derivation) and then compiled to fast NumPy functions with
``lambdify`` so it can be evaluated on a dense grid for plotting/animation.

Inlet / outlet angle
---------------------
``theta`` sweeps the spiral from ``theta_min`` ("inlet") to ``theta_max``
("outlet"). These two angles are exposed as the primary design parameters of
the channel: changing them changes how far around the spiral the channel
opening starts and ends, without altering the underlying spiral or twist law.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# 1-2. Symbolic spiral + Frenet frame (matches notebooks/05_log_spiral_v2_final.nb)
# ---------------------------------------------------------------------------

theta, w, t, a, b = sp.symbols("theta w t a b", real=True)

# Logarithmic spiral, a=0.2, b=0.2 in the source notebook.
r_expr = sp.Matrix([a * sp.exp(b * theta) * sp.cos(theta),
                     a * sp.exp(b * theta) * sp.sin(theta)])

rp_expr = r_expr.diff(theta)
speed_expr = sp.sqrt((rp_expr.T * rp_expr)[0])
T_expr = sp.simplify(rp_expr / speed_expr)
Tp_expr = T_expr.diff(theta)
Tp_norm = sp.sqrt(sp.simplify((Tp_expr.T * Tp_expr)[0]))
N_expr = sp.simplify(Tp_expr / Tp_norm)  # "NormalVec" in the notebook (N is reserved there)

# ---------------------------------------------------------------------------
# 3-4. Fiber direction and ribbon surface
# ---------------------------------------------------------------------------

Z_expr = sp.Matrix([0, 0, 1])
beta_expr = sp.pi / 4 - theta / 8  # twist angle law from the notebook

N3_expr = sp.Matrix([N_expr[0], N_expr[1], 0])
F_expr = -sp.sin(beta_expr) * N3_expr + sp.cos(beta_expr) * Z_expr

r3_expr = sp.Matrix([r_expr[0], r_expr[1], 0])
X_expr = r3_expr + (w / sp.cos(beta_expr)) * F_expr

# ---------------------------------------------------------------------------
# 5. Wall thickness via the surface normal
# ---------------------------------------------------------------------------

Xu_expr = X_expr.diff(theta)
Xw_expr = X_expr.diff(w)
cross_expr = Xu_expr.cross(Xw_expr)
cross_norm = sp.sqrt(sp.simplify((cross_expr.T * cross_expr)[0]))
SurfaceNormal_expr = cross_expr / cross_norm

Xthick_expr = X_expr + t * SurfaceNormal_expr

# ---------------------------------------------------------------------------
# Compile to fast numeric functions
# ---------------------------------------------------------------------------

_X_num = sp.lambdify((theta, w, a, b), X_expr, modules="numpy")
_Xthick_num = sp.lambdify((theta, w, t, a, b), Xthick_expr, modules="numpy")
_SurfaceNormal_num = sp.lambdify((theta, w, a, b), SurfaceNormal_expr, modules="numpy")


@dataclass(frozen=True)
class SpiralParams:
    """Physical/geometric parameters of the log-spiral channel."""

    a: float = 0.2          # spiral scale
    b: float = 0.2          # spiral growth rate
    w_min: float = -0.1     # channel cross-section extent (inner edge)
    w_max: float = 0.1      # channel cross-section extent (outer edge)
    wall_thickness: float = 0.05  # h in the notebook


DEFAULT_PARAMS = SpiralParams()

# Baseline inlet/outlet sweep used throughout the notebook's final renders.
BASELINE_INLET_DEG = 0.0
BASELINE_OUTLET_DEG = 1.9 * 180.0  # 1.9*pi rad, expressed in degrees


def channel_surface(theta_min_deg: float, theta_max_deg: float,
                     params: SpiralParams = DEFAULT_PARAMS,
                     n_theta: int = 140, n_w: int = 24):
    """Evaluate the mid-surface X(theta, w) on a grid.

    Parameters
    ----------
    theta_min_deg, theta_max_deg:
        Inlet / outlet angle of the channel sweep, in degrees.
    """
    th = np.linspace(np.deg2rad(theta_min_deg), np.deg2rad(theta_max_deg), n_theta)
    wv = np.linspace(params.w_min, params.w_max, n_w)
    TH, W = np.meshgrid(th, wv, indexing="ij")
    X = np.asarray(_X_num(TH, W, params.a, params.b), dtype=float)
    # lambdify on a Matrix broadcasts to shape (3, 1, n_theta, n_w); squeeze it.
    X = X.reshape(3, n_theta, n_w)
    return X[0], X[1], X[2], th, wv


def channel_shell(theta_min_deg: float, theta_max_deg: float,
                   sign: float, params: SpiralParams = DEFAULT_PARAMS,
                   n_theta: int = 140, n_w: int = 24):
    """Evaluate one of the two finite-thickness shells (sign = +1 or -1)."""
    th = np.linspace(np.deg2rad(theta_min_deg), np.deg2rad(theta_max_deg), n_theta)
    wv = np.linspace(params.w_min, params.w_max, n_w)
    TH, W = np.meshgrid(th, wv, indexing="ij")
    t_val = sign * params.wall_thickness / 2.0
    Xt = np.asarray(_Xthick_num(TH, W, t_val, params.a, params.b), dtype=float)
    Xt = Xt.reshape(3, n_theta, n_w)
    return Xt[0], Xt[1], Xt[2]


if __name__ == "__main__":
    x, y, z, th, wv = channel_surface(BASELINE_INLET_DEG, BASELINE_OUTLET_DEG)
    print("Surface grid:", x.shape)
    print("theta range (deg):", np.rad2deg(th[0]), "to", np.rad2deg(th[-1]))
    print("x range:", x.min(), x.max())
    print("y range:", y.min(), y.max())
    print("z range:", z.min(), z.max())
