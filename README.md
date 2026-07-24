# Single-Channel Pump Skeleton: A Log-Spiral, Twisted-Ribbon Construction

**Author:** Panagiotis Zavitsanos
**Status:** geometric skeleton complete; fluid-dynamic characterization not yet performed (see [Limitations](#8-limitations-and-future-work))

<p align="center">
  <img src="media/hero_channel.png" alt="Single-channel pump skeleton, baseline inlet/outlet angles" width="640">
</p>

## Abstract

This repository documents the derivation and computational implementation of a parametric surface
intended as the geometric skeleton of a single-channel pump: a spiral channel of finite wall
thickness, swept along a logarithmic-spiral base curve with a slowly varying cross-sectional twist.
The construction generalizes the classical Möbius strip's rotating-fiber idea — normally a uniform
half-twist applied around a circle — to a non-uniform twist applied around a logarithmic spiral,
producing an open (non-closed) ribbon whose start angle ("inlet") and end angle ("outlet") are free
design parameters. The derivation proceeded through five iterations, worked out symbolically and
implemented in Mathematica, and is reproduced here in Python together with an animation illustrating
the surface's sensitivity to small perturbations of the inlet and outlet angles.

## Table of contents

1. [Motivation](#1-motivation)
2. [Development progression](#2-development-progression)
3. [Mathematical formulation](#3-mathematical-formulation)
4. [Inlet/outlet angle parametrization](#4-inletoutlet-angle-parametrization)
5. [Animation](#5-animation)
6. [Repository structure](#6-repository-structure)
7. [Reproducing the results](#7-reproducing-the-results)
8. [Limitations and future work](#8-limitations-and-future-work)
9. [Acknowledgments](#9-acknowledgments)

## 1. Motivation

A single-channel pump moves fluid through one continuous spiral passage rather than through several
blade-separated channels, as in a conventional multi-vane impeller. Before any fluid-dynamic design
can begin, that passage needs a well-defined geometric skeleton: a surface that sweeps from an inlet
opening to an outlet opening along a continuously curving path, with a controllable cross-section and
wall thickness. This project develops that skeleton from first principles in differential geometry,
rather than starting from a CAD primitive, so that every design parameter (spiral growth rate, twist
law, wall thickness, inlet angle, outlet angle) is traceable to an explicit formula.

## 2. Development progression

The construction was developed iteratively across five Mathematica notebooks, worked out through a
sequence of exchanges recorded (and lightly deduplicated) in
[`docs/derivation_notes.md`](docs/derivation_notes.md). Each notebook is kept in
[`notebooks/`](notebooks/), numbered in the order below.

| Step | Notebook | Date | Contribution |
|---|---|---|---|
| 1 | [`01_homotopy_example.nb`](notebooks/01_homotopy_example.nb) | 2025-11-19 | Establishes the basic tool used throughout: a homotopy `H(s,t) = (1-t) f0(s) + t f1(s)` between two 3D curves, animated with a parameter slider. This is the warm-up exercise that motivates thinking of a family of curves as a continuously deformable surface. |
| 2 | [`02_mobius_animations.nb`](notebooks/02_mobius_animations.nb) | 2025-11-22 | Builds the classical Möbius strip, `(1 + w cos(θ/2)) (cos θ, sin θ, · ) , w sin(θ/2)`, and animates the rotating cross-sectional fiber to make the half-twist mechanism visible: at each angle θ, the fiber direction is a rotation of the radial vector `R(θ)` and the vertical axis `Z` by `θ/2`. |
| 3 | [`03_log_spiral_v1.nb`](notebooks/03_log_spiral_v1.nb) | 2025-11-24 | Replaces the circle with a logarithmic spiral `r(θ) = a e^{bθ}(cos θ, sin θ)` and derives its Frenet frame (tangent `T(θ)`, normal `N(θ)`) symbolically. The rotating-fiber idea from Step 2 is reapplied using `N(θ)` in place of the circle's radial vector, producing the first spiral ribbon surface `X(θ, w)`. |
| 4 | [`04_modified_mobius_curvature.nb`](notebooks/04_modified_mobius_curvature.nb) | 2025-12-01 | A parallel branch on the *circular* Möbius strip (not the spiral): adds a warping term `α v² sin(u)` to introduce non-zero Gaussian curvature, and confirms the curvature formula analytically. This established that the flat ribbon construction can be deliberately curved — a fact used qualitatively in later design discussion, though the curvature term itself was not carried into the spiral channel. |
| 5 | [`05_log_spiral_v2_final.nb`](notebooks/05_log_spiral_v2_final.nb) | 2026-01-08 | The final iteration and the basis of this repository. Revises Step 3's twist law to a slower, partial rotation `β(θ) = π/4 − θ/8` (rather than a full `θ/2` half-twist), rescales the offset so the channel width stays consistent as the twist angle varies, and — new in this step — gives the ribbon finite wall thickness by offsetting along its own surface normal, then samples and exports the resulting point cloud to CSV for CAD (Autodesk Inventor) import. [`data/spiral_3D_points.xlsx`](data/spiral_3D_points.xlsx) is that exported point cloud. |

Steps 1–3 and 5 form a single connected line of development (homotopy → rotating fiber → spiral
generalization → finite-thickness channel with tunable sweep). Step 4 is a side investigation into
intrinsic curvature that informed the design conceptually but is architecturally independent of the
final channel surface.

## 3. Mathematical formulation

All symbols below are implemented symbolically in [`src/channel_model.py`](src/channel_model.py),
mirroring the derivation in Step 5.

**Base curve.** A logarithmic spiral in the plane,

```
r(θ) = a e^{bθ} (cos θ, sin θ),        a = 0.2, b = 0.2
```

**Frenet frame.** The unit tangent and (principal) unit normal, obtained by differentiation and
normalization:

```
T(θ) = r'(θ) / |r'(θ)|
N(θ) = T'(θ) / |T'(θ)|
```

**Twist law.** A partial-twist angle, much slower than the classical Möbius strip's `θ/2`:

```
β(θ) = π/4 − θ/8
```

**Fiber direction.** A unit vector blending the in-plane normal (lifted to 3D with a zero
z-component, written `N̂`) with the vertical axis `Ẑ = (0,0,1)`:

```
F(θ) = −sin(β(θ)) N̂(θ) + cos(β(θ)) Ẑ
```

**Channel mid-surface.** The spiral, lifted to 3D (`r̂`, zero z-component), swept along `F(θ)` and
offset by the cross-channel coordinate `w`:

```
X(θ, w) = r̂(θ) + (w / cos β(θ)) F(θ)
```

**Wall thickness.** The surface normal `n(θ, w) = (X_θ × X_w)/|X_θ × X_w|` gives two finite-thickness
shells at offset `t = ±h/2`:

```
X_thick(θ, w, t) = X(θ, w) + t · n(θ, w),      h = 0.05
```

The channel cross-section spans `w ∈ [−0.1, 0.1]` in the visualizations in this repository (the
original CAD export in Step 5 instead used the one-sided range `w ∈ [0, 0.3]`, which produces the
same surface shape, just parametrized from one edge rather than the centerline).

## 4. Inlet/outlet angle parametrization

`θ` sweeps the channel from `θ_inlet` to `θ_outlet`. In Step 5's final render these were `0°` and
`342°` (`1.9π` rad) respectively — just short of a full turn, leaving a visible gap between the
channel's two open ends. Because the spiral and twist law are defined for all `θ`, these two bounds
can be moved independently to widen, narrow, rotate, or reposition the channel opening without
touching the underlying surface formula. This is the one design freedom explored quantitatively in
this repository (see below); wall thickness, spiral rate, and twist law are held fixed at their
Step 5 values.

## 5. Animation

[`scripts/make_animation.py`](scripts/make_animation.py) renders
[`media/inlet_outlet_sweep.gif`](media/inlet_outlet_sweep.gif), which perturbs the inlet and outlet
angles by small, independent, out-of-phase sinusoids around the baseline (±12° on the inlet, ±15° on
the outlet) and recomputes the full surface at each frame. The surface is colored by θ (a proxy for
position along the flow path) to make the sweep direction legible. The on-screen title reports the
exact inlet/outlet angle for the current frame.

<p align="center">
  <img src="media/inlet_outlet_sweep.gif" alt="Inlet/outlet angle sensitivity animation" width="480">
</p>

[`scripts/render_hero.py`](scripts/render_hero.py) renders the static, finite-thickness banner image
at the top of this document, at the baseline inlet/outlet angles.

## 6. Repository structure

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── notebooks/            original Mathematica notebooks, in development order
├── docs/
│   └── derivation_notes.md   deduplicated transcript of the symbolic derivation
├── src/
│   └── channel_model.py      SymPy/NumPy re-implementation of the Step 5 geometry
├── scripts/
│   ├── render_hero.py        static banner render (finite thickness)
│   └── make_animation.py     inlet/outlet angle-sensitivity animation
├── data/
│   └── spiral_3D_points.xlsx original CAD-export point cloud from Step 5
└── media/
    ├── hero_channel.png
    └── inlet_outlet_sweep.gif
```

## 7. Reproducing the results

```bash
python -m venv .venv && source .venv/bin/activate   # optional
pip install -r requirements.txt

python src/channel_model.py        # sanity-check the geometry (prints grid ranges)
python scripts/render_hero.py      # -> media/hero_channel.png
python scripts/make_animation.py   # -> media/inlet_outlet_sweep.gif
```

To explore a different inlet/outlet angle by hand:

```python
from src.channel_model import channel_surface

x, y, z, theta, w = channel_surface(theta_min_deg=-10, theta_max_deg=330)
```

## 8. Limitations and future work

This repository establishes the **geometric skeleton only**. No fluid-mechanical modeling — flow
rate, head, velocity/pressure field, Reynolds number, cavitation margin, or any CFD analysis — has
been performed on this surface. The wall-thickness offset and CAD export in Step 5 are the only steps
taken toward a manufacturable part; the "pump" framing describes the intended end use of the
geometry, not a validated pump design.

Natural next steps:

- Fluid-dynamic simulation (e.g. CFD) over the channel interior to characterize flow behavior.
- A structured sweep over spiral parameters `(a, b)` and twist law, not just inlet/outlet angle.
- Solid-model export (STEP/STL) directly from `channel_model.py` rather than via the CSV point-cloud
  round trip used in Step 5.
- Reconciling the two `w` conventions (centered `[-0.1, 0.1]` vs. one-sided `[0, 0.3]`) into a single
  parametrization.

## 9. Acknowledgments

The symbolic derivation underlying this geometry was worked out in conversation with ChatGPT; the
deduplicated transcript is preserved in [`docs/derivation_notes.md`](docs/derivation_notes.md) for
provenance. All Mathematica code was written and run by the author; the Python reimplementation in
this repository was produced with the assistance of Claude.

## License

MIT — see [`LICENSE`](LICENSE).
