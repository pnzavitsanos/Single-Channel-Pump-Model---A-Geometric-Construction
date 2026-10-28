# Single-Channel Pump Blade: a Ruled Surface Between Two Parallel Planes

**Author:** Panagiotis Zavitsanos  
**Status:** geometry complete and verified; no fluid-dynamic analysis yet (see [Limitations](#9-limitations-and-future-work))

<p align="center">
  <img src="media/hero_channel.png" alt="Log-spiral blade with wall thickness between two parallel planes" width="640">
</p>

## Abstract

The blade of a single-channel pump has to wind around the axis and sit between two parallel discs (the hub and the shroud). This project builds such a blade from differential geometry rather than from a CAD primitive.

The surface is **ruled**: it is a one-parameter family of straight segments (a foliation of the surface by lines) attached along a logarithmic spiral. Each segment lies in the normal plane of the spiral and is turned about the spiral's tangent by a lean angle $\beta(\theta)$. This is the same rotating-fibre mechanism as the Möbius strip, transported from a circle to a spiral. Each segment is then stretched by $1/\cos\beta$, which makes every one of them start on the plane $z=0$ and end exactly on the plane $z=H$, whatever its lean. The blade therefore fits between two parallel planes by construction.

Because the construction is explicit, its geometry can be worked out in closed form:

- the spiral's Frenet frame is explicit, with $\kappa s = 1$;
- the height of every point equals its parameter, $z = w$;
- each horizontal slice of the blade is an offset curve of the spiral;
- the Gaussian curvature has a closed form, and the blade is developable (it can be rolled from flat sheet) **if and only if the lean angle is constant**.

Everything is implemented in [`src/channel_model.py`](src/channel_model.py) and checked against an independent symbolic derivation.

## Where the idea came from

The starting point was a drawing from Frederic Schuller's lecture course *Geometric Anatomy of Theoretical Physics*, in the part on connections and parallel transport. It shows a curve in a base manifold, a fibre standing over each point of the curve, and the curve lifted up through those fibres to a new height, where the lifted points are joined into a new curve (the horizontal lift).

This project asks what that picture looks like when everything is concrete and Euclidean:

- the base is the plane $z = 0$, and the curve is a logarithmic spiral;
- the fibre over each point is a real straight segment, a ruling of the blade;
- the lifted curve is the top edge of the blade, at height $z = H$.

In a bundle, a connection decides how a point is carried from fibre to fibre along the curve. Here that choice is made by hand: each ruling is turned about the spiral's tangent by a lean angle $\beta(\theta)$, the way the fibre of a Möbius strip turns as it goes around the circle. The two-plane property of [§3](#3-the-two-plane-property) is what makes "lifted to a new height and joined there" exact.

## Table of contents

1. [The idea in one animation](#1-the-idea-in-one-animation)
2. [Construction](#2-construction)
3. [The two-plane property](#3-the-two-plane-property)
4. [Two foliations: rulings and horizontal slices](#4-two-foliations-rulings-and-horizontal-slices)
5. [Curvature and developability](#5-curvature-and-developability)
6. [Wall thickness](#6-wall-thickness)
7. [The lean law as a design control](#7-the-lean-law-as-a-design-control)
8. [Development history and reproduction](#8-development-history-and-reproduction)
9. [Limitations and future work](#9-limitations-and-future-work)

## 1. The idea in one animation

A single straight segment travels along the spiral. At every point it stands in the plane spanned by the normal $N$ and the vertical $Z$, leans away from the vertical by $\beta(\theta)$, and always reaches from the bottom plane to the top plane. The surface it sweeps out is the blade.

<p align="center">
  <img src="media/ruling_sweep.gif" alt="A ruling swept along the logarithmic spiral" width="520">
</p>

## 2. Construction

### 2.1 Base curve: the logarithmic spiral

```math
r(\theta) = a\,e^{b\theta}\,(\cos\theta,\ \sin\theta,\ 0), \qquad a = 0.2,\quad b = 0.2 .
```

The logarithmic spiral is the curve whose tangent makes a constant angle with the radius vector, which is why it is the classical shape of a pump blade with constant blade angle. Writing $\cot\varphi = b$, the Frenet data come out in closed form:

```math
\begin{aligned}
s(\theta) &= |r'(\theta)| = a\sqrt{1+b^2}\,e^{b\theta},\\
T(\theta) &= \bigl(\cos(\theta+\varphi),\ \sin(\theta+\varphi),\ 0\bigr),\\
N(\theta) &= \bigl(-\sin(\theta+\varphi),\ \cos(\theta+\varphi),\ 0\bigr),\\
\kappa(\theta)\,s(\theta) &= 1 .
\end{aligned}
```

So the tangent turns at exactly unit rate in $\theta$, and the radius of curvature $1/\kappa$ equals the speed $s(\theta)$. With $b = 0.2$ the blade angle, measured from the circumferential direction, is $\arctan b \approx 11.3^\circ$ along the whole blade. $(T, N, Z)$ is a right-handed orthonormal frame, with $T \times N = Z$ and $N \times Z = T$.

The normal is computed as the derivative of the unit tangent, $N = T'/|T'|$, not as $r''$. The second derivative depends on how fast the curve is parametrised, while $T'$ only records turning. Because $T \cdot T = 1$, differentiating gives $T \cdot T' = 0$, so $T'$ is automatically perpendicular to the curve.

### 2.2 Rulings: rotation about the tangent

At each point of the spiral a unit direction is chosen in the normal plane $\mathrm{span}\lbrace N, Z\rbrace$:

```math
F(\theta) = -\sin\beta(\theta)\,N(\theta) + \cos\beta(\theta)\,Z .
```

This is the vertical $Z$ rotated about the tangent $T$ by the angle $\beta$:

```math
R_T(\beta)\,Z = \cos\beta\,Z + \sin\beta\,(T \times Z) = \cos\beta\,Z - \sin\beta\,N = F .
```

It is exactly the Möbius-strip mechanism. For the Möbius strip the base curve is a circle and the fibre $F = \cos(\theta/2)\,R + \sin(\theta/2)\,Z$ turns about the tangent at rate $1/2$, half a turn in one loop. Here the base curve is a spiral and the lean follows the law of the final notebook:

```math
\beta(\theta) = \frac{\pi}{4} - \frac{\theta}{8}, \qquad \theta \in [0,\ 1.9\pi],
```

so the blade leans outward by $45^\circ$ at the inlet and straightens to about $2.25^\circ$ at the outlet.

**The rulings only turn.** Since $F \cdot F = 1$, differentiating gives $F \cdot F' = 0$. As the foot of a ruling moves along the spiral, its direction rotates but never stretches, by the same argument that shows the Frenet vectors only rotate. The only change of length in the whole construction is the deliberate factor $1/\cos\beta$ of §2.3.

### 2.3 The blade surface

```math
X(\theta, w) = r(\theta) + \frac{w}{\cos\beta(\theta)}\,F(\theta), \qquad 0 \le w \le H = 0.3 .
```

Substituting $F$ gives the explicit form used in the code:

```math
X(\theta, w) = r(\theta) - w\,\tan\beta(\theta)\,N(\theta) + w\,Z .
```

## 3. The two-plane property

> **Proposition.** For every $\theta$ and $w$, the height of $X(\theta, w)$ is exactly $w$. Each ruling ($\theta$ fixed) is a straight segment of length $H/\cos\beta$ that starts on the plane $z=0$ and ends on the plane $z=H$.

*Proof.* $r$ and $N$ are horizontal, so the only vertical contribution is $\frac{w}{\cos\beta}\cos\beta = w$. $\square$

This is what the factor $1/\cos\beta$ is for. Without it, a ruling of fixed length $\ell$ would reach only the height $\ell\cos\beta$, so the top edge of the blade would sag wherever it leans more and could not meet a flat shroud. With it, the blade fits between two parallel planes for any lean law.

The same formula sets the limit of the construction. A ruling can reach $z=H$ only while $\cos\beta > 0$, so every lean law must satisfy

```math
-90^\circ < \beta(\theta) < 90^\circ .
```

## 4. Two foliations: rulings and horizontal slices

The blade is covered by two transverse families of curves:

- **Rulings** ($\theta$ fixed) are straight segments from plane to plane.
- **Horizontal slices** ($w = c$ fixed) are the intersections of the blade with the planes $z = c$:

```math
X(\theta, c) = r(\theta) - c\,\tan\beta(\theta)\,N(\theta) + c\,Z .
```

Every slice is therefore a **variable-offset curve** of the base spiral: it is the spiral pushed sideways by the distance $c\tan\beta(\theta)$. If the lean is constant the slices are ordinary parallel curves of the spiral.

<p align="center">
  <img src="media/geometry.png" alt="Rulings, horizontal slices and Gaussian curvature of the blade" width="100%">
</p>

**Regularity.** The tangent vectors are

```math
\begin{aligned}
X_\theta &= \bigl(s + w\tan\beta\bigr)\,T \;-\; \frac{w\,\beta'}{\cos^2\beta}\,N,\\
X_w &= -\tan\beta\,N + Z .
\end{aligned}
```

$N$ points toward the centre of curvature. When the blade leans outward ($\beta > 0$, as in the notebook), the slices move away from the centre of curvature and stay smooth. When it leans inward, the offset distance must stay below the radius of curvature $s(\theta)$:

```math
s(\theta) + w\tan\beta(\theta) > 0 .
```

Otherwise the slices pass the evolute of the spiral (itself a logarithmic spiral) and fold back on themselves.

## 5. Curvature and developability

The unnormalised normal is

```math
X_\theta \times X_w = -\frac{w\,\beta'}{\cos^2\beta}\,T \;-\; \bigl(s + w\tan\beta\bigr)\,N \;-\; \bigl(s + w\tan\beta\bigr)\tan\beta\,Z .
```

For a ruled surface $X = r + vF$, the Gaussian curvature is $K = -\det(r', F, F')^2 / |X_\theta \times X_v|^4$. Here

```math
\det\bigl(r',\, F,\, F'\bigr) = s\,\beta', \qquad |F'|^2 = \sin^2\beta + \beta'^2 ,
```

and changing back to the height coordinate $w = v\cos\beta$ gives the closed form

```math
\boxed{\,K(\theta, w) = -\,\frac{\bigl(s\,\beta'\bigr)^2}{\Bigl[\bigl(s + w\tan\beta\bigr)^2 + \bigl(w\,\beta'/\cos\beta\bigr)^2\Bigr]^2}\,}
```

Three consequences follow.

1. **$K \le 0$ everywhere.** The blade is saddle-shaped wherever it is curved at all.
2. **Developability.** $K \equiv 0$ if and only if $\beta' \equiv 0$. A constant lean gives a surface of constant slope over the spiral, which can be rolled from flat sheet metal. Any change of lean along the blade makes it doubly curved, so it has to be pressed or machined.
3. **Along the bottom edge** ($w = 0$) the curvature is minus the square of the twist rate per unit length, where $\sigma$ is arc length:

```math
K\big|_{w=0} = -\frac{\beta'^2}{s^2} = -\left(\frac{d\beta}{d\sigma}\right)^2 .
```

For the notebook's law ($\beta' = -1/8$), $K$ ranges from $-0.376$ at the inlet foot to $-0.0099$ at the top of the outlet. The formula agrees with a symbolic computation from the fundamental forms to $10^{-16}$; run `python src/channel_model.py`.

## 6. Wall thickness

The notebook gives the blade a thickness $h = 0.05$ by offsetting along the unit normal, $X \pm \tfrac{h}{2}\,\hat n$, and exports that point cloud for Autodesk Inventor ([`data/spiral_3D_points.xlsx`](data/spiral_3D_points.xlsx)). Because $\hat n$ has a vertical component wherever $\beta \neq 0$, this offset pushes the edges off the two planes: the faces reach $z \in [-0.018,\ 0.318]$ instead of $[0,\ 0.3]$, most at the $45^\circ$ inlet.

`channel_shell(..., offset="horizontal")` offsets each slice inside its own plane instead, along the in-plane normal of the slice. Both faces then stay exactly between $z=0$ and $z=H$, and the thickness is measured horizontally. The hero image above uses this option.

## 7. The lean law as a design control

The lean law is the free design function of the construction. [`src/channel_model.py`](src/channel_model.py) generalises the notebook's law to a linear law between an inlet lean and an outlet lean:

```math
\beta(\theta) = \beta_{\text{in}} + \bigl(\beta_{\text{out}} - \beta_{\text{in}}\bigr)\,\frac{\theta}{\theta_{\max}}, \qquad -90^\circ < \beta_{\text{in}},\ \beta_{\text{out}} < 90^\circ .
```

The $1/\cos\beta$ stretch is kept, so every member of the family fits between the two planes. Values outside $(-90^\circ, 90^\circ)$ raise an error, because a horizontal ruling can never reach the top plane.

The animation fixes $\beta_{\text{in}} = 45^\circ$ and swings $\beta_{\text{out}}$ from $45^\circ$ to $-45^\circ$ and back. At $\beta_{\text{out}} = 45^\circ$ the lean is constant and the blade is developable ($K = 0$). As soon as the lean varies, the blade becomes saddle-shaped, most strongly at the inlet foot, where $|K| = \beta'^2/s^2$ and $s$ is smallest.

<p align="center">
  <img src="media/control_function_sweep.gif" alt="Lean law sweep and its effect on Gaussian curvature" width="760">
</p>

## 8. Development history and reproduction

The construction was developed in Mathematica; the notebooks are kept in [`notebooks/`](notebooks/).

| Step | Notebook | Date | Contribution |
|---|---|---|---|
| 1 | [`01_homotopy_example.nb`](notebooks/01_homotopy_example.nb) | 2025-11-19 | Homotopy $H(s,t) = (1-t)f_0(s) + t f_1(s)$ between two space curves: a family of curves seen as a surface. |
| 2 | [`02_mobius_animations.nb`](notebooks/02_mobius_animations.nb) | 2025-11-22 | Möbius strip as a circle carrying a rotating fibre $F = \cos(\theta/2)R + \sin(\theta/2)Z$, animated to show the rotation. |
| 3 | [`03_log_spiral_v1.nb`](notebooks/03_log_spiral_v1.nb) | 2025-11-24 | Logarithmic spiral and its Frenet frame derived symbolically. First spiral ribbon, with the fibre $N + Z$ rotated about the vertical axis. |
| 4 | [`04_modified_mobius_curvature.nb`](notebooks/04_modified_mobius_curvature.nb) | 2025-12-01 | Side study: a warping term $\alpha v^2 \sin u$ added to the Möbius strip to give it intrinsic curvature. |
| 5 | [`05_log_spiral_v2_final.nb`](notebooks/05_log_spiral_v2_final.nb) | 2026-01-08 | Final construction: the fibre rotated inside the normal plane (about $T$), lean law $\pi/4 - \theta/8$, the $1/\cos\beta$ stretch so rulings span $0 \le z \le 0.3$, wall thickness $h = 0.05$, and point-cloud export for Inventor. |

The decisive change between steps 3 and 5 was the axis of rotation. Turning the fibre about the vertical axis moves it out of the normal plane of the spiral. Turning it about the tangent keeps it in the normal plane, which is what the Möbius strip does and what makes the two-plane property possible.

In the notebooks the normal is called `NormalVec`, because `N` is a protected built-in function in Mathematica, and planar vectors are lifted to 3D with `Append[..., 0]` before they are combined with $Z$.

The Python code in this repository adds the closed forms, the curvature and developability analysis, the general lean law, and the horizontal wall offset.

```
.
├── src/channel_model.py                    closed-form geometry, curvature, symbolic self-check
├── scripts/
│   ├── render_hero.py                      media/hero_channel.png
│   ├── render_geometry.py                  media/geometry.png
│   ├── make_ruling_animation.py            media/ruling_sweep.gif
│   ├── make_control_function_animation.py  media/control_function_sweep.gif
│   └── plot_utils.py                       shared styling
├── notebooks/                              original Mathematica notebooks
├── data/spiral_3D_points.xlsx              CAD point cloud exported from notebook 5
└── media/
```

```bash
pip install -r requirements.txt

python src/channel_model.py                        # closed forms vs symbolic derivation
python scripts/render_hero.py
python scripts/render_geometry.py
python scripts/make_ruling_animation.py
python scripts/make_control_function_animation.py
```

```python
import sys; sys.path.insert(0, "src")
from channel_model import surface, gaussian_curvature, linear_lean

lean = linear_lean(theta_max_deg=342, beta_in_deg=45, beta_out_deg=10)
X = surface(1.0, 0.15, lean)              # point at theta = 1 rad, height 0.15
K = gaussian_curvature(1.0, 0.15, lean)
```

## 9. Limitations and future work

This repository establishes the geometry only. No flow analysis (head, flow rate, velocity and pressure fields, cavitation) has been carried out, and "pump" describes the intended use rather than a validated design.

Natural next steps:

- Choose the lean law by optimisation, for example to meet given inlet and outlet angles while keeping $\int |K|\,dA$ small, which makes the blade easier to manufacture.
- Make the spiral parameters $(a, b)$ and the plane spacing $H$ part of the design study.
- Export a closed solid (STEP or STL) directly from the closed forms instead of the CSV point cloud.
- Run CFD on the channel between the blade and the two planes.

## License

Public domain, under [the Unlicense](LICENSE). Anyone may copy, modify, use and distribute this work for any purpose, without asking permission.
