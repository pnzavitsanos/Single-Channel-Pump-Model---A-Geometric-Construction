"""
export_cad.py
==============

Writes the blade and the two plates as closed STL solids to cad/, ready to
import into CAD (Autodesk Inventor, Fusion, FreeCAD, ...).

Usage:
    python scripts/export_cad.py                              # the notebook design
    python scripts/export_cad.py --spiral archimedean --b 0.0771
    python scripts/export_cad.py --spiral "0.2 + 0.05*theta + 0.01*theta**2"
    python scripts/export_cad.py --inlet 30 --outlet 10 --height 0.25
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from blade_geometry import (NOTEBOOK, LogSpiral, ArchimedeanSpiral, CustomSpiral,
                            pump_meshes, write_stl, summary)

OUT = os.path.join(os.path.dirname(__file__), "..", "cad")


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spiral", default="log",
                   help='"log", "archimedean", or a formula in theta, e.g. "0.2*exp(0.15*theta)"')
    p.add_argument("--a", type=float, default=0.2, help="start radius (log and Archimedean)")
    p.add_argument("--b", type=float, default=None,
                   help="opening rate (default 0.2 for log, 0.0771 for Archimedean)")
    p.add_argument("--extent", type=float, default=NOTEBOOK.wrap_deg, help="angular extent Phi, degrees")
    p.add_argument("--inlet", type=float, default=NOTEBOOK.inlet_angle_deg, help="inlet angle, degrees")
    p.add_argument("--outlet", type=float, default=NOTEBOOK.outlet_angle_deg, help="outlet angle, degrees")
    p.add_argument("--height", type=float, default=NOTEBOOK.height, help="distance H between the planes")
    p.add_argument("--thickness", type=float, default=NOTEBOOK.thickness, help="wall thickness h")
    return p.parse_args(argv)


def main(argv):
    a = parse(argv)
    if a.spiral == "log":
        spiral = LogSpiral(a.a, 0.2 if a.b is None else a.b)
    elif a.spiral == "archimedean":
        spiral = ArchimedeanSpiral(a.a, 0.0771 if a.b is None else a.b)
    else:
        spiral = CustomSpiral(a.spiral)
    b = NOTEBOOK.with_(spiral=spiral, wrap_deg=a.extent, inlet_angle_deg=a.inlet,
                       outlet_angle_deg=a.outlet, height=a.height, thickness=a.thickness)
    os.makedirs(OUT, exist_ok=True)
    for name, tris in pump_meshes(b).items():
        path = os.path.join(OUT, f"{name}.stl")
        write_stl(tris, path, name)
        print(f"wrote {path} ({len(tris)} triangles)")
    summary(b)


if __name__ == "__main__":
    main(sys.argv[1:])
