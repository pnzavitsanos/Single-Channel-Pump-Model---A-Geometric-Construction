"""
export_cad.py
==============

Writes the impeller as three closed STL solids (blade, hub disc, shroud disc)
to cad/, ready to import into CAD (Autodesk Inventor, Fusion, FreeCAD, ...).

Usage:
    python scripts/export_cad.py                 # the notebook design
    python scripts/export_cad.py 15 28 45 10     # beta1 beta2 lean1 lean2 (degrees)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from impeller_geometry import NOTEBOOK, impeller_meshes, write_stl, summary

OUT = os.path.join(os.path.dirname(__file__), "..", "cad")


def main(argv):
    d = NOTEBOOK
    if argv:
        b1, b2, l1, l2 = map(float, argv)
        d = d.with_(beta1_deg=b1, beta2_deg=b2, lean1_deg=l1, lean2_deg=l2)
    os.makedirs(OUT, exist_ok=True)
    for name, tris in impeller_meshes(d).items():
        path = os.path.join(OUT, f"impeller_{name}.stl")
        write_stl(tris, path, name)
        print(f"wrote {path} ({len(tris)} triangles)")
    summary(d)


if __name__ == "__main__":
    main(sys.argv[1:])
