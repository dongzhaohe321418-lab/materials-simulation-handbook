#!/usr/bin/env python3
"""Simplified FCC first Brillouin zone (truncated octahedron) with high-symmetry points."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

setup()

W_perms = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        for s3 in [1, -1]:
            for perm in [(1, 0.5, 0), (1, 0, 0.5), (0.5, 1, 0), (0, 1, 0.5),
                         (0.5, 0, 1), (0, 0.5, 1)]:
                W_perms.append((s1*perm[0], s2*perm[1], s3*perm[2]))
W_perms = np.unique(np.round(W_perms, 6), axis=0) * 2*np.pi

hull = ConvexHull(W_perms)

fig = plt.figure(figsize=(7, 6.2))
ax = fig.add_subplot(111, projection='3d')

polys = [W_perms[s] for s in hull.simplices]
ax.add_collection3d(Poly3DCollection(polys, alpha=0.18, facecolor=COLOURS["blue"],
                                     edgecolor=COLOURS["blue"], linewidths=0.9))

points = {
    r"$\Gamma$": np.array([0, 0, 0]),
    "X": np.array([0, 2*np.pi, 0]),
    "L": np.array([np.pi, np.pi, np.pi]),
    "K": np.array([1.5*np.pi, 1.5*np.pi, 0]),
    "U": np.array([2*np.pi, 0.5*np.pi, 0.5*np.pi]),
    "W": np.array([2*np.pi, np.pi, 0]),
}
for label, p in points.items():
    ax.scatter([p[0]], [p[1]], [p[2]], color=COLOURS["red"], s=60, zorder=10)
    ax.text(p[0]*1.10+0.3, p[1]*1.10+0.3, p[2]*1.10+0.5, label,
            fontsize=14, color="k", weight="bold")

path = np.array([points[r"$\Gamma$"], points["X"], points["W"],
                 points["L"], points[r"$\Gamma$"], points["K"]])
ax.plot(path[:,0], path[:,1], path[:,2], color=COLOURS["orange"], lw=2.5,
        label=r"path $\Gamma$-X-W-L-$\Gamma$-K")

ax.set_xlabel(r"$k_x$"); ax.set_ylabel(r"$k_y$"); ax.set_zlabel(r"$k_z$")
ax.set_title("FCC first Brillouin zone (truncated octahedron)\nwith high-symmetry points")
lim = 7.5
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.view_init(elev=22, azim=30)
ax.legend(frameon=False, loc="upper left")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03/fig_brillouin_zone_fcc.png")
save(fig, out)
