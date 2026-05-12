"""Common style settings for handbook figures (quality-upgrade compliant)."""
import os
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

# Wong colour-blind-safe palette
COLOURS = {
    "blue":   "#0173b2",
    "orange": "#de8f05",
    "green":  "#029e73",
    "purple": "#cc78bc",
    "tan":    "#ca9161",
    "red":    "#d55e00",
    "yellow": "#ece133",
    "skyblue":"#56b4e9",
    "grey":   "#555555",
}
PALETTE = [COLOURS["blue"], COLOURS["orange"], COLOURS["green"],
           COLOURS["purple"], COLOURS["tan"], COLOURS["red"]]


def setup():
    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "lines.linewidth": 2,
        "axes.linewidth": 1.2,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 100,
        "savefig.dpi": 300,
    })


def save(fig, path):
    """Save a figure with the standard quality settings and verify it."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    plt.close(fig)
    # Verify
    img = Image.open(path)
    w, h = img.size
    size = os.path.getsize(path)
    print(f"wrote {path}  ({w}x{h}, {size/1024:.1f} KB)")
    if w < 1500:
        print(f"  WARNING: width {w} < 1500")
    if size > 500_000:
        print(f"  WARNING: file size {size} > 500 kB")
    return w, h, size
