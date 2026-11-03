"""
figure1_assemble.py

Stitches figure1a, figure1b, figure1c, figure1d_mets_alluvial side by side
into results/figures/figure1_combined.png.

Run AFTER all four panel scripts have been executed.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

BASE_DIR = "/home/lorenzo/Microbesomics"
PANELS   = [
    os.path.join(BASE_DIR, "results/figures/figure1a.png"),
    os.path.join(BASE_DIR, "results/figures/figure1b.png"),
    os.path.join(BASE_DIR, "results/figures/figure1c.png"),
    os.path.join(BASE_DIR, "results/figures/figure1d_mets_alluvial.png"),
]
OUT = os.path.join(BASE_DIR, "results/figures/figure1_combined.png")

imgs = [mpimg.imread(p) for p in PANELS]

# pad shorter panels with white rows so all heights match (centres vertically)
max_h = max(img.shape[0] for img in imgs)
def pad_to_height(img, target_h):
    """Pad equally top and bottom so panel is centred vertically."""
    h, w = img.shape[:2]
    if h >= target_h:
        return img
    top    = (target_h - h) // 2
    bottom = target_h - h - top
    pad    = np.ones((1, w, img.shape[2]), dtype=img.dtype)
    return np.vstack([np.repeat(pad, top, axis=0),
                      img,
                      np.repeat(pad, bottom, axis=0)])

imgs = [pad_to_height(img, max_h) for img in imgs]

heights  = [img.shape[0] for img in imgs]
widths   = [img.shape[1] for img in imgs]

TARGET_H = max(heights) / 300.0 * 0.82
w_ratios = [w / max(heights) * TARGET_H for w in widths]
total_w  = sum(w_ratios)

fig, axes = plt.subplots(
    1, 4,
    figsize=(total_w, TARGET_H),
    gridspec_kw={"width_ratios": w_ratios},
)
fig.patch.set_facecolor("white")

labels = ["a", "b", "c", "d"]
for ax, img, label in zip(axes, imgs, labels):
    ax.imshow(img, aspect="auto", interpolation="lanczos")
    ax.axis("off")
    ax.text(0.01, 0.99, label, transform=ax.transAxes,
            fontsize=4, fontweight="bold", va="top", ha="left",
            color="#111")

plt.subplots_adjust(wspace=0.03, left=0, right=1, top=1, bottom=0)
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
