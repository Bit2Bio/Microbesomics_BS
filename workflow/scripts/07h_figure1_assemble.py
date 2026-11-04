"""
Assembles all panels into a single combined figure:
  Row 1 (a–d): ridge plots (clinical) + MetS alluvial
  Row 2 (e–h): PCoA Bray-Curtis, BC distances, PCoA Unweighted UniFrac, UUF distances
  Row 3 (i–j): Phylum composition stacked barplot, CLR diverging barplot
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

BASE_DIR = "/home/lorenzo/Microbesomics"

ROW1 = [
    os.path.join(BASE_DIR, "results/figures/figure1a.png"),
    os.path.join(BASE_DIR, "results/figures/figure1b.png"),
    os.path.join(BASE_DIR, "results/figures/figure1c.png"),
    os.path.join(BASE_DIR, "results/figures/figure1d_mets_alluvial.png"),
]
ROW2 = [
    os.path.join(BASE_DIR, "results/figures/figure2a_pcoa_braycurtis.png"),
    os.path.join(BASE_DIR, "results/figures/figure2b_bc_distances.png"),
    os.path.join(BASE_DIR, "results/figures/figure2c_pcoa_unifrac.png"),
    os.path.join(BASE_DIR, "results/figures/figure2d_uuf_distances.png"),
]
ROW3 = [
    os.path.join(BASE_DIR, "results/figures/figure1i_composition.png"),
    os.path.join(BASE_DIR, "results/figures/figure1j_phylum_clr.png"),
]
LABELS_R1 = ["a", "b", "c", "d"]
LABELS_R2 = ["e", "f", "g", "h"]
LABELS_R3 = ["i", "j"]
OUT = os.path.join(BASE_DIR, "results/figures/figure1_combined.png")
DPI = 300


def read_rgba(path):
    img = mpimg.imread(path)
    if img.ndim == 2:
        img = np.stack([img]*3 + [np.ones_like(img)], axis=-1)
    elif img.shape[2] == 3:
        img = np.concatenate([img, np.ones((*img.shape[:2], 1), dtype=img.dtype)], axis=-1)
    return img.astype(np.float32)


def pad_height(img, target_h):
    h, w, c = img.shape
    if h >= target_h:
        return img
    top    = (target_h - h) // 2
    bottom = target_h - h - top
    white  = np.ones((1, w, c), dtype=img.dtype)
    return np.vstack([np.repeat(white, top, axis=0),
                      img,
                      np.repeat(white, bottom, axis=0)])


def pad_width(img, target_w):
    h, w, c = img.shape
    if w >= target_w:
        return img
    left  = (target_w - w) // 2
    right = target_w - w - left
    white = np.ones((h, 1, c), dtype=img.dtype)
    return np.hstack([np.repeat(white, left, axis=1),
                      img,
                      np.repeat(white, right, axis=1)])


def add_label(img, label):
    h, w = img.shape[:2]
    fig, ax = plt.subplots(figsize=(w/DPI, h/DPI), dpi=DPI)
    ax.imshow(img, aspect="auto", interpolation="lanczos")
    ax.text(0.012, 0.985, label, transform=ax.transAxes,
            fontsize=4, fontweight="bold", va="top", ha="left", color="#111")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.canvas.draw()
    buf = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8)
    buf = buf.reshape(fig.canvas.get_width_height()[::-1] + (4,))
    plt.close(fig)
    return buf[:, :, [1, 2, 3, 0]].astype(np.float32) / 255.0


# Pad to row height first, then add label — so the letter sits at the
# absolute top-left of the final padded image, aligned across all panels.
imgs_r1 = [read_rgba(p) for p in ROW1]
imgs_r2 = [read_rgba(p) for p in ROW2]
imgs_r3 = [read_rgba(p) for p in ROW3]

max_h_r1 = max(img.shape[0] for img in imgs_r1)
max_h_r2 = max(img.shape[0] for img in imgs_r2)
max_h_r3 = max(img.shape[0] for img in imgs_r3)
row1 = np.concatenate([add_label(pad_height(img, max_h_r1), lbl)
                        for img, lbl in zip(imgs_r1, LABELS_R1)], axis=1)
row2 = np.concatenate([add_label(pad_height(img, max_h_r2), lbl)
                        for img, lbl in zip(imgs_r2, LABELS_R2)], axis=1)
row3 = np.concatenate([add_label(pad_height(img, max_h_r3), lbl)
                        for img, lbl in zip(imgs_r3, LABELS_R3)], axis=1)

target_w = max(row1.shape[1], row2.shape[1], row3.shape[1])
row1 = pad_width(row1, target_w)
row2 = pad_width(row2, target_w)
row3 = pad_width(row3, target_w)

combined = np.vstack([row1, row2, row3])

fig, ax = plt.subplots(figsize=(combined.shape[1]/DPI, combined.shape[0]/DPI), dpi=DPI)
fig.patch.set_facecolor("white")
ax.imshow(combined, aspect="auto", interpolation="lanczos")
ax.axis("off")
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.savefig(OUT, dpi=DPI, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
