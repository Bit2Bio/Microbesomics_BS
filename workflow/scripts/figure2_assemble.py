"""
Assembles the four beta diversity panels into results/figures/figure2_combined.png.
Panels are centre-padded vertically to a common height before concatenation.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

BASE_DIR = "/home/lorenzo/Microbesomics"

PANELS = [
    os.path.join(BASE_DIR, "results/figures/figure2a_pcoa_braycurtis.png"),
    os.path.join(BASE_DIR, "results/figures/figure2b_bc_distances.png"),
    os.path.join(BASE_DIR, "results/figures/figure2c_pcoa_unifrac.png"),
    os.path.join(BASE_DIR, "results/figures/figure2d_uuf_distances.png"),
]
LABELS = ["a", "b", "c", "d"]
OUT    = os.path.join(BASE_DIR, "results/figures/figure2_combined.png")
DPI    = 300


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


imgs  = [read_rgba(p) for p in PANELS]
imgs  = [add_label(img, lbl) for img, lbl in zip(imgs, LABELS)]
max_h = max(img.shape[0] for img in imgs)
row   = np.concatenate([pad_height(img, max_h) for img in imgs], axis=1)

fig, ax = plt.subplots(figsize=(row.shape[1]/DPI, row.shape[0]/DPI), dpi=DPI)
fig.patch.set_facecolor("white")
ax.imshow(row, aspect="auto", interpolation="lanczos")
ax.axis("off")
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.savefig(OUT, dpi=DPI, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
