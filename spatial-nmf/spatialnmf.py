#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 17:07:54 2025

"""

import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from PIL import Image


expr_path   = "./TNBC1_tumor_expr.csv"    # spots × genes
coords_path = "./TNBC1_tumor_coords.csv"  # pixel_x, pixel_y
he_path     = "./TNBC1_CN1_C1.jpg"  # H&E image


expr   = pd.read_csv(expr_path, index_col=0)
coords = pd.read_csv(coords_path, index_col=0)

# Drop mitochondrial genes (MT-)
expr = expr.loc[:, ~expr.columns.str.upper().str.startswith("MT-")]

#Drop ribosomal genes (RPS*, RPL*)
expr = expr.loc[:, ~expr.columns.str.upper()
                      .str.startswith(("RPS","RPL"))]

# Keep top 2,000 highly variable genes
hvg = expr.var(axis=0).nlargest(2000).index
expr_hvg = expr[hvg]

# Drop known housekeeping/lncRNA to avoid technical factors
expr_hvg = expr_hvg.drop(columns=["MALAT1", "ACTB", "GAPDH"], errors="ignore")
expr_log = np.log1p(expr_hvg)

# NMF at k=4
k = 4
model = NMF(n_components=k, init="nndsvda", random_state=42, max_iter=10000)
W     = model.fit_transform(expr_log)  # spots × factors
H     = model.components_             # factors × genes
err   = model.reconstruction_err_

resid_sq = err**2
total_sq = np.linalg.norm(expr_log.values, ord="fro")**2
var_exp  = 100 * (1 - resid_sq / total_sq)

print(f"Reconstruction error (k={k}): {err:.4f}")
print(f"Variance explained  (k={k}): {var_exp:.1f}%\n")


he_img = Image.open(he_path)

selected_idxs = [0, 1, 3]
labels = ["Factor 1", "Factor 2", "Factor 3"]

# fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
# for ax, idx, lbl in zip(axes, selected_idxs, labels):
#     rx, _ = spearmanr(W[:, idx], coords["pixel_x"])
#     ry, _ = spearmanr(W[:, idx], coords["pixel_y"])
#     axis = "pixel_y" if abs(ry) >= abs(rx) else "pixel_x"
#     order = np.argsort(coords[axis])
#     ax.plot(W[order, idx], lw=2)
#     ax.set_title(f"{lbl}\n(sorted by {axis})")
#     ax.set_xlabel(axis)
# axes[0].set_ylabel("Activity")
# plt.tight_layout()
# plt.show()


# 2D Scatter Plots
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True, sharey=True)
for ax, idx, lbl in zip(axes, selected_idxs, labels):
    sc = ax.scatter(coords["pixel_x"], coords["pixel_y"], c=W[:, idx], s=12)
    ax.set_title(lbl)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.axis("off")
    fig.colorbar(sc, ax=ax, shrink=0.7)
plt.tight_layout()
plt.show()

# Quiver Vector Fields
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
x = coords["pixel_x"].values
y = coords["pixel_y"].values
for ax, idx, lbl in zip(axes, selected_idxs, labels):
    z = W[:, idx]
    gx, gy = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
    gz = griddata((x, y), z, (gx, gy), method="cubic")
    gz = gaussian_filter(gz, sigma=2)
    dz_dx, dz_dy = np.gradient(gz)
    pcm = ax.imshow(gz.T, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower', alpha=0.8)
    ax.quiver(gx, gy, dz_dx, dz_dy, scale=30, width=0.002)
    ax.set_title(lbl)
    ax.invert_yaxis(); ax.axis("off")
    fig.colorbar(pcm, ax=ax, shrink=0.7)
plt.tight_layout()
plt.show()

#H&E Overlays
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, idx, lbl in zip(axes, selected_idxs, labels):
    z = W[:, idx]
    gx, gy = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
    gz = griddata((x, y), z, (gx, gy), method="cubic")
    gz = gaussian_filter(gz, sigma=2)
    dz_dx, dz_dy = np.gradient(gz)
    ax.imshow(he_img, extent=(x.min(), x.max(), y.max(), y.min()))
    ax.imshow(gz.T, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower',
              cmap="inferno", alpha=0.4)
    ax.quiver(gx, gy, dz_dx, -dz_dy, scale=30, color='cyan', width=0.003, headwidth=3, alpha=0.7)
    ax.set_title(f"{lbl} overlay")
    ax.axis("off")
plt.tight_layout()
plt.show()