import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

# Load data
scatter_df = pd.read_excel("./5d.xlsx")
x_vals = scatter_df["FC Lymphocyte vs TLS"]
other_fc_cols = [
    "FC Fat tissue vs TLS", "FC in situ vs TLS", "FC Lactiferous duct vs TLS",
    "FC Necrosis vs TLS", "FC Tumor vs TLS", "FC Stroma vs TLS", "FC Vessels vs TLS"
]
y_vals = scatter_df[other_fc_cols].min(axis=1)

# Define TLS signature gene set with categories and colors
tls_signature_genes = {
    "CD79A": ("B cell", "royalblue"), "CD79B": ("B cell", "royalblue"),
    "TNFRSF13C": ("B cell", "royalblue"), "BLK": ("B cell", "royalblue"),
    "CD22": ("B cell", "royalblue"), "CD37": ("B cell", "royalblue"),
    "MS4A1": ("B cell", "royalblue"), "NIBAN3": ("B cell", "royalblue"),
    "CD19": ("B cell", "royalblue"), "IKZF3": ("B cell", "royalblue"),
    "LINC00926": ("B cell", "royalblue"),

    "FCRLA": ("Immunoglobulin", "orange"), "VPREB3": ("Immunoglobulin", "orange"),
    "FCMR": ("Immunoglobulin", "orange"), "AL928768.3": ("Immunoglobulin", "orange"),

    "CXCR5": ("Lymphoid nodule initation", "crimson"), "LTB": ("Lymphoid nodule initation", "crimson"),
    "CCL19": ("Lymphoid nodule initation", "crimson"),"POU2AF1": ("Lymphoid nodule initation", "crimson"),
    "CXCL13": ("Lymphoid nodule initation", "crimson"),

    "TCL1A": ("Cell survival", "lightgreen"),

    "SELL": ("Lymphoid nodule structure", "green"),

    "RASGRP2": ("T cell", "gold"), "TCF7": ("T cell", "gold"),
    "RIPOR2": ("T cell", "gold"),

    "RAC2": ("Immune", "mediumorchid"), "IL16": ("Immune", "mediumorchid"),
    "CCR7": ("Immune", "mediumorchid"), "CD52": ("Immune", "mediumorchid"),

    "ATP2A3": ("Contraction", "pink")
}

# Create plot
fig, ax = plt.subplots(figsize=(8, 8))  # Only keep this line

# Convert data to log2 space
x_log2 = np.log2(x_vals)
y_log2 = np.log2(y_vals)
ax.scatter(x_log2, y_log2, color='dodgerblue', alpha=0.1, s=10)

# Set axis limits and ticks
ax.set_xlim(-2, 5.3)
ax.set_ylim(-2, 3.3)
ax.set_xticks([-2, -1, 0, 1, 2, 3, 4, 5])
ax.set_yticks([-2, -1, 0, 1, 2, 3])
ax.set_xticklabels(["0","0.5", "1", "2", "4", "8", "16", "32"])
ax.set_yticklabels(["0", "0.5", "1", "2", "4", "8"])

# Reference fold-change lines
ax.axhline(0, color='gray', linewidth=0.8)
ax.axvline(0, color='gray', linewidth=0.8)

# Labels
ax.set_xlabel("Fold–change TLS vs. lymphocyte compartment")
ax.set_ylabel("Min fold–change TLS vs. another compartment")

# Ensure transform is initialized
fig.canvas.draw()
transform = ax.transAxes

# --- X axis (Top Arrows) ---
# Left: Lymphocytes
ax.annotate("", xy=(0, 1.03), xytext=(0.27, 1.03),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5))
ax.text(0.15, 1.05, "Lymphocytes", transform=transform, fontsize=12, ha='center')

# Right: TLS
ax.annotate("", xy=(1, 1.03), xytext=(0.28, 1.03),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5))
ax.text(0.65, 1.05, "TLS", transform=transform, fontsize=12, ha='center')

# --- Y axis (Right Arrows) ---
# Bottom: Other non-lymphocytes
ax.annotate("", xy=(1.03, 0), xytext=(1.03, 0.37),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5))
ax.text(1.05, 0.18, "Other non-lymphocytes", transform=transform, fontsize=12, va='center', rotation=270)

# Top: TLS
ax.annotate("", xy=(1.03, 1), xytext=(1.03, 0.38),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5))
ax.text(1.05, 0.65, "TLS", transform=transform, fontsize=12, va='center', rotation=270)

# Plot TLS signature genes
for gene, (category, color) in tls_signature_genes.items():
    match = scatter_df[scatter_df["Gene"] == gene]
    if not match.empty:
        x = np.log2(match["FC Lymphocyte vs TLS"].values[0])
        y = np.log2(match[other_fc_cols].min(axis=1).values[0])
        ax.plot(x, y, marker='o', markersize=8, markeredgecolor=color, markerfacecolor='none')
        ax.annotate(gene, (x, y), textcoords="offset points", xytext=(5, 5), ha='left',
                    fontsize=8, color=color, arrowprops=dict(arrowstyle='-', color=color, lw=0.8))

# Add legend
legend_elements = []
seen = set()
for _, (category, color) in tls_signature_genes.items():
    if category not in seen:
        legend_elements.append(Patch(facecolor=color, edgecolor='black', label=category))
        seen.add(category)

ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(0.6, 0.05), fontsize=10)
plt.tight_layout()
plt.show()