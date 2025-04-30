import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.colors import to_rgb
from matplotlib.colors import LinearSegmentedColormap

# Load your real Excel data
file_path = "./4e.xlsx"
df_full = pd.read_excel(file_path, index_col=0)

# Extract meta and data
meta = df_full.iloc[:, :2]  # 'Stroma subtype' and 'PB subtype'
data = df_full.iloc[:, 2:]  # Actual expression data

# Your desired sample IDs
sample_ids_m = [7, 8, 21, 22, 24, 25, 36, 68, 79, 81, 96]
sample_ids_msl = [2, 12, 16, 19, 23, 32, 34, 47, 48, 49, 53, 66, 71, 80, 85, 86]
sample_ids = sample_ids_m + sample_ids_msl

# Your desired feature names
signatures = ['IL-iCAF', 'IFNγ-iCAF', 'Inflammatory response', 'IL6 JAK STAT3 signaling', 'Immune2',
              'TLS Cabrita', 'TLS Meylan', 'TLS Lundeberg', 'TLS ST', 'Apoptosis', 'VCpred TN']
cell_types = ['B-cells', 'Memory B-cells', 'pro B-cells', 'CD4+ naive T-cells', 'aDC', 'pDC']
single_genes = ['VISTA', 'CCL5', 'CD134', 'CD27', 'CD40L', 'CSF1R', 'FOXP3', 'PD1', 'TIGIT']
feature_names = signatures + cell_types + single_genes

# Transpose so samples become columns
data = data.T
data = data.loc[feature_names, sample_ids]
df_scaled = (data - data.mean(axis=1).values[:, None]) / data.std(axis=1).values[:, None]

# Layout
fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(12, 14, wspace=0.05, hspace=0.05)

# Top Labels
ax_label = fig.add_subplot(gs[0, 1:-2])
ax_label.axis('off')
ax_label.text(5.5/26, 0.8, 'M-M$_{\\mathrm{STROMA}}$', ha='center', va='bottom', fontsize=14, color='black', transform=ax_label.transAxes)
ax_label.text((11+7.5)/26, 0.8, 'M-MSL$_{\\mathrm{STROMA}}$', ha='center', va='bottom', fontsize=14, color='#C34E79', transform=ax_label.transAxes)

# Heatmap
colors = ['blue', 'gray', 'yellow']
cmap = LinearSegmentedColormap.from_list('custom_blue_gray_yellow', colors, N=256)
ax_heatmap = fig.add_subplot(gs[1:10, 1:-2])
sns.heatmap(
    df_scaled,
    cmap = cmap,
    center=0,
    xticklabels=True,
    yticklabels=True,
    cbar=False,
    ax=ax_heatmap
)

# White Lines
signature_end = len(signatures)
celltype_end = signature_end + len(cell_types)
midpoint = len(sample_ids_m)

ax_heatmap.hlines(signature_end, *ax_heatmap.get_xlim(), colors='white', linewidth=2)
ax_heatmap.hlines(celltype_end, *ax_heatmap.get_xlim(), colors='white', linewidth=2)
ax_heatmap.vlines(midpoint, *ax_heatmap.get_ylim(), colors='white', linewidth=2)

# X ticks
ax_heatmap.set_xticks(np.arange(len(sample_ids)) + 0.5)
ax_heatmap.set_xticklabels(sample_ids, rotation=0, fontsize=6)
ax_heatmap.tick_params(axis='x', length=0)

# Y ticks
ax_heatmap.set_yticks(np.arange(len(feature_names)) + 0.5)
ax_heatmap.set_yticklabels(feature_names, rotation=0, fontsize=10)
ax_heatmap.tick_params(axis='y', length=0)

# Right group colorbar
ax_colors = fig.add_subplot(gs[1:10, -2])
group_colors_text = ['purple'] * len(signatures) + ['gray'] * len(cell_types) + ['green'] * len(single_genes)
group_colors_rgb = np.array([to_rgb(color) for color in group_colors_text]).reshape(len(group_colors_text), 1, 3)
ax_colors.imshow(group_colors_rgb, aspect='auto')
ax_colors.axis('off')

# Right group text
ax_grouptext = fig.add_subplot(gs[1:10, -1])
ax_grouptext.axis('off')
ax_grouptext.text(0, (len(signatures)/2)/len(feature_names), 'Signatures', rotation=270, va='center', ha='center', fontsize=12, color='black')
ax_grouptext.text(0, (len(signatures)+len(cell_types)/2)/len(feature_names), 'Cell types', rotation=270, va='center', ha='center', fontsize=12, color='black')
ax_grouptext.text(0, (len(signatures)+len(cell_types)+len(single_genes)/2)/len(feature_names), 'Single genes', rotation=270, va='center', ha='center', fontsize=12, color='black')
ax_grouptext.invert_yaxis()

# Bottom Colorbar
ax_cb = fig.add_subplot(gs[11, 3:-3])
norm = plt.Normalize(vmin=df_scaled.values.min(), vmax=df_scaled.values.max())
sm = plt.cm.ScalarMappable(cmap= cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=ax_cb, orientation='horizontal')
cb.set_ticks([df_scaled.values.min(), df_scaled.values.max()])
cb.set_ticklabels(['Low', 'High'])
cb.ax.tick_params(labelsize=10)

plt.show()