import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import kruskal, ranksums
from statsmodels.stats.multitest import multipletests

# Load the data
df = pd.read_excel("./5f.xlsx")
plot_df = df[['TLS ST signature', 'Subtype']].dropna()

# Define subtype colors
subtype_colors = {
    'LAR': '#E74C3C',       # red
    'MSL': '#F39C12',       # orange
    'M': '#F1C40F',         # yellow
    'BL': '#3498DB',        # blue
    'IM': '#27AE60'         # green
}

# Kruskal–Wallis test across all subtypes
groups = [plot_df[plot_df['Subtype'] == subtype]['TLS ST signature'] for subtype in subtype_colors]
kw_stat, kw_p = kruskal(*groups)

# Wilcoxon rank-sum tests: each subtype vs. all others
p_vals = []
for subtype in subtype_colors:
    this_group = plot_df[plot_df['Subtype'] == subtype]['TLS ST signature']
    other_group = plot_df[plot_df['Subtype'] != subtype]['TLS ST signature']
    stat, p = ranksums(this_group, other_group)
    p_vals.append(p)

# FDR correction
_, fdr_vals, _, _ = multipletests(p_vals, method='fdr_bh')

# Significance levels
def get_star(fdr):
    if fdr < 0.0001:
        return '****'
    elif 0.0001 <= fdr and fdr < 0.001:
        return '***'
    elif 0.001 <= fdr and fdr < 0.01:
        return '**'
    elif 0.01 <= fdr and fdr < 0.05:
        return '*'
    else:
        return ''

# Format p-value as a × 10⁻⁸
def format_pval(p):
    exponent = int(f"{p:.1e}".split('e')[1])
    base = float(f"{p:.1e}".split('e')[0])
    return f"{base:.1f} × 10$^{{{exponent}}}$"

# Start plotting
plt.figure(figsize=(6, 4))
sns.set(style="white", context="paper")

# Plot KDE and means
for i, (subtype, color) in enumerate(subtype_colors.items()):
    subset = plot_df[plot_df['Subtype'] == subtype]
    sns.kdeplot(data=subset, x='TLS ST signature', color=color, linewidth=1.8, clip=(-0.47, 0.95))
    mean_val = subset['TLS ST signature'].mean()
    plt.axvline(mean_val, color=color, linestyle='--', linewidth=1)
    
    star = get_star(fdr_vals[i])
    if star:
        plt.text(mean_val, 1.48, star, fontsize=12, ha='center')

# Custom legend using patches
legend_elements = [
    Patch(facecolor=color, edgecolor='black', label=label)
    for label, color in subtype_colors.items()
]
plt.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.15),
           ncol=5, frameon=False, fontsize=10, handleheight=1.2, handlelength=1.2)

# Labels and formatting
plt.xlim(-0.5, 1.0)
plt.ylim(-0.1, 1.6)
plt.xticks([-0.5, 0.0, 0.5, 1.0], fontsize=10)
plt.yticks([0.0, 0.5, 1.0, 1.5], fontsize=10)
plt.xlabel("TLS signature", fontsize=11)
plt.ylabel("Density", fontsize=11)
plt.text(-0.55, 1.9, "TNBC molecular subtypes", fontsize=13, weight='bold')
plt.text(-0.55, -0.5, f"Kruskall–Wallis, p = {format_pval(kw_p)}", fontsize=10)
plt.tight_layout()
plt.show()