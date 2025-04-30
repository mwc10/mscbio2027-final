import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import kruskal, ranksums
from statsmodels.stats.multitest import multipletests

# Load the data
df = pd.read_excel("./5g.xlsx")
plot_df = df[['TLS ST signature', 'TIME']].dropna()

# Define colors for TIME classifications
time_colors = {
    'FI': '#D7301F',  # strong red
    'SR': '#EF6548',  # orange-red
    'MR': '#FDBB84',  # light brown
    'ID': '#BDBDBD'   # grey
}

# Run Kruskal–Wallis test
groups = [plot_df[plot_df['TIME'] == time]['TLS ST signature'] for time in time_colors]
kw_stat, kw_p = kruskal(*groups)

# Wilcoxon rank-sum tests: each subtype vs. all others
p_vals = []
for time in time_colors:
    this_group = plot_df[plot_df['TIME'] == time]['TLS ST signature']
    other_group = plot_df[plot_df['TIME'] != time]['TLS ST signature']
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

# Plot KDE and vertical means
for i, (time_class, color) in enumerate(time_colors.items()):
    subset = plot_df[plot_df['TIME'] == time_class]
    sns.kdeplot(data=subset, x='TLS ST signature', color=color, linewidth=1.8, clip=(-0.8, 0.95))
    mean_val = subset['TLS ST signature'].mean()
    plt.axvline(mean_val, color=color, linestyle='--', linewidth=1)

    star = get_star(fdr_vals[i])
    if star:
        plt.text(mean_val, 2, star, fontsize=12, ha='center')

# Custom legend using square patches
legend_elements = [
    Patch(facecolor=color, edgecolor='black', label=label)
    for label, color in time_colors.items()
]
plt.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.15),
    ncol=4, frameon=False, fontsize=10, handleheight=1.2, handlelength=1.2)

# Axes and ticks
plt.xlim(-0.85, 1)
plt.ylim(-0.1, 2.2)
plt.xticks([-0.5, 0.0, 0.5, 1.0], fontsize=10)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0], fontsize=10)
plt.xlabel("TLS signature", fontsize=11)
plt.ylabel("Density", fontsize=11)
plt.text(-0.9, 2.6, "TIME classification", fontsize=13, weight='bold')
plt.text(-0.9, -0.7, f"Kruskall–Wallis, p = {format_pval(kw_p)}", fontsize=10)
plt.tight_layout()
plt.show()