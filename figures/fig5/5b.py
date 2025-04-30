import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# Load the data
df = pd.read_excel("./5b.xlsx")

# Define ordered cell types and label colors
ordered_cell_types = [
    "B-cells", "Class-switched memory B-cells", "Memory B-cells", "naive B-cells", "pro B-cells",
    "CD4+ T-cells", "CD4+ memory T-cells", "CD4+ naive T-cells", "CD4+ Tcm", "CD4+ Tem", "Th1 cells",
    "DC", "aDC", "Macrophages", "Mast cells", "Monocytes", "Endothelial cells", "Myocytes", "MSC", "Neurons"
]

label_colors = [
    "goldenrod", "goldenrod", "goldenrod", "goldenrod", "goldenrod",
    "mediumseagreen", "mediumseagreen", "mediumseagreen", "mediumseagreen", "mediumseagreen", "mediumseagreen",
    "red", "red", "crimson", "crimson", "crimson",
    "orange", "orange", "purple", "skyblue"
]

# Statistical testing and reshape
plot_data = []
p_values = []

for cell_type in ordered_cell_types:
    tls_vals = df[df['Compartment'] == 'TLS'][cell_type].values
    lymph_vals = df[df['Compartment'] == 'Lymphocytes'][cell_type].values
    _, p = ranksums(tls_vals, lymph_vals)
    p_values.append(p)

    for val in tls_vals:
        plot_data.append({'Cell Type': cell_type, 'Compartment': 'TLS', 'Score': val})
    for val in lymph_vals:
        plot_data.append({'Cell Type': cell_type, 'Compartment': 'Lymphocyte compartment', 'Score': val})

# FDR correction
_, fdr_corrected_p, _, _ = multipletests(p_values, method='fdr_bh')

# Create DataFrame
plot_df = pd.DataFrame(plot_data)
plot_df["Cell Type"] = pd.Categorical(plot_df["Cell Type"], categories=ordered_cell_types, ordered=True)

# Plot
plt.figure(figsize=(10, 6))
ax = sns.violinplot(
    data=plot_df, x="Score", y="Cell Type", hue="Compartment",
    palette={"TLS": "cornflowerblue", "Lymphocyte compartment": "lightsalmon"},
    split=True, inner=None, linewidth=1.2, scale="width", cut=0)

# Add significance asterisks
for i, p in enumerate(fdr_corrected_p):
    if p < 0.05:
        if p < 0.0001:
            stars = "****"
        elif p < 0.001:
            stars = "***"
        elif p < 0.01:
            stars = "**"
        else:
            stars = "*"
        plt.text(-0.01, i, stars, va='center', ha='right', fontsize=8)

# Annotate median lines
for i, cell_type in enumerate(ordered_cell_types):
    for j, compartment in enumerate(['TLS', 'Lymphocyte compartment']):
        median_val = plot_df[(plot_df['Cell Type'] == cell_type) & (plot_df['Compartment'] == compartment)]['Score'].median()
        ax.plot([median_val, median_val], [i - 0.3 + j*0.6, i - 0.1 + j*0.6], color='black', linewidth=1)

# Color the y-axis tick labels
y_labels = ax.get_yticklabels()
for tick_label, color in zip(y_labels, label_colors):
    tick_label.set_color(color)

# Set x-axis range and ticks to match example
ax.set_xlim(-0.08, 1.75)
ax.set_xticks([0.0, 0.5, 1.0, 1.5])
plt.xlabel("xCell enrichment score", fontsize=12)
plt.ylabel("")
plt.title("")
plt.legend(
    title=None,
    loc='lower left',
    bbox_to_anchor=(0.55, 0.65),  # adjust to fine-tune position
    frameon=False,
    fontsize=12
)
plt.tight_layout()
plt.show()