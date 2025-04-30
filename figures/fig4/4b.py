import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# Load your mapping data
file_path = "./4b.xlsx"
df = pd.read_excel(file_path)
df = df.dropna(subset=['TIME'])  # Drop missing TIME samples

# Column order
columns = ['PB', 'Tumor', 'Stroma', 'TIME']

# Color map
color_map = {
    'LAR': 'indianred',
    'MSL': 'orange',
    'M': 'yellow',
    'BL': 'royalblue',
    'IM': 'mediumseagreen',
    'ID': 'gray',
    'MR': 'lightcoral',
    'SR': 'salmon',
    'FI': 'red'
}

# Define the desired order of labels for each column
category_order_map = {
    'PB': ['LAR', 'MSL', 'M', 'BL', 'IM'],
    'Tumor': ['LAR', 'MSL', 'M', 'BL', 'IM'],
    'Stroma': ['LAR', 'MSL', 'M', 'BL', 'IM'],
    'TIME': ['ID', 'MR', 'SR', 'FI']
}

# Layout constants
x_gap = 2.0
column_x = {col: i * x_gap for i, col in enumerate(columns)}
box_width = 0.2
fig, ax = plt.subplots(figsize=(10, 8))
ax.axis('off')

# Scale factor
total_count = len(df)
scale = 1.0 / total_count

# Collect labels and assign block heights
label_positions = {}
block_heights = {}

for col in columns:
    label_positions[col] = {}
    block_heights[col] = {}
    y = 1.0
    ordered = [lbl for lbl in category_order_map[col] if lbl in df[col].values]
    for lbl in ordered:
        height = df[df[col] == lbl].shape[0] * scale
        label_positions[col][lbl] = (y - height, y)
        block_heights[col][lbl] = height
        rect = patches.Rectangle((column_x[col] - box_width/2, y - height), box_width, height,
                                 facecolor=color_map.get(lbl, 'lightgray'), edgecolor='black')
        ax.add_patch(rect)
        ax.text(column_x[col], y - height / 2, lbl, ha='center', va='center', fontsize=10)
        y -= height

# Draw flows
for i in range(len(columns) - 1):
    src_col = columns[i]
    tgt_col = columns[i + 1]

    # Calculate flow counts
    flow_data = df.groupby([src_col, tgt_col]).size().reset_index(name='count')
    flow_data = flow_data.sort_values(by=[src_col], key=lambda x: x.map({k: i for i, k in enumerate(category_order_map[src_col])}))

    # Track offsets within each label block
    src_offsets = {lbl: label_positions[src_col][lbl][0] for lbl in label_positions[src_col]}
    tgt_offsets = {lbl: label_positions[tgt_col][lbl][0] for lbl in label_positions[tgt_col]}

    for _, row in flow_data.iterrows():
        src_lbl, tgt_lbl = row[src_col], row[tgt_col]
        count = row['count']
        thickness = count * scale

        src_x = column_x[src_col] + box_width / 2
        tgt_x = column_x[tgt_col] - box_width / 2

        src_y0 = src_offsets[src_lbl]
        src_y1 = src_y0 + thickness
        src_offsets[src_lbl] = src_y1

        tgt_y0 = tgt_offsets[tgt_lbl]
        tgt_y1 = tgt_y0 + thickness
        tgt_offsets[tgt_lbl] = tgt_y1

        # Build curved path
        path_data = [
            (Path.MOVETO, (src_x, src_y0)),
            (Path.CURVE4, (src_x + 0.5, src_y0)),
            (Path.CURVE4, (tgt_x - 0.5, tgt_y0)),
            (Path.CURVE4, (tgt_x, tgt_y0)),

            (Path.LINETO, (tgt_x, tgt_y1)),
            (Path.CURVE4, (tgt_x - 0.5, tgt_y1)),
            (Path.CURVE4, (src_x + 0.5, src_y1)),
            (Path.CURVE4, (src_x, src_y1)),

            (Path.CLOSEPOLY, (src_x, src_y0))
        ]

        codes, verts = zip(*path_data)
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor=color_map.get(src_lbl, 'lightgray'), edgecolor='none', alpha=0.4)
        ax.add_patch(patch)

# Add column titles
titles = {
    'PB': 'Global\npseudobulk',
    'Tumor': 'Tumor\npseudobulk',
    'Stroma': 'Stroma\npseudobulk',
    'TIME': 'TIME\nclassification'
}
for col in columns:
    ax.text(column_x[col], 1.08, titles[col], ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_xlim(-0.5, x_gap * (len(columns) - 1) + 0.5)
ax.set_ylim(0, 1.1)
plt.tight_layout()
plt.show()
