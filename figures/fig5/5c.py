import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load and clean data
df = pd.read_excel("../selected_GO.xlsx")
df["GO"] = df["GO"].str.replace("GO_", "", regex=False).str.replace("_", " ").str.lower().str.title()

# Calculate 1/FDR for bar length (so bars grow with significance)
df["Bar Length"] = 1 / df["FDR TLS > Lymphocyte"]

# Sort by FDR
df_sorted = df.sort_values("FDR TLS > Lymphocyte", ascending=True)

# Plot
plt.figure(figsize=(12, 10))
sns.barplot(data=df_sorted, y="GO", x="Bar Length", color="skyblue")

# Set x-axis to log scale and inverted direction
plt.xscale("log")
plt.xlim(1, 1e15) 
plt.xlabel("FDR")
plt.ylabel("")
plt.xticks([1, 1e5, 1e10, 1e15], ["1", "10⁻⁵", "10⁻¹⁰", "10⁻¹⁵"])  # Custom tick labels
plt.tight_layout()
plt.show()
