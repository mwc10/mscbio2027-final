import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# Load your data
file_path = "./4f.xlsx"
df = pd.read_excel(file_path)
df = df.dropna(subset=['DRFS.time', 'DRFS.event'])

# Split groups
group1 = df[df['Class'] == 'M-Mstroma']
group2 = df[df['Class'] == 'M-MSLstroma']

# Plot setup
fig, ax = plt.subplots(figsize=(6, 5))
kmf1 = KaplanMeierFitter()
kmf2 = KaplanMeierFitter()

kmf1.fit(group1['DRFS.time'], group1['DRFS.event'], label='M-M$_{\\mathrm{STROMA}}$')
kmf2.fit(group2['DRFS.time'], group2['DRFS.event'], label='M-MSL$_{\\mathrm{STROMA}}$')

kmf1.plot_survival_function(ax=ax, ci_show=False, show_censors=True,
                            censor_styles={'marker': '|', 'ms': 8}, color='black', linewidth=1.5)
kmf2.plot_survival_function(ax=ax, ci_show=False, show_censors=True,
                            censor_styles={'marker': '|', 'ms': 8}, color='#C34E79', linewidth=1.5)

# Axes formatting
ax.set_xlabel('Time (yr)', fontsize=12)
ax.set_ylabel('DRFS', fontsize=12)
ax.set_xlim(0, 10.2)
ax.set_ylim(0, 1.02)
ax.set_title('M subtype in global PB', fontsize=14, weight='bold')
ax.legend().remove()

# Permutation log-rank test
T = pd.concat([group1['DRFS.time'], group2['DRFS.time']])
E = pd.concat([group1['DRFS.event'], group2['DRFS.event']])
labels = np.array([1]*len(group1) + [0]*len(group2))

# Observed stat
obs = logrank_test(T[labels == 1], T[labels == 0], E[labels == 1], E[labels == 0])
obs_stat = obs.test_statistic

# Permutation
n_perms = 1000
perm_stats = []
for _ in range(n_perms):
    perm = np.random.permutation(labels)
    stat = logrank_test(T[perm == 1], T[perm == 0], E[perm == 1], E[perm == 0]).test_statistic
    perm_stats.append(stat)

p_perm = np.mean(np.array(perm_stats) >= obs_stat)

# Add permutation p-value to plot
ax.text(7, 0.2, f'$p$ = {p_perm:.3f}', fontsize=12, style='italic')

# Number at risk table
times = np.arange(0, 11, 2)
at_risk_1 = [np.sum(group1['DRFS.time'] >= t) for t in times]
at_risk_2 = [np.sum(group2['DRFS.time'] >= t) for t in times]

# Below plot
for i, (n1, n2) in enumerate(zip(at_risk_1, at_risk_2)):
    xpos = times[i]
    # numbers stay centered under time points
    ax.text(xpos, -0.25, str(n1), ha='center', va='center', fontsize=10, color='black')
    ax.text(xpos, -0.35, str(n2), ha='center', va='center', fontsize=10, color='black')

# Then separately plot the group labels at the far left
ax.text(-1.2, -0.25, 'M-M$_{\\mathrm{STROMA}}$', ha='right', va='center', fontsize=11, color='black', weight='bold')
ax.text(-1.2, -0.35, 'M-MSL$_{\\mathrm{STROMA}}$', ha='right', va='center', fontsize=11, color='#C34E79', weight='bold')

# Adjust layout to make room
plt.subplots_adjust(left=0.3, bottom=0.3)
plt.show()