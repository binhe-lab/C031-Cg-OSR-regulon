# %%

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2_unweighted

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
axes[0].set_title('H2O2 and Menadione stress targets', fontsize=(20))
axes[1].set_title('Non-Yap1Skn7 dependent targets', fontsize=(20))

df = pd.read_csv(
    '01-gene-expression/data/normalized_data_oxidative_stress_Candida_glabrata.txt.csv')

cols_wanted = []
h2o2_stress_unstress = 'E01-07-WT.gpr;E01-09-WT.gpr'
cols_wanted.append(h2o2_stress_unstress)

h2o2_yap1skn7_stress = 'E01-08-dd.gpr;E01-10-dd.gpr'
cols_wanted.append(h2o2_yap1skn7_stress)

mena_stress_unstress = 'mena_wt.gpr;mena_wt22.gpr'
cols_wanted.append(mena_stress_unstress)

mena_yap1skn7_stress = 'mena_sy.gpr;mena_sy22.gpr'
cols_wanted.append(mena_yap1skn7_stress)

df = df.drop(0)
df.index = df['Scan REF']

df = df.loc[:, cols_wanted]

df.columns = ['h2o2_stress_unstress', 'h2o2_yap1skn7_stress',
              'mena_stress_unstress', 'mena_yap1skn7_stress']

df = df.astype(float)

h2o2_targets = df[df.iloc[:, 0] > 1]
h2o2_targets = list(h2o2_targets.index)

mena_targets = df[df.iloc[:, 2] > 1]
mena_targets = list(mena_targets.index)

h2o2_mena_overlap_targets = list(set(h2o2_targets) & set(mena_targets))

venn1 = venn2_unweighted(subsets=(len(h2o2_targets), len(mena_targets),
                         len(h2o2_mena_overlap_targets)), ax=axes[0],
                         set_labels=('H2O2 stress targets', 'Menadione stress targets'))

for text in venn1.set_labels:
    text.set_fontsize(15)

for text in venn1.subset_labels:
    text.set_fontsize(20)

h2o2_targets_df = df.loc[h2o2_targets, :]
h2o2_non_yap1skn7_dependent_targets = h2o2_targets_df[abs(h2o2_targets_df.iloc[:, 1]) < 0.3]
h2o2_non_yap1skn7_dependent_targets = list(h2o2_non_yap1skn7_dependent_targets.index)

mena_targets_df = df.loc[mena_targets, :]
mena_non_yap1skn7_dependent_targets = mena_targets_df[abs(mena_targets_df.iloc[:, 1]) < 0.3]
mena_non_yap1skn7_dependent_targets = list(mena_non_yap1skn7_dependent_targets.index)

h2o2_mena_non_yap1skn7_dependent_targets_overlap = list(
    set(h2o2_non_yap1skn7_dependent_targets)
    & set(mena_non_yap1skn7_dependent_targets))

venn2 = venn2_unweighted(subsets=(len(h2o2_non_yap1skn7_dependent_targets),
                         len(mena_non_yap1skn7_dependent_targets),
                         len(h2o2_mena_non_yap1skn7_dependent_targets_overlap)),
                         ax=axes[1],
                         set_labels=('H2O2 targets non-dependent', 'Menadione targets non-dependent'))

for text in venn2.set_labels:
    text.set_fontsize(15)

for text in venn2.subset_labels:
    text.set_fontsize(20)

# Uncomment below to save
# fig.savefig('01-gene-expression/output/fig2c.png')

# %%

table = pd.read_csv('01-gene-expression/analysis/meta-table-sample.txt')

chip_cols = table[table.loc[:, 'Lelandais log2 Chip-Seq Rep1']]
print(chip_cols)
