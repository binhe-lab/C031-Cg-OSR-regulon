import pandas as pd

df = pd.read_csv(
    '../data/normalized_data_oxidative_stress_Candida_glabrata.txt.csv')


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

h2o2_targets_df = df.loc[h2o2_targets, :]
h2o2_non_yap1skn7_dependent_targets = h2o2_targets_df[abs(h2o2_targets_df.iloc[:, 1]) < 0.1]
h2o2_non_yap1skn7_dependent_targets = list(h2o2_non_yap1skn7_dependent_targets.index)

mena_targets_df = df.loc[mena_targets, :]
mena_non_yap1skn7_dependent_targets = mena_targets_df[abs(mena_targets_df.iloc[:, 1]) < 0.1]
mena_non_yap1skn7_dependent_targets = list(mena_non_yap1skn7_dependent_targets.index)

h2o2_mena_non_yap1skn7_dependent_targets_overlap = list(
    set(h2o2_non_yap1skn7_dependent_targets)
    & set(mena_non_yap1skn7_dependent_targets))

print(len(mena_non_yap1skn7_dependent_targets))
