from source import sampling as spl
import sys
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 500)


model = spl.load_pgm()

print(model.edges())
nms = spl.get_var_nms()


sample_size = 10000
df_sim = spl.random_sample(model, nms, sample_size)
#print(df_sim.iloc[:10])
df_sim.to_csv('~/Documents/workspace/causal-eval-rct/neuropathic-pain-10000-2.csv', index=False)
"""
# print(df_sim.shape)

df_sim_num = spl.nam2num(df_sim)
df_sim_num.iloc[:10]

df_sim_mcar = spl.add_missing_data(df_sim, mode='mcar', seed=10, mcar_p=0.1)
#print(df_sim_mcar.iloc[:10])
print(df_sim_mcar.shape)
#for col in df_sim.columns:
#    print(col)

df_sim_sb = spl.add_selection_bias(df_sim, prob=0.1, seed=10)
df_sim_sb.iloc[:10]
print(df_sim_sb.shape)

df_sim_uc = spl.add_confounder(df_sim)
df_sim_uc.iloc[:10]
print(df_sim_uc.shape)

sel_var = [35, 36, 37, 38, 73, 74, 75, 76]
sel_list = np.array(df_sim.iloc[:, sel_var])#.sum(axis=1) > 6)  # True will be selected
print(sel_list[:10])

sel_list = np.array(df_sim.iloc[:, sel_var].sum(axis=1) > 4)  # True will be selected
print(sel_list[:10])

nrow = 100
prob = 0.9
del_list = np.random.choice([0, 1], size=nrow, p=[prob, 1 - prob]) == 1  # 0 missing false,that will not be selected
print(del_list)
"""