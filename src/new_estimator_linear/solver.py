import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

df = pd.read_csv('test_observations2', sep=' ')
i_list = df['E_I'].tolist()
n_list = df['N_shared'].tolist()
d_list = df['E_D'].tolist()
s_list = df['E_S'].tolist()
f_A_mut_list = df['f_A_mut'].tolist()
k_list = df['k'].tolist()
L_list = df['L'].tolist()
p_s_list = df['p_s'].tolist()
L2_list = df['L2'].tolist()

"""
solutions = []
d_value_to_solutions={}
for i,h,d in list(zip(i_list, n_list, d_list)):
    if d not in d_value_to_solutions.keys():
        d_value_to_solutions[d] = []
    try:
        d_value_to_solutions[d].append( 1.0*i/( h*(k-1) - i ) )
    except:
        break
for d in d_value_to_solutions.keys():
    print(d, np.mean(d_value_to_solutions[d]))
"""

mut_rate_to_fraction_1_values = {}
for s,d,f_A_mut,L,mut_rate in list( zip(s_list, d_list, f_A_mut_list, L_list, p_s_list) ):
    if mut_rate not in mut_rate_to_fraction_1_values.keys():
        mut_rate_to_fraction_1_values[mut_rate] = []
    try:
        fraction_1 = 4.0 * s * f_A_mut / ( 4.0 * s + 3.0 * d )
        mut_rate_to_fraction_1_values[mut_rate].append( fraction_1 )
    except:
        break

print('fraction_1')
L = 10000
f_A = 0.2 * L
for mut_rate in mut_rate_to_fraction_1_values.keys():
    p_s = mut_rate
    p_d = mut_rate
    d = mut_rate
    true_fraction_value = 4.0 * (p_s) * ( f_A * (1 - p_s - p_d) + (L - f_A) * p_s/3.0 + L * d / 4.0 ) / (4.0 * p_s + 3.0 * p_d)
    print( mut_rate, np.mean(mut_rate_to_fraction_1_values[mut_rate]), true_fraction_value, np.mean(mut_rate_to_fraction_1_values[mut_rate])/true_fraction_value )

mut_rate_to_fraction_2_values = {}
for s,d,f_A_mut,L,mut_rate,L2 in list( zip(s_list, d_list, f_A_mut_list, L_list, p_s_list, L2_list) ):
    if mut_rate not in mut_rate_to_fraction_2_values.keys():
        mut_rate_to_fraction_2_values[mut_rate] = []
    try:
        fraction_2 = 1.0 * L2 * s / (4.0 * s + 3.0 * d)
        mut_rate_to_fraction_2_values[mut_rate].append( fraction_2 )
    except:
        break

print('fraction_2')
L = 10000
f_A = 0.2 * L
for mut_rate in mut_rate_to_fraction_1_values.keys():
    p_s = mut_rate
    p_d = mut_rate
    d = mut_rate
    true_fraction_value = p_s * L * (1.0 - p_d + d) / (4.0 * p_s + 3.0 * p_d)
    print( mut_rate, np.mean(mut_rate_to_fraction_2_values[mut_rate]), true_fraction_value, np.mean(mut_rate_to_fraction_2_values[mut_rate])/true_fraction_value )
