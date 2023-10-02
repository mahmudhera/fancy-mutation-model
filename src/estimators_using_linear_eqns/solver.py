import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

range_low = 20
range_high = 51

if __name__ == '__main__':
    df = pd.read_csv('observations_more_freq_nts', delimiter=' ')

    ps_est_data = {}
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high) ]:
        ps_est_data[f_A] = {}
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            ps_est_data[f_A][f_C] = []

    print(ps_est_data.keys())

    ps_true = 0.01
    df = df[ df['p_s']==ps_true ]
    df = df[ df['k']==21 ]
    df = df[ df['p_d']==0.01 ]
    df = df[ df['d']==0.01 ]

    for index, row in df.iterrows():
        mut_freqs = [row['f_A_mut'], row['f_C_mut'], row['f_G_mut'], row['f_T_mut']]
        orig_freqs = [row['f_A_orig'], row['f_C_orig'], row['f_G_orig'], row['f_T_orig']]

        f_A = -1
        f_A_mut = -1
        for i in range(4):
            if orig_freqs[i] > f_A:
                f_A = orig_freqs[i]
                f_A_mut = mut_freqs[i]

        L = row['L']
        L2 = row['L2']
        D = row['E_D']
        S = row['E_S']
        f_A_true = row['freq_A']
        f_C_true = row['freq_C']
        try:
            p_s_est = (f_A_mut - f_A + L/4.0 - L2/4.0) / ( (L - 4 * f_A) * (1.0/3.0 + D/(4.0 * S) ) )
            ps_est_data[f_A_true][f_C_true].append(p_s_est)
        except:
            p_s_est = None

    for f_A in [ 1.0*i/100 for i in range(range_low, range_high) ]:
        pdf_filename = f'test_plots/plot_f_a_{f_A}.pdf'
        Z1 = []
        Z2 = []
        Z3 = []
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            Z1.append( np.mean(ps_est_data[f_A][f_C]) )
            Z2.append( np.mean(ps_est_data[f_A][f_C]) + np.var(ps_est_data[f_A][f_C])**0.5 )
            Z3.append( np.mean(ps_est_data[f_A][f_C]) - np.var(ps_est_data[f_A][f_C])**0.5 )
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z1)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z2)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z3)
        plt.savefig(pdf_filename)
        plt.clf()
