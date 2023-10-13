import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

range_low = 20
range_high = 51

def solve_using_half_kmers(k, S_k_half, S_k, K1):
    #print(k, S_k_half, S_k, K1)
    #exit(-1)
    try:
        p_s_est = 4.0 * k * S_k_half**2 / ((k+1)**2 * S_k * K1)
    except:
        p_s_est = None
    return p_s_est

def solve_using_linear_equations(f_A, f_A_mut, L, L2, D, S):
    try:
        p_s_est = (f_A_mut - f_A + L/4.0 - L2/4.0) / ( (L - 4 * f_A) * (1.0/3.0 + D/(4.0 * S) ) )
    except:
        p_s_est = None
    return p_s_est

if __name__ == '__main__':
    df = pd.read_csv('observations_more_freq_nts', delimiter=' ')

    ps_est_data = {}
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
        ps_est_data[f_A] = {}
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
            ps_est_data[f_A][f_C] = []

    print(ps_est_data.keys())

    ps_true = 0.1
    df = df[ df['p_s']==ps_true ]
    df = df[ df['k']==21 ]
    df = df[ df['p_d']==0.1 ]
    df = df[ df['d']==0.1 ]

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
        S_k_half = row['S_sp']
        S_k = S
        k = row['k']
        K1 = row['K1']

        #ps_est_data[f_A_true][f_C_true].append( solve_using_linear_equations(f_A, f_A_mut, L, L2, D, S) )
        ps_est_data[f_A_true][f_C_true].append( solve_using_half_kmers(k, S_k_half, S_k, K1) )


    for f_A in [ 1.0*i/100 for i in range(range_low, range_high) ]:
        pdf_filename = f'test_plots/plot_f_a_{f_A}.pdf'
        Z1 = []
        Z2 = []
        Z3 = []
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            Z1.append( np.mean(ps_est_data[f_A][f_C]) )
            Z2.append( np.mean(ps_est_data[f_A][f_C]) + np.var(ps_est_data[f_A][f_C])**0.5 )
            Z3.append( np.mean(ps_est_data[f_A][f_C]) - np.var(ps_est_data[f_A][f_C])**0.5 )

        fig, ax = plt.subplots()
        ax.fill_between([ 1.0*i/100 for i in range(range_low, range_high) ], Z3, Z2, color='blue', alpha=0.2)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z1, color='blue')
        #plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z2, color='blue', alpha=0.5)
        #plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z3, color='blue', alpha=0.5)
        #plt.plot( [range_low, range_high], [ps_true, ps_true], linestyle='--', color='grey' )
        plt.xlabel('Frequency of Cs in orig string')
        plt.ylabel('p_s_est')
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            for val in ps_est_data[f_A][f_C]:
                plt.scatter([f_C], [val], color='blue')
        #plt.ylim(0.09, 0.11)
        plt.savefig(pdf_filename)
        plt.clf()
