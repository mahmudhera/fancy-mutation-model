import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

range_low = 20
range_high = 51

def solve_using_polynomials(S, D, I, k, K1, S_smaller):
    k = int(k);
    S_norm = 1.0 * S / (K1 * k)
    D_norm = 1.0 * D / (K1 * k)
    I_norm = 1.0 * I / (K1 * k - K1)

    coeffs = [0 for i in range(k+1)]
    coeffs[0] = (S_norm + D_norm + I_norm) * D_norm**k
    coeffs[-1] = 1
    coeffs[-2] = -D_norm

    roots = np.polynomial.polynomial.polyroots(coeffs)
    p_d_ests = (D_norm - roots)/(S_norm + D_norm + I_norm)
    p_d_ests = [np.real(p_d_est) for p_d_est in p_d_ests if not np.iscomplex(p_d_est)]
    p_d_ests.sort()

    d_ests = [ (D_norm - (S_norm + D_norm) * p_d_est)/(D_norm - (S_norm + D_norm + I_norm) * p_d_est) - 1.0 for p_d_est in p_d_ests ]
    p_s_ests = [ (S_norm * p_d_est)/(D_norm) for p_d_est in p_d_ests ]

    all_solutions = list( zip(p_s_ests, p_d_ests, d_ests) )
    solution_ratio = 0.0
    solution = (None, None, None)
    for p_s_est, p_d_est, d_est in all_solutions:
        if p_s_est < 0 or p_d_est < 0 or d_est < 0:
            continue
        k2 = int(k/2)
        S_smaller_est = K1*(k2)*(1 - p_s_est - p_d_est) ** (k2-1) * p_s_est * (d_est + 1.0)**(-k2+1)
        candidate_ratio = 1.0 * S_smaller/S_smaller_est
        if candidate_ratio > 1:
            candidate_ratio = 1.0/candidate_ratio
        if candidate_ratio < 0.7:
            continue
        if candidate_ratio > solution_ratio:
            solution_ratio = candidate_ratio
            solution = (p_s_est, p_d_est, d_est)

    return solution[0]

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

    ps_est_data_linear = {}
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
        ps_est_data_linear[f_A] = {}
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
            ps_est_data_linear[f_A][f_C] = []
    ps_est_data_halfkmers = {}
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
        ps_est_data_halfkmers[f_A] = {}
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
            ps_est_data_halfkmers[f_A][f_C] = []
    ps_est_data_poly = {}
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
        ps_est_data_poly[f_A] = {}
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high+1) ]:
            ps_est_data_poly[f_A][f_C] = []

    ps_true = 0.01
    df = df[ df['p_s']==ps_true ]
    df = df[ df['k']==31 ]
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
        I = row['E_I']
        f_A_true = row['freq_A']
        f_C_true = row['freq_C']
        S_k_half = row['S_sp']
        S_k = S
        k = row['k']
        K1 = row['K1']
        S_smaller = S_k_half

        ps_est_data_linear[f_A_true][f_C_true].append( solve_using_linear_equations(f_A, f_A_mut, L, L2, D, S) )
        ps_est_data_halfkmers[f_A_true][f_C_true].append( solve_using_half_kmers(k, S_k_half, S_k, K1) )
        ps_est = solve_using_polynomials(S, D, I, k, K1, S_smaller)
        if ps_est is not None:
            ps_est_data_poly[f_A_true][f_C_true].append( ps_est )



    fig, ax = plt.subplots()
    for f_A in [ 1.0*i/100 for i in range(range_low, range_high) ]:
        pdf_filename = f'test_plots/plot_f_a_{f_A}.pdf'
        Z1 = []
        Z2 = []
        Z3 = []
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            Z1.append( np.mean(ps_est_data_linear[f_A][f_C]) )
            Z2.append( np.mean(ps_est_data_linear[f_A][f_C]) + np.var(ps_est_data_linear[f_A][f_C])**0.5 )
            Z3.append( np.mean(ps_est_data_linear[f_A][f_C]) - np.var(ps_est_data_linear[f_A][f_C])**0.5 )

        plt.fill_between([ 1.0*i/100 for i in range(range_low, range_high) ], Z3, Z2, color='blue', alpha=0.2)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z1, color='blue', label='Linear')

        Z1 = []
        Z2 = []
        Z3 = []
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            Z1.append( np.mean(ps_est_data_halfkmers[f_A][f_C]) )
            Z2.append( np.mean(ps_est_data_halfkmers[f_A][f_C]) + np.var(ps_est_data_halfkmers[f_A][f_C])**0.5 )
            Z3.append( np.mean(ps_est_data_halfkmers[f_A][f_C]) - np.var(ps_est_data_halfkmers[f_A][f_C])**0.5 )

        plt.fill_between([ 1.0*i/100 for i in range(range_low, range_high) ], Z3, Z2, color='red', alpha=0.2)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z1, color='red', label='Half-k-mers')

        Z1 = []
        Z2 = []
        Z3 = []
        for f_C in [ 1.0*i/100 for i in range(range_low, range_high) ]:
            Z1.append( np.mean(ps_est_data_poly[f_A][f_C]) )
            Z2.append( np.mean(ps_est_data_poly[f_A][f_C]) + np.var(ps_est_data_poly[f_A][f_C])**0.5 )
            Z3.append( np.mean(ps_est_data_poly[f_A][f_C]) - np.var(ps_est_data_poly[f_A][f_C])**0.5 )

        plt.fill_between([ 1.0*i/100 for i in range(range_low, range_high) ], Z3, Z2, color='green', alpha=0.2)
        plt.plot([ 1.0*i/100 for i in range(range_low, range_high) ], Z1, color='green', label='Polynomial')

        plt.plot( [1.0*range_low/100, 1.0*(range_high-1)/100], [ps_true, ps_true], color='grey', linestyle='--' )

        plt.legend()
        plt.xlabel('Frequency of Cs in orig string')
        plt.ylabel('p_s_est')
        #plt.ylim(0.09, 0.11)
        plt.xlim(0.2, 0.49)
        plt.title(f'Frequency of As set to {f_A}')
        plt.savefig(pdf_filename)
        plt.clf()
