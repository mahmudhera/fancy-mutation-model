import pandas as pd
from tqdm import tqdm
import numpy as np

def solve_using_polynomials(S, D, I, k, K1, S_smaller):
    k = int(k)
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
        k2 = int( (k+1)/2 )
        S_smaller_est = K1*(k2)*(1 - p_s_est - p_d_est) ** (k2-1) * p_s_est * (d_est + 1.0)**(-k2+1)
        candidate_ratio = 1.0 * S_smaller/S_smaller_est
        if candidate_ratio > 1:
            candidate_ratio = 1.0/candidate_ratio
        if candidate_ratio < 0.7:
            continue
        if candidate_ratio > solution_ratio:
            solution_ratio = candidate_ratio
            solution = (p_s_est, p_d_est, d_est)

    return solution[0], solution[1], solution[2]

def solve_using_half_kmers(k, S_k_half, S_k, K1, D_k_half, D_k, L2, L):
    #print(k, S_k_half, S_k, K1)
    #exit(-1)
    try:
        p_s_est = 4.0 * k * S_k_half**2 / ((k+1)**2 * S_k * K1)
        p_d_est = 4.0 * k * D_k_half**2 / ((k+1)**2 * D_k * K1)
        d_est   = 1.0*L2/L - 1.0 + p_d_est
    except:
        p_s_est = None
        p_d_est = None
        d_est = None
    return p_s_est, p_d_est, d_est

def solve_using_linear_equations1(f_A, f_A_mut, L, L2, D, S):
    try:
        p_s_est = (f_A_mut - f_A + L/4.0 - L2/4.0) / ( (L - 4 * f_A) * (1.0/3.0 + D/(4.0 * S) ) )
        p_d_est = 1.0 * D * p_s_est / S
        d_est   = 1.0*L2/L - 1.0 + p_d_est
    except:
        p_s_est = None
        p_d_est = None
        d_est = None
    return p_s_est, p_d_est, d_est

def solve_using_linear_equations2(f_A, f_A_mut, L, L2, I, N_shared, k):
    try:
        d_est = I / (N_shared * (k-1) - I)
        p_d_est = 1 + d_est - 1.0 * L2 / L
        a_constant = (L - 4 * f_A) / 3.0
        p_s_est = ( f_A_mut - f_A + f_A * p_d_est - L * d_est / 4.0 ) / a_constant
    except:
        p_s_est = None
        p_d_est = None
        d_est = None
    return p_s_est, p_d_est, d_est

if __name__ == '__main__':
    file_path = 'observations.gz'
    df = pd.read_csv(file_path, compression='gzip', header=0, sep=' ')
    df_out = pd.DataFrame(df)
    df_out['p_s_est_halfkmer'] = None
    df_out['p_s_est_poly'] = None
    df_out['p_s_est_linear1'] = None
    df_out['p_s_est_linear2'] = None
    df_out['p_d_est_halfkmer'] = None
    df_out['p_d_est_poly'] = None
    df_out['p_d_est_linear1'] = None
    df_out['p_d_est_linear2'] = None
    df_out['d_est_halfkmer'] = None
    df_out['d_est_poly'] = None
    df_out['d_est_linear1'] = None
    df_out['d_est_linear2'] = None

    # following are the header columns
    # freq_A freq_C freq_G freq_T p_s p_d d k L L2 K1 K2 E_S E_D E_I f_A_orig f_A_mut f_C_orig f_C_mut f_G_orig f_G_mut f_T_orig f_T_mut S_sp D_sp I_sp N_shared
    # iterate over all rows in the dataframe
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        # get the values from the row
        freq_A = row['freq_A']
        freq_C = row['freq_C']
        freq_G = row['freq_G']
        freq_T = row['freq_T']
        p_s = row['p_s']
        p_d = row['p_d']
        d = row['d']
        k = row['k']
        L = row['L']
        L2 = row['L2']
        K1 = row['K1']
        K2 = row['K2']
        E_S = row['E_S']
        E_D = row['E_D']
        E_I = row['E_I']
        f_A_orig = row['f_A_orig']
        f_A_mut = row['f_A_mut']
        f_C_orig = row['f_C_orig']
        f_C_mut = row['f_C_mut']
        f_G_orig = row['f_G_orig']
        f_G_mut = row['f_G_mut']
        f_T_orig = row['f_T_orig']
        f_T_mut = row['f_T_mut']
        S_sp = row['S_sp']
        D_sp = row['D_sp']
        I_sp = row['I_sp']
        N_shared = row['N_shared']

        # solve using half-kmers
        S_k_half = S_sp
        S_k = E_S 
        D_k_half = D_sp
        D_k = E_D
        p_s_est, p_d_est, d_est = solve_using_half_kmers(k, S_k_half, S_k, K1, D_k_half, D_k, L2, L)
        df_out.at[index, 'p_s_est_halfkmer'] = p_s_est
        df_out.at[index, 'p_d_est_halfkmer'] = p_d_est
        df_out.at[index, 'd_est_halfkmer'] = d_est

        # solve using polynomials
        p_s_est, p_d_est, d_est = solve_using_polynomials(E_S, E_D, E_I, k, K1, S_sp)
        df_out.at[index, 'p_s_est_poly'] = p_s_est
        df_out.at[index, 'p_d_est_poly'] = p_d_est
        df_out.at[index, 'd_est_poly'] = d_est

        # solve using linear equations
        p_s_est, p_d_est, d_est = solve_using_linear_equations1(f_A_orig, f_A_mut, L, L2, E_D, E_S)
        df_out.at[index, 'p_s_est_linear1'] = p_s_est
        df_out.at[index, 'p_d_est_linear1'] = p_d_est
        df_out.at[index, 'd_est_linear1'] = d_est
        
        # solve using new linear equations
        p_s_est, p_d_est, d_est = solve_using_linear_equations2(f_A_orig, f_A_mut, L, L2, E_I, N_shared, k)
        df_out.at[index, 'p_s_est_linear2'] = p_s_est
        df_out.at[index, 'p_d_est_linear2'] = p_d_est
        df_out.at[index, 'd_est_linear2'] = d_est

    # only keep the columns freq_A, freq_C, p_s, p_d, d, k, L, p_s_est_halfkmer, p_s_est_poly, p_s_est_linear1, p_s_est_linear2, p_d_est_halfkmer, p_d_est_poly, p_d_est_linear1, p_d_est_linear2, d_est_halfkmer, d_est_poly, d_est_linear1, d_est_linear2
    df_out = df_out[['freq_A', 'freq_C', 'p_s', 'p_d', 'd', 'k', 'L', 'p_s_est_halfkmer', 'p_s_est_poly', 'p_s_est_linear1', 'p_s_est_linear2', 'p_d_est_halfkmer', 'p_d_est_poly', 'p_d_est_linear1', 'p_d_est_linear2', 'd_est_halfkmer', 'd_est_poly', 'd_est_linear1', 'd_est_linear2']]
    df_out.to_csv('observations_with_estimates_of_d.gz', compression='gzip', sep=' ', index=False)