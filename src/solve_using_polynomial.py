import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.read_csv('observations.csv', delimiter=' ', header=None)
    pss = df[0].tolist()
    pds = df[1].tolist()
    ds  = df[2].tolist()
    Ss  = df[3].tolist()
    Ds  = df[4].tolist()
    Is  = df[5].tolist()
    K1s = df[6].tolist()
    ks  = df[7].tolist()
    S_smallers = df[8].tolist()

    for ps, pd, d, S, D, I, K1, k, S_smaller in zip(pss, pds, ds, Ss, Ds, Is, K1s, ks, S_smallers):
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

        print(ps, pd, d, solution[0], solution[1], solution[2])
