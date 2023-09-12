import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.read_csv('observations.csv', delimiter=' ')

    p_s_list = df['p_s'].tolist()
    p_d_list = df['p_d'].tolist()
    d_list = df['d'].tolist()
    k_list = df['k'].tolist()
    L_list = df['L'].tolist()
    L2_list = df['L2'].tolist()
    K1_list = df['K1'].tolist()
    K2_list = df['K2'].tolist()
    E_S_list = df['E_S'].tolist()
    E_D_list = df['E_D'].tolist()
    E_I_list = df['E_I'].tolist()
    f_A_list = df['f_A_orig'].tolist()
    f_A_mu_list = df['f_A_mut'].tolist()

    est_vals = []
    for p_s, p_d, d, k, L, L2, K1, K2, E_S, E_D, E_I, f_A, f_A_mut in list( zip(p_s_list, p_d_list, d_list, k_list, L_list, L2_list, K1_list, K2_list, E_S_list, E_D_list, E_I_list, f_A_list, f_A_mu_list) ):
        p_s_est = (f_A_mut - f_A + L/4.0 - L2/4.0) / ( (L - 4 * f_A) * (1.0/3.0 + E_D/(4.0 * E_S) ) )
        print(p_s, p_s_est)
        print(L * (1 - p_d + d), L2)
        print(f_A * (1 - p_s - p_d) + (L - f_A)*p_s/3.0 + L * d / 4.0, f_A_mut)
        print(E_S/E_D, p_s/p_d)
        est_vals.append(p_s_est)
        if len(est_vals) > 18:
            break

    print(est_vals)
    print(np.mean(est_vals))
