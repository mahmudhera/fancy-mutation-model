import pandas as pd
import numpy as np

ps_filter = 0.1
pd_filter = 0.1
d_filter = 0.1
L_filter = 1000000

if __name__ == '__main__':
    df = pd.read_csv('observations.csv', delimiter=' ')

    df = df[ df['p_s']==ps_filter ]
    df = df[ df['p_d']==pd_filter ]
    df = df[ df['d']==d_filter ]
    #df = df[  df['L']==L_filter ]

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

    for p_s, p_d, d, k, L, L2, K1, K2, S, D, I, f_A, f_A_mut in list( zip(p_s_list, p_d_list, d_list, k_list, L_list, L2_list, K1_list, K2_list, E_S_list, E_D_list, E_I_list, f_A_list, f_A_mu_list) ):
        try:
            p_s_est = (f_A_mut - f_A + L/4.0 - L2/4.0) / ( (L - 4 * f_A) * (1.0/3.0 + D/(4.0 * S) ) )
        except ZeroDivisionError:
            p_s_est = None
        print(f'Solution using linear stuff: p_s={p_s_est}, True p_s={p_s}')
        est_vals.append(p_s_est)

    print(np.mean(est_vals), np.var(est_vals), np.mean(est_vals)/p_s)
