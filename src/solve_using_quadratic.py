import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.read_csv('observations_with_k2.csv', delimiter=' ', header=None)

    print(df)

    #ps_test = 0.03
    #pd_test = 0.01
    #d_test = 0.01

    #df = df[ df[0] == ps_test ]
    #df = df[ df[1] == pd_test ]
    #df = df[ df[2] == d_test ]

    pss = df[0].tolist()
    pds = df[1].tolist()
    ds  = df[2].tolist()
    Ss  = df[3].tolist()
    Ds  = df[4].tolist()
    Is  = df[5].tolist()
    K1s = df[6].tolist()
    K2s = df[7].tolist()
    ks  = df[8].tolist()
    S_smallers = df[9].tolist()

    for ps, pd, d, S, D, I, K1, K2, k, S_smaller in zip(pss, pds, ds, Ss, Ds, Is, K1s, K2s, ks, S_smallers):
        S_norm = 1.0 * S / (K1 * k)
        D_norm = 1.0 * D / (K1 * k)
        I_norm = 1.0 * I / (K1 * k - K1)
        observed_pd_diff_d = 1.0 - 1.0 * (K2 + k - 1) / (K1 + k - 1)

        a = S_norm + D_norm + I_norm
        b = I_norm - observed_pd_diff_d * (S_norm + D_norm + I_norm) - D_norm
        c = observed_pd_diff_d * D_norm

        pd_est = (-b + (b*b - 4*a*c)**0.5)/(2*a)
        ps_est = 1.0 * S * pd_est / D
        d_est = pd_est - observed_pd_diff_d

        if np.iscomplex(pd_est):
            pd_est, ps_est, d_est = None, None, None

        print(ps, pd, d, ps_est, pd_est, d_est)
