import pandas as pd
import numpy as np

if __name__ == '__main__':
    df = pd.read_csv('observations_with_k1_k2_ka_kc.csv', delimiter=' ', header=None)

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
    K1s  = df[3].tolist()
    KAs  = df[4].tolist()
    KCs  = df[5].tolist()
    K2s = df[6].tolist()
    YAs = df[7].tolist()
    YCs  = df[8].tolist()

    k = 21

    for ps, pd, d, K1, KA, KC, K2, YA, YC in zip(pss, pds, ds, K1s, KAs, KCs, K2s, YAs, YCs):
        pd_est = - 1.0 * (K2*KA - K2*KC - K1*YA + 4*KC*YA + K1*YC - 4*KA*YC) / ( (-1 + k)*(KA - KC) )
        ps_est = 1.0 * ( -3*(KC*(-1 + k + K2 - 4*YA) + (-1 + k + K1)*(YA - YC)) + 3*KA*(-1 + k + K2 - 4*YC) ) / (4*(-1 + k)*(KA - KC))
        d_est = ( (K1*(KC*(-1 + k + K2 - 4*YA) - KA*(-1 + k + K2 - 4*YC) + (-1 + k)*(YA - YC)) + K1**2 * (YA - YC) - 4*(-1 + k)*(KC*YA - KA*YC)) ) / ((-1 + k)*(-1 + k + K1)*(KA - KC))

        print(ps, pd, d, ps_est, pd_est, d_est)
