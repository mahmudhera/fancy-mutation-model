from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

width = 0.001
df = pd.read_csv('results_using_polynomials', delimiter=' ', header=None)
df.columns = ['p_s', 'p_d', 'd', 'p_s_est', 'p_d_est', 'd_est']
#print(df)

file_serializer = 0
for p_s in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    for p_d in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
        f_name = f'plots_polynomials/plots_d_est/{file_serializer:02d}_ps_{p_s}_pd_{p_d}.pdf'
        file_serializer += 1

        df2 = df[ df['p_s'] == p_s ]
        df2 = df2[ df2['p_d'] == p_d ]

        ds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        d_ests = []
        d_est_stdvs = []

        for d in ds:
            d_est_list = df2[ df2['d'] == d ]['d_est'].tolist()
            d_est_list = [ float(d_est) for d_est in d_est_list if str(d_est) != 'None']
            try:
                d_est_avg = np.mean(d_est_list)
                d_est_var = np.var(d_est_list)
            except:
                print(d_est_list)
                exit(-1)
            d_ests.append(d_est_avg)
            d_est_stdvs.append(d_est_var ** 0.5)
        plt.scatter(ds, d_ests, color='purple')
        for d, d_est, d_est_stddev in list(zip(ds, d_ests, d_est_stdvs)):
            plt.plot( [d, d], [d_est-d_est_stddev, d_est+d_est_stddev], color='purple' )
            plt.plot( [d-width, d+width], [d_est-d_est_stddev, d_est-d_est_stddev], color='purple' )
            plt.plot( [d-width, d+width], [d_est+d_est_stddev, d_est+d_est_stddev], color='purple' )

        plt.plot( [0.01, 0.1], [0.01, 0.1], linestyle='--', color='grey' )
        plt.title(f'p_s = {p_s}, p_d = {p_d}')
        plt.savefig(f_name)
        #plt.show()
        plt.clf()

print('Checkpoint 1')
#exit(0)

file_serializer = 0
for d in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    for p_d in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
        f_name = f'plots_polynomials/plots_ps_est/{file_serializer:02d}_ps_{p_s}_pd_{p_d}.pdf'
        file_serializer += 1

        df2 = df[ df['d'] == d ]
        df2 = df2[ df2['p_d'] == p_d ]

        p_ss = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        p_s_ests = []
        p_s_est_stdvs = []

        for p_s in p_ss:
            p_s_est_list = df2[ df2['p_s'] == p_s ]['p_s_est'].tolist()
            p_s_est_list = [ float(p_s_est) for p_s_est in p_s_est_list if str(p_s_est)!='None' ]
            p_s_est_avg = np.mean(p_s_est_list)
            p_s_est_var = np.var(p_s_est_list)
            p_s_ests.append(p_s_est_avg)
            p_s_est_stdvs.append(p_s_est_var ** 0.5)
        plt.scatter(p_ss, p_s_ests, color='purple')
        for p_s, p_s_est, p_s_est_stddev in list(zip(p_ss, p_s_ests, p_s_est_stdvs)):
            plt.plot( [p_s, p_s], [p_s_est-p_s_est_stddev, p_s_est+p_s_est_stddev], color='purple' )
            plt.plot( [p_s-width, p_s+width], [p_s_est-p_s_est_stddev, p_s_est-p_s_est_stddev], color='purple' )
            plt.plot( [p_s-width, p_s+width], [p_s_est+p_s_est_stddev, p_s_est+p_s_est_stddev], color='purple' )

        plt.plot( [0.01, 0.1], [0.01, 0.1], linestyle='--', color='grey' )
        plt.title(f'd = {d}, p_d = {p_d}')
        plt.savefig(f_name)
        plt.clf()

print('Checkpoint 2')
file_serializer = 0
for d in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    for p_s in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
        f_name = f'plots_polynomials/plots_pd_est/{file_serializer:02d}_ps_{p_s}_pd_{p_d}.pdf'
        file_serializer += 1

        df2 = df[ df['d'] == d ]
        df2 = df2[ df2['p_s'] == p_s ]

        p_ds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
        p_d_ests = []
        p_d_est_stdvs = []

        for p_d in p_ds:
            p_d_est_list = df2[ df2['p_d'] == p_d ]['p_d_est'].tolist()
            p_d_est_list = [ float(p_d_est) for p_d_est in p_d_est_list if str(p_d_est)!='None' ]
            p_d_est_avg = np.mean(p_d_est_list)
            p_d_est_var = np.var(p_d_est_list)
            p_d_ests.append(p_d_est_avg)
            p_d_est_stdvs.append(p_d_est_var ** 0.5)
        plt.scatter(p_ds, p_d_ests, color='purple')
        for p_d, p_d_est, p_d_est_stddev in list(zip(p_ds, p_d_ests, p_d_est_stdvs)):
            plt.plot( [p_d, p_d], [p_d_est-p_d_est_stddev, p_d_est+p_d_est_stddev], color='purple' )
            plt.plot( [p_d-width, p_d+width], [p_d_est-p_d_est_stddev, p_d_est-p_d_est_stddev], color='purple' )
            plt.plot( [p_d-width, p_d+width], [p_d_est+p_d_est_stddev, p_d_est+p_d_est_stddev], color='purple' )

        plt.plot( [0.01, 0.1], [0.01, 0.1], linestyle='--', color='grey' )
        plt.title( f'd = {d}, p_s = {p_s}' )
        plt.savefig(f_name)
        plt.clf()
