from mutation_model_simulator import mutation_model
import argparse
import math
from tqdm import tqdm
from scipy.optimize import fsolve

def string_to_kmers(input_str, k):
    kmer_list = []
    for i in range(len(input_str)-k+1):
        kmer_list.append( input_str[i : i+k] )
    return kmer_list


mut_rates = [0.01]
frequency_distr = [0.4, 0.2, 0.3, 0.1]

if __name__ == '__main__':
    obs_filename = f'observations.csv'
    f = open(obs_filename, 'w')
    f.write('p_s p_d d k L L2 K1 K2 E_S E_D E_I f_A_orig f_A_mut\n')
    f.flush()

    num_runs = 1
    seed_to_generate = 0

    for k in [21, 31, 41]:
        for str_len in [10000]:
            for p_s in tqdm(mut_rates, desc='p_s progress'):
                for p_d in tqdm(mut_rates, desc='p_d progress'):
                    for d in mut_rates:
                        mm = mutation_model(seed_to_generate, str_len, p_s, p_d, d)
                        str_orig = mm.generate_random_string(frequency_distr)
                        kmers_in_orig_str = string_to_kmers(str_orig, k)
                        K1 = len(kmers_in_orig_str)
                        f_A = str_orig.count('A')
                        L = str_len
                        for i in range(num_runs):
                            mutated_string, num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion = mm.mutate_string(k)
                            kmers_in_mutated_str = string_to_kmers(mutated_string, k)
                            K2 = len(kmers_in_mutated_str)
                            L2 = len(mutated_string)
                            S_obs, I_obs, D_obs = num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion
                            f_A_dash = mutated_string.count('A')
                            f.write( f'{p_s} {p_d} {d} {k} {L} {L2} {K1} {K2} {S_obs} {D_obs} {I_obs} {f_A} {f_A_dash}\n' )
                            f.flush()
    f.close()
