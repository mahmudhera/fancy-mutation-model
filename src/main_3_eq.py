from mutation_model_simulator import mutation_model
import argparse
import pygtrie as trie
import math
from tqdm import tqdm
from scipy.optimize import fsolve


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('length', type=int, help='Length of the original string')
    parser.add_argument('k', type=int, help='Length of the kmers')
    parser.add_argument('p_s', type=float, help='Rate of substitution')
    parser.add_argument('p_d', type=float, help='Rate of deletion')
    parser.add_argument('d', type=float, help='Average length of insertion')

    args = parser.parse_args()
    return args.length, args.p_s, args.p_d, args.d, args.k


def string_to_kmers(input_str, k):
    kmer_list = []
    for i in range(len(input_str)-k+1):
        kmer_list.append( input_str[i : i+k] )
    return kmer_list


def observe_kmers_only_one_deletion(kmers_in_orig_str, kmers_in_mutated_str):
    # build prefix tree with kmers in mutated str -> tree1
    tree1 = trie.CharTrie()
    k = len(kmers_in_orig_str[0])
    for kmer in kmers_in_mutated_str:
        for i in range(1, k+1):
            substring = kmer[:i]
            tree1[substring] = substring
            #print(f'Inserting {substring} in tree1')

    # build prefix tree with reversed kmers in mutated str -> tree2
    tree2 = trie.CharTrie()
    for kmer in kmers_in_mutated_str:
        kmer_reversed = kmer[::-1]
        for i in range(1, k+1):
            substring = kmer_reversed[:i]
            tree2[substring] = substring
            #print(f'Inserting {substring} in tree2')

    num_observed_only_one_deletion = 0
    observed_kmers_only_one_deletion = []
    # for all kmers in orig str:
    for kmer in kmers_in_orig_str:
        print(f'Working with kmer: {kmer}')
        # len1 = find len of longest prefix of kmer in tree1
        longest_pref = tree1.longest_prefix(kmer).key
        print(f'Longest prefix in tree1: {longest_pref}')
        len1 = 0
        if longest_pref is not None:
            len1 = len(longest_pref)

        # len2 = find len of longest prefix of rev_kmer in tree2
        kmer_reversed = kmer[::-1]
        longest_pref = tree2.longest_prefix(kmer_reversed).key
        if longest_pref is not None:
            print(f'Longest suffix in tree2: {longest_pref[::-1]}')
        else:
            print(f'Longest suffix in tree2: {longest_pref}')
        len2 = 0
        if longest_pref is not None:
            len2 = len(longest_pref)

        # if len1 + len2 = k - 1 then mark it as one observed
        if len1 == k-1 or len2 == k-1 or len1 + len2 == k - 1:
            num_observed_only_one_deletion += 1
            observed_kmers_only_one_deletion.append(kmer)

    return num_observed_only_one_deletion, observed_kmers_only_one_deletion, tree1, tree2

if __name__ == '__main__':
    #str_len, p_s, p_d, d, k = parse_arguments()
    seed = 0
    k = 21
    str_len = 10000

    # for multiple times, mutate string randomly, and estimate the parameters
    # repeat for a lot of parameters
    # store results in a  file
    f = open('estimation_records_3_eq_approach.csv', 'w')
    num_runs = 40
    for p_s in tqdm([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], desc='p_s progress'):
        for p_d in tqdm([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], desc='p_d progress'):
            for d in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
                mm = mutation_model(seed, str_len, p_s, p_d, d)
                str_orig = mm.generate_random_string()
                kmers_in_orig_str = string_to_kmers(str_orig, k)
                K1 = len(kmers_in_orig_str)

                for i in range(num_runs):
                    mutated_string, num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion = mm.mutate_string(k)
                    kmers_in_mutated_str = string_to_kmers(mutated_string, k)

                    K2 = len(kmers_in_mutated_str)
                    S, I, D = num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion

                    def eqn(inputs):
                        p_s, p_d, d = inputs
                        v1 = k * K1 * (1 - p_s - p_d) ** (k-1) * p_d * (d+1) ** (1.0-k) - D
                        v2 = (k-1) * K1 * (1 - p_s - p_d) ** k * d * (d+1) ** (1.0-k) - I
                        v3 = k * K1 * (1 - p_s - p_d) ** (k-1) * p_s * (d+1) ** (1.0-k) - S
                        return [v1, v2, v3]

                    # solve here
                    p_s_est, p_d_est, d_est = fsolve(eqn, (0,0,0))
                    #print(p_s_est, p_d_est, d_est)
                    f.write( f'{p_s} {p_d} {d} {p_s_est} {p_d_est} {d_est}\n' )
    f.close()

    # generate kmers from these two strings
    # observe three quantities
    # solve and find three parameters

    # TEST
    '''
    str_orig       = 'ACGTGCAGCT'
    mutated_string = 'ACGTGAGCT'
    k              =  5


    kmers_in_orig_str = string_to_kmers(str_orig, k)
    kmers_in_mutated_str = string_to_kmers(mutated_string, k)


    num_observed_only_one_deletion, observed_kmers_only_one_deletion, t1, t2 = observe_kmers_only_one_deletion(kmers_in_orig_str, kmers_in_mutated_str)
    print(num_observed_only_one_deletion)
    print(observed_kmers_only_one_deletion)
    '''
