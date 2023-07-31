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

def count_num_kmers_with_starting_character(kmer_list, start_char):
    count = 0
    for kmer in kmer_list:
        if start_char == kmer[0]:
            count += 1
    return count

if __name__ == '__main__':
    #str_len, p_s, p_d, d, k = parse_arguments()
    seed = 0
    k = 21
    str_len = 10000
    nucleotide_frequency = [0.1, 0.2, 0.3, 0.4]

    # for multiple times, mutate string randomly, and estimate the parameters
    # repeat for a lot of parameters
    # store results in a  file
    #f = open('observations_with_k2_and_special_khalved.csv', 'w')
    f = open('observations_with_k1_k2_ka_kc.csv', 'w')

    num_runs = 20
    for p_s in tqdm([0.01*(i+1) for i in range(5)], desc='p_s progress'):
        for p_d in tqdm([0.01*(i+1) for i in range(5)], desc='p_d progress'):
            for d in [0.01*(i+1) for i in range(5)]:
                mm = mutation_model(seed, str_len, p_s, p_d, d)
                str_orig = mm.generate_random_string(nucleotide_frequency)
                kmers_in_orig_str = string_to_kmers(str_orig, k)
                K1 = len(kmers_in_orig_str)

                KA = count_num_kmers_with_starting_character(kmers_in_orig_str, 'A')
                KC = count_num_kmers_with_starting_character(kmers_in_orig_str, 'C')

                for i in range(num_runs):
                    mutated_string, num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion, num_k_1_mers_single_substitution, len_mutated_str = mm.mutate_string(k)
                    kmers_in_mutated_str = string_to_kmers(mutated_string, k)

                    K2 = len(kmers_in_mutated_str)
                    YA = count_num_kmers_with_starting_character(kmers_in_mutated_str, 'A')
                    YC = count_num_kmers_with_starting_character(kmers_in_mutated_str, 'C')
                    S, I, D, S_smaller = num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion, num_k_1_mers_single_substitution

                    f.write( f'{p_s} {p_d} {d} {K1} {KA} {KC} {K2} {YA} {YC}\n' )
                    f.flush()
    f.close()
