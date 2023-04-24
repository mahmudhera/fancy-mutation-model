from mutation_model_simulator import mutation_model
import argparse

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


if __name__ == '__main__':
    str_len, p_s, p_d, d, k = parse_arguments()
    seed = 0

    # generate mutation model, generate random str and mutate it randomly
    mm = mutation_model(seed, str_len, p_s, p_d, d)
    str_orig = mm.generate_random_string()
    mutated_string = mm.mutate_string()

    # generate kmers from these two strings
    kmers_orig_set = string_to_kmers(str_orig, k)
    kmers_mutated_set = string_to_kmers(mutated_string, k)

    # observe three quantities
    # ???? HOW ????

    # solve and find three parameters
