import numpy as np
import argparse

alphabet = ['A', 'C', 'G', 'T']
substitute_dic = {  'A':['C', 'G', 'T'],
                    'C':['A', 'G', 'T'],
                    'G':['C', 'A', 'T'],
                    'T':['C', 'G', 'A']}

class mutation_model:
    def __init__(self, seed, orig_length, p_s, p_d, d):
        self.seed = seed
        self.orig_length = orig_length
        self.p_s = p_s
        self.p_d = p_d
        self.d = d
        np.random.seed(seed)

    def generate_random_string(self):
        len = self.orig_length
        self.orig_string = ''.join( np.random.choice(alphabet, size=len) )
        return self.orig_string

    def mutate_string(self):
        new_string = list(self.orig_string)
        actions = ['sbst', 'delt', 'stay']

        actions_chosen = np.random.choice( actions, size=self.orig_length, p=[self.p_s, self.p_d, 1.0-self.p_s-self.p_d] )
        insert_lengths = np.random.geometric(1.0/(1+self.d), size=self.orig_length) - 1
        strings_to_insert = [ ''.join( np.random.choice(alphabet, size=len) ) for len in insert_lengths ]

        for i in range( self.orig_length ):
            if actions_chosen[i] == 'sbst':
                new_string[i] = np.random.choice( substitute_dic[new_string[i]] )
            if actions_chosen[i] == 'delt':
                new_string[i] = ''

        final_str_list = []
        for i in range( self.orig_length ):
            final_str_list.append( new_string[i] )
            final_str_list.append( strings_to_insert[i] )

        return ''.join(final_str_list)

def parse_args():
    parser = argparse.ArgumentParser(description='A fancy mutation model supporting insertion, deletion and substitution.')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--len', type=int, default=10000, help='length of the original string')
    parser.add_argument('--ps', type=float, default=0.1, help='substitution probability')
    parser.add_argument('--pd', type=float, default=0.1, help='deletion probability')
    parser.add_argument('--d', type=float, default=1.0, help='average length to insert')
    parser.add_argument('--num_sim', type=int, default=100, help='num of simulations to run')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    seed = args.seed
    orig_len = args.len
    p_s = args.ps
    p_d = args.pd
    d = args.d
    num_simulations = args.num_sim

    mm = mutation_model(seed, orig_len, p_s, p_d, d)

    print('The original string is:')
    print( mm.generate_random_string() )

    print('Now simulating....')
    lengths = []
    for i in range(num_simulations):
        lengths.append( len(mm.mutate_string()) )

    print('Average length of the mutated string:')
    print( 1.0*sum(lengths)/len(lengths) )

    print('By formula:')
    exp_len = orig_len * (1.0 - p_d + d)
    print(exp_len)
