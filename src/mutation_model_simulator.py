import numpy as np

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

if __name__ == '__main__':
    mm = mutation_model(2, 10, 0.2, 0.2, 1)
    print( mm.generate_random_string() )
    print( mm.mutate_string() )
