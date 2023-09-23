#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <unordered_map>
#include <tuple>
#include <fstream>

using namespace std;

std::vector<char> alphabet = {'A', 'C', 'G', 'T'};

std::unordered_map<char, std::vector<char>> substitute_dic = {
    {'A', {'C', 'G', 'T'}},
    {'C', {'A', 'G', 'T'}},
    {'G', {'C', 'A', 'T'}},
    {'T', {'C', 'G', 'A'}}
};

class mutation_model {
public:
    mutation_model(int seed, size_t orig_length, double p_s, double p_d, double d)
        : seed(seed), orig_length(orig_length), p_s(p_s), p_d(p_d), d(d) {
        srand(seed);
    }

    string generate_random_string(const vector<double>& frequency = {}) {
        size_t len = orig_length;
        string orig_string;

        orig_string.reserve(len);
        if (frequency.empty()) {
            for (size_t i = 0; i < len; i++) {
                orig_string.push_back(alphabet[std::rand() % alphabet.size()]);
            }
        } else {
            for (size_t i = 0; i < len; i++) {
                double r = static_cast<double>(rand()) / RAND_MAX;
                double cumulative_prob = 0.0;

                for (int j = 0; j < alphabet.size(); j++) {
                    cumulative_prob += frequency[j];
                    if (r < cumulative_prob) {
                        orig_string.push_back(alphabet[j]);
                        break;
                    }
                }
            }
        }

        orig_string.shrink_to_fit();
        this->orig_string = orig_string;
        return orig_string;
    }

    tuple<string, size_t, size_t, size_t, size_t, size_t, size_t> mutate_string(int k = -1) {
        string new_string = string(this->orig_string);
        vector<int> actions_chosen(orig_length);
        vector<int> insert_lengths(orig_length);
        vector<string> strings_to_insert(orig_length);

        // Choose actions for each position in the string
        for (size_t i = 0; i < orig_length; i++) {
            actions_chosen[i] = choose_action();
            insert_lengths[i] = choose_insert_length();
            strings_to_insert[i] = generate_random_string(insert_lengths[i]);
        }

        // Apply chosen actions to mutate the string
        for (size_t i = 0; i < orig_length; i++) {
            if (actions_chosen[i] == 2) {
                char orig_char = new_string[i];
                new_string[i] = substitute_dic[orig_char][std::rand() % 3];
            } else if (actions_chosen[i] == 3) {
                new_string[i] = ' ';
            }
        }

        // Combine the mutated string and inserted substrings
        std::string final_str = combine_strings(new_string, strings_to_insert);

        // Calculate mutation counts if k is provided
        if (k == -1) {
            return std::make_tuple(final_str, 0, 0, 0, 0, 0, 0);
        }

        size_t num_kmers_single_insertion = 0;
        size_t num_kmers_single_deletion = 0;
        size_t num_kmers_single_substitution = 0;

        for (size_t i = 0; i < orig_length - k + 1; i++) {
            int sum1 = 0;
            int sum2 = 0;
            for (int j = i; j < i + k; j++) {
                sum1 += actions_chosen[j];
                sum2 += insert_lengths[j];
            }
            if (sum1 == 3 && sum2 == 0) {
                num_kmers_single_deletion++;
            }
            if (sum1 == 2 && sum2 == 0) {
                num_kmers_single_substitution++;
            }
            if (sum1 == 0 && sum2 == 1) {
                num_kmers_single_insertion++;
            }
        }

        size_t num_kmers_single_insertion_special = 0;
        size_t num_kmers_single_deletion_special = 0;
        size_t num_kmers_single_substitution_special = 0;

        int short_ksize = (int)( (k+1)/2 );

        for (size_t i = 0; i < orig_length - short_ksize + 1; i++) {
            int sum1 = 0;
            int sum2 = 0;
            for (int j = i; j < i + short_ksize; j++) {
                sum1 += actions_chosen[j];
                sum2 += insert_lengths[j];
            }
            if (sum1 == 3 && sum2 == 0) {
                num_kmers_single_deletion_special++;
            }
            if (sum1 == 2 && sum2 == 0) {
                num_kmers_single_substitution_special++;
            }
            if (sum1 == 0 && sum2 == 1) {
                num_kmers_single_insertion_special++;
            }
        }

        return make_tuple(final_str, num_kmers_single_substitution, num_kmers_single_insertion, num_kmers_single_deletion,
            num_kmers_single_substitution_special, num_kmers_single_insertion_special, num_kmers_single_deletion_special);
    }

    void test() {
        cout << this->orig_string << endl;
    }

private:
    int seed;
    size_t orig_length;
    double p_s;
    double p_d;
    double d;
    string orig_string;

    int choose_action() {
        double rand_num = static_cast<double>(rand()) / RAND_MAX;
        if (rand_num < p_s) {
            return 2; // Substitution
        } else if (rand_num < p_s + p_d) {
            return 3; // Deletion
        } else {
            return 0; // Stay (no mutation)
        }
    }

    // Choose an insertion length following a geometric distribution
    int choose_insert_length() {
        double rand_num = static_cast<double>(rand()) / RAND_MAX;
        return static_cast<int>(log(rand_num) / log(d/(1.0+d)) );
    }

    // Generate a random string of a given length
    string generate_random_string(int length) {
        string random_str;
        random_str.reserve(length);
        for (int i = 0; i < length; i++) {
            random_str.push_back(alphabet[rand() % alphabet.size()]);
        }
        return random_str;
    }

    // Combine the mutated string with inserted substrings
    string combine_strings(const string& base, const vector<string>& inserts) {
        string final_str;
        final_str.reserve(base.size() + inserts.size());
        for (int i = 0; i < orig_length; i++) {
            if (base[i] != ' ')
                final_str += base[i];
            final_str += inserts[i];
        }
        return final_str;
    }
};


int main(int argc, char* argv[])
{
    cout << "Usage: ./program_name output_filename" << endl;
    if (argc < 2) {
        cout << "Output filename not given" << endl;
        return 0;
    }
    string filename = string(argv[1]);
    ofstream fout(filename);

    fout << "freq_A freq_C freq_G freq_T p_s p_d d k L L2 K1 K2 E_S E_D E_I f_A_orig f_A_mut f_C_orig f_C_mut f_G_orig f_G_mut f_T_orig f_T_mut S_sp D_sp I_sp" << endl;
    size_t num_runs_each_setting = 20;
    int seed = 0;
    int ksizes[] = {21, 31, 41};
    size_t str_lengths[] = {10000, 100000, 1000000};
    vector<double> mutation_rates = {};
    for (int i = 1; i <= 15; i+=2) {
        mutation_rates.push_back(1.0*i/100);
    }

    for (double freq_A = 0.1; freq_A <= 0.8; freq_A += 0.1) {
        for (double freq_C = 0.1; freq_A <= 0.8; freq_A += 0.1) {
            if (freq_A + freq_C >= 1.0) {
                continue;
            }
            double freq_G = (1.0 - freq_A - freq_C)/2.0;
            double freq_T = (1.0 - freq_A - freq_C)/2.0;
            vector<double> nucloetide_frequencies = {freq_A, freq_C, freq_G, freq_T};

            for (auto str_len : str_lengths) {
                for (auto p_s : mutation_rates) {
                    for (auto p_d : mutation_rates) {
                        for (auto d : mutation_rates) {
                            mutation_model model(seed, str_len, p_s, p_d, d);
                            string orig_string = model.generate_random_string(nucloetide_frequencies);
                            size_t f_A_orig = count(orig_string.begin(), orig_string.end(), 'A');
                            size_t f_C_orig = count(orig_string.begin(), orig_string.end(), 'C');
                            size_t f_G_orig = count(orig_string.begin(), orig_string.end(), 'G');
                            size_t f_T_orig = count(orig_string.begin(), orig_string.end(), 'T');
                            size_t L = str_len;
                            for (auto ksize : ksizes) {
                                for (int i = 0; i < num_runs_each_setting; i++) {
                                    auto mut_res = model.mutate_string(ksize);
                                    size_t L2 = get<0>(mut_res).length();
                                    size_t S = get<1>(mut_res);
                                    size_t I = get<2>(mut_res);
                                    size_t D = get<3>(mut_res);
                                    size_t S_sp = get<4>(mut_res);
                                    size_t I_sp = get<5>(mut_res);
                                    size_t D_sp = get<6>(mut_res);

                                    size_t f_A_mut = count(get<0>(mut_res).begin(), get<0>(mut_res).end(), 'A');
                                    size_t f_C_mut = count(get<0>(mut_res).begin(), get<0>(mut_res).end(), 'C');
                                    size_t f_G_mut = count(get<0>(mut_res).begin(), get<0>(mut_res).end(), 'G');
                                    size_t f_T_mut = count(get<0>(mut_res).begin(), get<0>(mut_res).end(), 'T');

                                    fout << freq_A << ' ' << freq_C << ' ' << freq_G << ' ' << freq_T << ' '
                                    << p_s << ' ' << p_d << ' ' << d << ' ' << ksize << ' ' << L << ' ' <<
                                    L2 << ' ' << L - ksize + 1 << ' ' << L2 - ksize + 1 << ' ' << S << ' ' << D
                                    << ' ' << I << ' ' << f_A_orig << ' ' << f_A_mut << ' ' << f_C_orig
                                    << ' ' << f_C_mut << ' ' << f_G_orig << ' ' << f_G_mut << ' '
                                    << f_T_orig << ' ' << f_T_mut << ' ' << S_sp << ' ' << D_sp << ' ' << I_sp << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}
