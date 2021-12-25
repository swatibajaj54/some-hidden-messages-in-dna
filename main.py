from ori_finder import get_skew_diag_data, get_hamming_mismatch, get_pattern_match, minimum_skew_value, get_approximate_pattern_match, frequent_words_with_mismatches, get_approximate_pattern_count
from neighbors import neighbors
from ori_finder import frequent_words_with_mismatches_and_rc
if __name__ == '__main__':
    with open('aprox_frequent_rc.txt') as f:
        Text = f.readline()
    k = 5
    d = 2
    print(frequent_words_with_mismatches_and_rc(Text, k, d))

    # print(d1_neighbours("ATC"))
    # print(neighbors("ACG", 1))
    # with open('approx_pattern_count_data.txt') as f:
    #      genome = f.readline()
    # pattern = 'TCGGA'
    # mismatch_count = 3
    # print(get_approximate_pattern_count(genome, pattern, mismatch_count))
    # pattern = 'ACGAGCATA'
    # mismatch_count = 4
    # res = get_approximate_pattern_match(genome, pattern, mismatch_count)
    # for item in res:
    #     print(item, end=" ")

    # with open('hamming_mismatch_genome1_data.txt') as f:
    #     genome1 = f.readline()
    # with open('hamming_mismatch_genome2_data.txt') as f:
    #     genome2 = f.readline()
    # print(get_hamming_mismatch(genome1, genome2))
    # with open('data.txt') as f:
    #     genome = f.readline()
    # print(minimum_skew_value(genome))hamming_mismatch_genome1_data
    # print(get_skew_diag_data(genome))
   # print(get_hamming_mismatch('GGGCCGTTGGT', 'GGACCGTTGAC'))

   # print(get_hamming_mismatch('abscabaxab', 'ab'))
   # print(get_pattern_match('abscabaxab', 'ab'))
   #print(get_approximate_pattern_count('TTTAGAGCCTTCAGAGG', 'GAGG', 2))



