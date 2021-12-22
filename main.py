from ori_finder import get_skew_diag_data, get_hamming_mismatch, get_pattern_match, minimum_skew_value, get_approximate_pattern_match, get_approximate_pattern_count

GENOME = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"

if __name__ == '__main__':
    with open('data.txt') as f:
        genome = f.readline()
    print(minimum_skew_value(genome))
    # print(get_skew_diag_data(genome))
   # print(get_hamming_mismatch('GGGCCGTTGGT', 'GGACCGTTGAC'))

   # print(get_hamming_mismatch('abscabaxab', 'ab'))
   # print(get_pattern_match('abscabaxab', 'ab'))
   #print(get_approximate_pattern_count('TTTAGAGCCTTCAGAGG', 'GAGG', 2))



