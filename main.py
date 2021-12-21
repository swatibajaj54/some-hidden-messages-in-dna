from ori_finder import get_skew_diag_data, get_hamming_mismatch, get_pattern_match, get_approximate_pattern_match

GENOME = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"

if __name__ == '__main__':
    # data = get_skew_diag_data(GENOME)
   # print(get_hamming_mismatch('GGGCCGTTGGT', 'GGACCGTTGAC'))

   # print(get_hamming_mismatch('abscabaxab', 'ab'))
   # print(get_pattern_match('abscabaxab', 'ab'))
   print(get_approximate_pattern_match('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 'ATTCTGGA', 3))