from ori_finder import get_skew_diag_data, get_hamming_mismatch, get_pattern_match, minimum_skew_value, get_approximate_pattern_match, frequent_words_with_mismatches, get_approximate_pattern_count
from neighbors import neighbors
from motif_finder import motif_matrix, score_matrix, greedy_motif_search, median_string, motif_enumeration1, profile_most_probable_k_mer, entropy_matrix, distance_between_pattern_and_string
from ori_finder import frequent_words_with_mismatches_and_rc
if __name__ == '__main__':
 dna = ['GGCGGGAT', 'TTCAAAGA', 'AATGGGAT', 'AAATCCGG']
 t = 4
 k = 3
 print(greedy_motif_search(dna, k, t))

 # text = 'abcdefghi'
 # k = 3
 # matrix_profile = 'abcd'
 # print(profile_most_probable_k_mer(text, k, matrix_profile))


 # dna = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']
 # k = 3
 # print(median_string(dna, k))

 # dna = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
 # pattern = 'AAA'
 # print(distance_between_pattern_and_string(pattern, dna))


    # input = ["TCATATTTTT",
    #  "CCCTATCCAC",
    #  "GGGGGGGGGG",
    #  "GGGGGGGGGG",
    #  "GTGGGGGGGG",
    #  "GGGGGGGGGT",
    #  "GAAAAAAAAA",
    #  "TCTCCCTTAT",
    #  "TTTTTTTTCA",
    #  "TTTTTTCCTA",
    #  "TATTCCACAC",
    #  "TCCTCCTTCC"]
    # print(entropy_matrix(input))
    # Dna = ["CTCGAACGTGTACTTCCTTCTCCCC", "CCCCCAATCTATTTCGGCCCCGGGC", "CGTCATCGCCCAACTCCTCCCGATT", "ATAGCTGACTGTTTGACCTGGCTCC", "CGCCAACACCGCCTAGCTAAGTGTG", "ACTGTATATCGGCAACACCCACCCC"]
    # Dna = ['AAAAA', 'AAAAA', 'AAAAA']
    # Dna = ['ACGT', 'ACGT', 'ACGT']
    # Dna =  ["TCTGAGCTTGCGTTATTTTTAGACC" ,"GTTTGACGGGAACCCGACGCCTATA" ,"TTTTAGATTTCCTCAGTCCACTATA" ,"CTTACAATTTCGTTATTTATCTAAT" ,"CAGTAGGAATAGCCACTTTGTTGTA" ,"AAATCCATTAAGGAAAGACGACCGT"]
    # k = 5
    # d = 2
    # res = list(motif_enumeration1(Dna, k, d))
    # for item in res:
    #     print(item, end=" ")


    # genome = open('Salmonella.txt').read().replace('\n', '')

    # # print(minimum_skew_value(genome))
    # Text = genome[3764856: 3764856 + 500]
    # print(Text)
    # k = 9
    # d = 1
    # print(frequent_words_with_mismatches_and_rc(Text, k, d))


# pattern = 'GCGCCGTAC'
# d = 3
# result = neighbors(pattern, d)
# for item in result:
#     print(item, end=" ")


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



