import math
from math import log2




def exists(pattern, dna_list):
    for dnaString in dna_list:
        if pattern not in dnaString:
            return False
    return True


def variantExists(pattern, dna_list, d):
    variants = list(neighbors(pattern, d))
    for variant in variants:
        if exists(variant, dna_list):
            return True
    return False



def entropy_column(input):
    nucleotides = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for item in input:
        nucleotides[item] = nucleotides[item]+1
    for key in nucleotides:
        temp_res = nucleotides[key]/len(input)
        if temp_res > 0:
            nucleotides[key] = temp_res * abs(log2(temp_res))
        else:
            continue
    sum = 0
    for key in nucleotides:
        sum = sum + nucleotides[key]
    # print(nucleotides)
    return sum

def entropy_matrix(twoD_input):
    sum = 0
    for item in twoD_input:
        print(item)
        sum = sum + entropy_column(item)
    return sum



# def motif_enumeration(Dna, k, d):
#     patterns = set()
#     for item in Dna:
#         i = 0
#         while i <= len(item)-k:
#             pattern = item[i:i+k]
#             neighborhood = list(neighbors(pattern, d))
#             print(pattern, ' ', neighborhood)
#             for approx_pattern in neighborhood:
#                 if variantExists(approx_pattern, Dna, d):
#                     patterns.add(approx_pattern)
#             i = i+1
#     return patterns


def motif_enumeration1(seq, k, d):
    kmer_list = [set() for _ in seq] # Creating a list of sets length len(seq)
    for pos, pattern in enumerate(seq):
        for k_pos in range(len(pattern) - k + 1):
            # Generate neighbors for all kmers in all strings and
            patn = pattern[k_pos:k_pos+k]
            # add them to a set() in a kmer_list
            kmer_list[pos].update(neighbors(patn, d))

    patterns = kmer_list[0]
    # for k_set in kmer_list:
    #     res = set1 & set2 & set3 & set4
    res = kmer_list[0]
    for pos in range(1, len(kmer_list)):
        res = res & kmer_list[pos]

    return res



def distance_between_pattern_and_string(pattern, dna):
    k = len(pattern)
    distance = 0
    for string in dna:
        hamming_distance = math.inf
        i = 0
        while i <= len(string) - k:
            k_mer = string[i:i + k]
            if hamming_distance > get_hamming_mismatch(pattern, k_mer):
                hamming_distance = get_hamming_mismatch(pattern, k_mer)
            i = i + 1
        distance = distance + hamming_distance
    return distance


# MedianString(Dna, k)
#     distance ← ∞
#     for each k-mer Pattern from AA…AA to TT…TT
#         if distance > d(Pattern, Dna)
#              distance ← d(Pattern, Dna)
#              Median ← Pattern
#     return Median


def median_string(dna, k):
    median = ''
    distance = math.inf
    for string in dna:
        i = 0
        k_mer = string[i:i + k]
        while i <= len(string) - k:
            if distance > distance_between_pattern_and_string(k_mer, dna):
                distance = distance_between_pattern_and_string(k_mer, dna)
                print(distance)
                median = k_mer
            i = i + 1
    return median


def profile_most_probable_k_mer(text, k, matrix_profile):
    with open('input.txt') as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    text = lines[0]
    k = int(lines[1])
    probability = {
        'A': [float(i) for i in lines[2].split()],
        'C': [float(i) for i in lines[3].split()],
        'G': [float(i) for i in lines[4].split()],
        'T': [float(i) for i in lines[5].split()]
    }
    i = 0
    output = 0
    most_probable_k_mer = ''
    while i <= len(text)-k:
        k_mer = text[i:i+k]
        current_value = 1
        for index, char in enumerate(k_mer):
            current_value = current_value * float(probability.get(char)[index])
        if current_value > output:
            output = current_value
            most_probable_k_mer = k_mer
        i = i+1
    return most_probable_k_mer
