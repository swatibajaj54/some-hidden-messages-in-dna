import math
from math import log2
from ori_finder import get_hamming_mismatch, neighbors


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


def median_string(dna, k):
    median = ''
    distance = math.inf
    for string in dna:
        i = 0
        temp_list = []
        while i <= len(string) - k:
            k_mer = string[i:i + k]
            neighbours = neighbors(k_mer, 1)
            for neighbouring_k_mer in neighbours:
                d = distance_between_pattern_and_string(neighbouring_k_mer, dna)
                temp_list.append(d)
                if distance > d:
                    distance = d
                    print(distance)
                    median = k_mer
            i = i + 1
    return median


def profile_most_probable_k_mer(text, k, probability):
    i = 0
    output = -1
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


def motif_matrix(motifs, k):
    nucleotides = {'A': [0]*k, 'T': [0]*k, 'C': [0]*k, 'G': [0]*k}
    t = len(motifs)
    for motif in motifs:
        for index, nucleotide in enumerate(motif):
            nucleotides[nucleotide][index] = nucleotides[nucleotide][index] + 1
    # print(nucleotides)
    for key in nucleotides:
        for index, nucleotide in enumerate(nucleotides[key]):
            nucleotides[key][index] = nucleotides[key][index]/t
    return nucleotides


def score_matrix(motifs, k):
    nucleotides = {'A': [0]*k, 'T': [0]*k, 'C': [0]*k, 'G': [0]*k}
    for motif in motifs:
        for index, nucleotide in enumerate(motif):
            nucleotides[nucleotide][index] = nucleotides[nucleotide][index] + 1
    i = 0
    matrix_score = 0
    while i < k:
        output = []
        column_score = 0
        for key in nucleotides:
            output.append(nucleotides[key][i])
        max_consumed = False
        max_item = max(output)
        for item in output:
            if item == max_item:
                if not max_consumed:
                    max_consumed = True
                    continue
                else:
                    column_score = column_score + item
            else:
                column_score = column_score+item
        matrix_score = matrix_score + column_score
        i = i + 1
    return matrix_score


def greedy_motif_search(dna, k, t):
    ist_k_mers = []
    result = ''
    for string in dna:
          ist_k_mers.append(string[:k])
    best_motifs = score_matrix(ist_k_mers, k)
    print(best_motifs)
    first_string = dna[0]
    for index, item in enumerate(first_string):
        if index <= len(first_string) - k:
            k_mer = first_string[index: index+k]
            motifs = [k_mer]
            input_matrix = motif_matrix(motifs, k)
            i = 1
            while i < t:
                motif_i = profile_most_probable_k_mer(dna[i], k, input_matrix)
                motifs.append(motif_i)
                input_matrix = motif_matrix(motifs, k)
                i=i+1
            new_matrix = score_matrix(motifs, k)
            if new_matrix < best_motifs:
                best_motifs = new_matrix
                result = motifs
    return result



def motif_matrix_with_pseudocounts(motifs, k):
    nucleotides = {'A': [1]*k, 'T': [1]*k, 'C': [1]*k, 'G': [1]*k}
    t = len(motifs)
    for motif in motifs:
        for index, nucleotide in enumerate(motif):
            nucleotides[nucleotide][index] = nucleotides[nucleotide][index] + 1
    # print(nucleotides)
    for key in nucleotides:
        for index, nucleotide in enumerate(nucleotides[key]):
            nucleotides[key][index] = nucleotides[key][index]/(t+4)
    return nucleotides



def greedy_motif_search_with_pseudocounts(dna, k, t):
    ist_k_mers = []
    result = ''
    for string in dna:
          ist_k_mers.append(string[:k])
    best_motifs_score = score_matrix(ist_k_mers, k)
    print(best_motifs_score)
    first_string = dna[0]
    for index, item in enumerate(first_string):
        if index <= len(first_string) - k:
            k_mer = first_string[index: index+k]
            motifs = [k_mer]
            input_matrix = motif_matrix_with_pseudocounts(motifs, k)
            i = 1
            while i < t:
                motif_i = profile_most_probable_k_mer(dna[i], k, input_matrix)
                motifs.append(motif_i)
                input_matrix = motif_matrix_with_pseudocounts(motifs, k)
                i=i+1
            new_matrix_score = score_matrix(motifs, k)
            if new_matrix_score < best_motifs_score:
                best_motifs_score = new_matrix_score
                result = motifs
    return result


# RandomizedMotifSearch(Dna, k, t)
#     randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
#     BestMotifs ← Motifs
#     while forever
#         Profile ← Profile(Motifs)
#         Motifs ← Motifs(Profile, Dna)
#         if Score(Motifs) < Score(BestMotifs)
#             BestMotifs ← Motifs
#         else
#             return BestMotifs


import random


def randomized_motif_search(dna, k):
    motifs = []
    for string in dna:
        random_index = random.randint(0, len(string)-k-1)
        k_mer = string[random_index:random_index+k]
        motifs.append(k_mer)
    score_best_motifs = score_matrix(motifs, k)
    same_score_count = 1000
    while True:
        profile = motif_matrix_with_pseudocounts(motifs, k)
        probable_motifs = []
        for string in dna:
            probable_motif = profile_most_probable_k_mer(string, k, profile)
            probable_motifs.append(probable_motif)
        motifs = probable_motifs
        score_motifs = score_matrix(motifs, k)
        if score_motifs < score_best_motifs:
            score_best_motifs = score_motifs
        else:
            same_score_count = same_score_count - 1
            print("Bad score: ", score_motifs)
        if same_score_count == 0:
            print("Score: ", score_best_motifs)
            return motifs


