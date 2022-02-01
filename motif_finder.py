import math
from math import log2
from ori_finder import get_hamming_mismatch, neighbors
import random

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
    """returns column entropy of entropy matrix.
    input is motifs"""
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
    """returns entropy of motif matrix"""
    sum = 0
    for item in twoD_input:
        print(item)
        sum = sum + entropy_column(item)
    return sum


def motif_enumeration1(seq, k, d):
    """A brute force approach for motif finding
    returns (k, d)-motifs in Dna, (Given a collection of strings Dna and an integer d,
    a k-mer (string) is a (k,d)-motif if it appears in every string from Dna with at most d mismatches.
    it is slow for large values of k and d)"""
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
    """returns the sum of distances between Pattern and each string in Dna
    it returns a distance"""
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
    """returns a k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers"""
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
    """returns the profile-most probable k-mers (""returns the profile-most probable k-mers in Text.
    k_mer having maximum probability) in Text."""
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


def probability_matrix(motifs, k):
    """returns probability matrix"""
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
    """returns matrix score formed from motifs"""
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
    """A collection of strings (k_mer from each DNA string),
    BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t)"""
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
            input_matrix = probability_matrix(motifs, k)
            i = 1
            while i < t:
                motif_i = profile_most_probable_k_mer(dna[i], k, input_matrix)
                motifs.append(motif_i)
                input_matrix = probability_matrix(motifs, k)
                i=i+1
            new_matrix = score_matrix(motifs, k)
            if new_matrix < best_motifs:
                best_motifs = new_matrix
                result = motifs
    return result



def probability_matrix_with_pseudocounts(motifs, k):
    """returns probability matrix applying Laplace's rule"""
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
    """returns a collection of strings BestMotifs resulting from
    applying GreedyMotifSearch(Dna, k, t) with pseudocounts"""
    ist_k_mers = []
    result = ''
    new_matrix_score = ''
    for string in dna:
          ist_k_mers.append(string[:k])
    best_motifs_score = score_matrix(ist_k_mers, k)
    print(best_motifs_score)
    first_string = dna[0]
    for index, item in enumerate(first_string):
        if index <= len(first_string) - k:
            k_mer = first_string[index: index+k]
            motifs = [k_mer]
            input_matrix = probability_matrix_with_pseudocounts(motifs, k)
            i = 1
            while i < t:
                motif_i = profile_most_probable_k_mer(dna[i], k, input_matrix)
                motifs.append(motif_i)
                input_matrix = probability_matrix_with_pseudocounts(motifs, k)
                i=i+1
            new_matrix_score = score_matrix(motifs, k)
            if new_matrix_score < best_motifs_score:
                best_motifs_score = new_matrix_score
                result = motifs
    return result, new_matrix_score




def best_score_from_random_motif(dna, k):
    """returns a collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) many times"""
    best_score = math.inf
    best_motifs = []
    for i in range(100):
        score, motifs = randomized_motif_search(dna, k)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    print("score: ", best_score)
    print("best motifs: ", best_motifs)
    return best_motifs


def randomized_motif_search(dna, k):
    """returns a collection of best motifs randomly selected from each string
    of Dna. May generate a rather poor set of motifs, thus run many times to generate the best set"""
    motifs = []
    for string in dna:
        random_index = random.randint(0, len(string)-k)
        k_mer = string[random_index:random_index+k]
        motifs.append(k_mer)
    score_best_motifs = score_matrix(motifs, k)
    score_not_improved_times = 0
    while score_not_improved_times < 20:
        profile = probability_matrix_with_pseudocounts(motifs, k)
        probable_motifs = []
        for string in dna:
            probable_motif = profile_most_probable_k_mer(string, k, profile)
            probable_motifs.append(probable_motif)
        motifs = probable_motifs
        score_motifs = score_matrix(motifs, k)
        if score_motifs < score_best_motifs:
            score_best_motifs = score_motifs
            score_not_improved_times = 0
        else:
            score_not_improved_times = score_not_improved_times + 1
    return score_best_motifs, motifs


def profile_randomly_generated_k_mer(string, k, probability):
    """generate random k_mer based on weight as per their probability"""
    i = 0
    k_mers_list = []
    while i <= len(string) - k:
        k_mer = string[i:i+k]
        k_mers_list.append(k_mer)
        i = i + 1
    k_mers_probability = []
    for k_mer in k_mers_list:
        current_value = 1
        for index, char in enumerate(k_mer):
            current_value = current_value * float(probability.get(char)[index])
        k_mers_probability.append(current_value)
    probability_sum = 0
    probability_weights = []
    for item in k_mers_probability:
        probability_sum = probability_sum + item
    for item in k_mers_probability:
        probability_weights.append(item/probability_sum)
    x = random.choices(k_mers_list, weights=probability_weights, k=1)
    return x[0]




def gibbs_sampler(dna, k):
    """GibbsSampler uses this random number generator to select a Profile-randomly generated k-mer at each step"""
    best_motifs = []
    t = len(dna)
    for string in dna:
        random_index = random.randint(0, len(string)-k)
        k_mer = string[random_index:random_index+k]
        best_motifs.append(k_mer)
    score_best_motifs = score_matrix(best_motifs, k)
    motifs = best_motifs
    score_not_improved_times = 0
    while score_not_improved_times < 200:
        random_dna_index = random.randint(0, t-1)
        del motifs[random_dna_index]
        profile = probability_matrix_with_pseudocounts(motifs, k)
        motif_i = profile_randomly_generated_k_mer(dna[random_dna_index], k, profile)
        motifs.insert(random_dna_index, motif_i)
        score_motifs = score_matrix(motifs, k)
        if score_motifs < score_best_motifs:
            score_best_motifs = score_motifs
            score_not_improved_times = 0
        else:
            score_not_improved_times = score_not_improved_times + 1
    return score_best_motifs, motifs


def many_runs_gibbs_sampler(dna, k):
    """returns the strings BestMotifs resulting from running GibbsSampler(Dna, k) with many random starts
    and consensus motif"""
    best_score = math.inf
    best_motifs = []
    for i in range(100):
        score, motifs = gibbs_sampler(dna, k)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    print("score: ", best_score)
    print("best motifs: ", best_motifs)
    consensus_motif_generated = consensus_motif(best_motifs, k)
    return best_motifs, consensus_motif_generated


def consensus_motif(motifs, k):
    """Returns consensus motif"""
    motif_matrix = {"A": [0]*k, "T": [0]*k, "G": [0]*k, "C": [0]*k}
    for motif in motifs:
        for index, char in enumerate(motif):
            motif_matrix[char][index] = motif_matrix[char][index]+1
    i = 0
    motif = []
    while i < k:
        current_value = 0
        max_valued_key = ''
        for key in motif_matrix:
            if motif_matrix[key][i] > current_value:
                current_value = motif_matrix[key][i]
                max_valued_key = key
        motif.append(max_valued_key)
        i = i + 1
    string_motif = ''.join(motif)
    return string_motif










