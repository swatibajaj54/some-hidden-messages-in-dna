


def count(text, pattern):
    """How many times the given pattern appears in the text"""
    i = 0
    n = 0
    while i <= len(text) - len(pattern):
        if text[i:i+len(pattern)] == pattern:
            n = n + 1
        i = i + 1
    return n



def reverse_compliment(pattern):
    """returns reverse complement for the input pattern"""
    compliment = {'C': 'G',
                  'G': 'C',
                  'A': 'T',
                  'T': 'A'}
    output = []
    for key in pattern:
        output.append(compliment[key])
    reverse_output = (output[::-1])
    output_str1 = ''.join(reverse_output)
    return output_str1


def frequent_words(text, k):
    """returns most frequent k_mers of length k in the text"""
    i = 0
    d = {}
    while i <= len(text) - k:
        pattern = text[i:i + k]
        d[pattern] = d.get(pattern, 0) + 1
        i = i + 1
    print(d)
    max_count = d[max(d, key=d.get)]
    print(max_count)
    frequent_pattern = []
    for key in d:
        if d[key] == max_count:
            frequent_pattern.append(key)
    return frequent_pattern




def frequency_table(text, k):
    """returns frequency table for k_mers of length k in the text"""
    d = {}
    i = 0
    while i <= len(text)-k:
        pattern = text[i:i+k]
        d[pattern] = d.get(pattern, 0) + 1
        i = i+1
    return d




def find_clumps(text, k, L, t):
    """returns all distinct k-mers forming (L, t)-clumps in text
    ( an interval of Genome of length L in which this k-mer appears at least t times)"""
    i = 0
    my_result = set()
    while i <= len(text) - L:
        window = text[i:i+L]
        freq_map = frequency_table(window, k)
        """FrequencyTable function will produce a frequency table for a given window of a string of length L"""
        for key in freq_map:
            if freq_map[key] >= t:
                my_result.add(key)
        i = i+1
    return list(my_result)




def get_skew_diag_data(genome):
    """returns a skew diagram"""
    output = 0
    result = []
    # print(output, end=' ')
    result.append(output)
    for i in genome:
        if i == 'C':
            output = output-1
            # print(output, end=' ')
            result.append(output)
        elif i == 'G':
            output = output+1
            # print(output, end=' ')
            result.append(output)
        else:
            # print(output, end=' ')
            result.append(output)
    return result


def minimum_skew_value(genome):
    """returns a position in a genome where the skew diagram attains a minimum"""
    skew_values = get_skew_diag_data(genome)
    min_val = min(skew_values)
    result = []
    for i, item in enumerate(skew_values):
        if item == min_val:
            result.append(i)
    return result




def get_hamming_mismatch(genome1, genome2):
    """returns the number of mismatches between strings"""
    count = 0
    for idx, item in enumerate(genome1):
        if genome1[idx] != genome2[idx]:
            count += 1
    return count


def get_pattern_match(genome, pattern):
    i = 0
    output = []
    while i < len(genome) - len(pattern)+1:
        if genome[i:i + len(pattern)] == pattern:
            output.append(i)
        i = i+1
    return output


def get_approximate_pattern_match(genome, pattern, mismatch_count):
    """returns all starting positions where Pattern appears as a substring of Text with at most d mismatches"""
    i = 0
    output = []
    while i < len(genome) - len(pattern) + 1:
        if get_hamming_mismatch(pattern, genome[i:i+len(pattern)]) <= mismatch_count:
            output.append(i)
        i += 1
    return output


def get_approximate_pattern_count(genome, pattern, mismatch_count):
    """ returns the total number of occurrences of Pattern in Text with at most d mismatches"""
    i = 0
    count = 0
    while i < len(genome) - len(pattern) + 1:
        if get_hamming_mismatch(pattern, genome[i:i+len(pattern)]) <= mismatch_count:
            count += 1
        i += 1
    return count




def suffix(pattern):
    return pattern[1:]


def neighbors(pattern, d):
    """returns set of all k-mers whose Hamming distance from Pattern does not exceed d."""
    nucleotide = ["A", "T", "C", "G"]
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'G', 'T', 'C'}
    neighborhood = set()
    suffix_pattern = suffix(pattern)
    suffix_neighbors = neighbors(suffix_pattern, d)
    for text in suffix_neighbors:
        if get_hamming_mismatch(suffix_pattern, text) < d:
            for x in nucleotide:
                neighborhood.add(x+text)
        else:
            neighborhood.add(pattern[0]+text)
    return neighborhood



def frequent_words_with_mismatches(Text, k, d):
    """All most frequent k-mers with up to d mismatches in Text"""
    patterns = []
    freqMap = {}
    i = 0
    while i <= len(Text)-k:
        pattern = Text[i:i + k]
        neighborhood = list(neighbors(pattern, d))
        # print(neighborhood)
        for item in neighborhood:
            freqMap[item] = freqMap.get(item, 0) + 1
        i = i+1
    max_count = freqMap[max(freqMap, key=freqMap.get)]
    print(max_count)
    print(freqMap)
    for key in freqMap:
        if freqMap[key] == max_count:
            patterns.append(key)
    return patterns


def frequent_words_with_mismatches_and_rc(text, k, d):
    """returns all k-mers Pattern maximizing the sum Count_d(Text, Pattern)
    + Count_d(Text, Pattern_rc) over all possible k-mers"""
    patterns = []
    freq_map = {}
    i = 0
    while i <= len(text)-k:
        pattern = text[i:i + k]
        forward_neighborhood = list(neighbors(pattern, d))
        reverse_neighbourhood = list(neighbors(reverse_compliment(pattern), d))
        for item in forward_neighborhood:
            freq_map[item] = freq_map.get(item, 0) + 1
        for item in reverse_neighbourhood:
            freq_map[item] = freq_map.get(item, 0) + 1
        i = i+1
    max_count = freq_map[max(freq_map, key=freq_map.get)]
    for key in freq_map:
        if freq_map[key] == max_count:
            patterns.append(key)
    print(max_count)
    return patterns
































