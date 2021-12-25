

def get_skew_diag_data(genome):
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
    skew_values = get_skew_diag_data(genome)
    min_val = min(skew_values)
    result = []
    for i, item in enumerate(skew_values):
        if item == min_val:
            result.append(i)
    return result

def get_hamming_mismatch(genome1, genome2):
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
    i = 0
    output = []
    while i < len(genome) - len(pattern) + 1:
        if get_hamming_mismatch(pattern, genome[i:i+len(pattern)]) <= mismatch_count:
            output.append(i)
        i += 1
    return output


def get_approximate_pattern_count(genome, pattern, mismatch_count):
    i = 0
    count = 0
    while i < len(genome) - len(pattern) + 1:
        if get_hamming_mismatch(pattern, genome[i:i+len(pattern)]) <= mismatch_count:
            count += 1
        i += 1
    return count

def frequent_words(text, k):
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

from neighbors import neighbors



def frequent_words_with_mismatches(Text, k, d):
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


def frequent_words_with_mismatches_and_rc(Text, k, d):
    patterns = []
    freqMap = {}
    i = 0
    while i <= len(Text)-k:
        pattern = Text[i:i + k]
        forward_neighborhood = list(neighbors(pattern, d))
        reverse_neighbourhood = list(neighbors(reverse_compliment(pattern), d))
        for item in forward_neighborhood:
            freqMap[item] = freqMap.get(item, 0) + 1
        for item in reverse_neighbourhood:
            freqMap[item] = freqMap.get(item, 0) + 1
        i = i+1
    max_count = freqMap[max(freqMap, key=freqMap.get)]
    for key in freqMap:
        if freqMap[key] == max_count:
            patterns.append(key)
    return patterns


def reverse_compliment(text):
    Compliment = {'C': 'G',
                  'G': 'C',
                  'A': 'T',
                  'T': 'A'}
    output = []
    for key in text:
        output.append(Compliment[key])
    reverse_output = (output[::-1])
    output_str1 = ''.join(reverse_output)
    return output_str1



