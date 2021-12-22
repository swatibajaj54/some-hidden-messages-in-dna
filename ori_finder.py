

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


def frequent_words_with_mismatches(Text, k, d):
    patterns = [0]
    freqMap = {}
    i = 0
    while i<len(Text)-k:
        neighborhood = neighbors(pattern, d)

