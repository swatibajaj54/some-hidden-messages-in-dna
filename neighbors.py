from ori_finder import get_hamming_mismatch


def suffix(pattern):
    return pattern[1:]


def neighbors(pattern, d):
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

