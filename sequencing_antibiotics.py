from ori_finder import reverse_compliment


def protein_translation(rna_sequence):
    genetic_codons = {
        "AAA": "K",
        "AAC": "N",
        "AAG": "K",
        "AAU": "N",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACU": "T",
        "AGA": "R",
        "AGC": "S",
        "AGG": "R",
        "AGU": "S",
        "AUA": "I",
        "AUC": "I",
        "AUG": "M",
        "AUU": "I",
        "CAA": "Q",
        "CAC": "H",
        "CAG": "Q",
        "CAU": "H",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCU": "P",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGU": "R",
        "CUA": "L",
        "CUC": "L",
        "CUG": "L",
        "CUU": "L",
        "GAA": "E",
        "GAC": "D",
        "GAG": "E",
        "GAU": "D",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCU": "A",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGU": "G",
        "GUA": "V",
        "GUC": "V",
        "GUG": "V",
        "GUU": "V",
        "UAA": "",
        "UAC": "Y",
        "UAG": "",
        "UAU": "Y",
        "UCA": "S",
        "UCC": "S",
        "UCG": "S",
        "UCU": "S",
        "UGA": "",
        "UGC": "C",
        "UGG": "W",
        "UGU": "C",
        "UUA": "L",
        "UUC": "F",
        "UUG": "L",
        "UUU": "F",
    }
    i = 0
    protein_seq = []
    while i <= len(rna_sequence) - 3:
        key = rna_sequence[i : i + 3]
        val = genetic_codons.get(key)
        protein_seq.append(val)
        i = i + 3
    return "".join(protein_seq)


genetic_codons = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "",
    "TAC": "Y",
    "TAG": "",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}

genetic_codons_combined = {
    " ": ["TAA", "TAG", "TGA"],
    "C": ["TGC", "TGT"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "G": ["GGA", "GGC", "GGG", "GGT"],
    "F": ["TTC", "TTT"],
    "E": ["GAA", "GAG"],
    "D": ["GAC", "GAT"],
    "K": ["AAA", "AAG"],
    "I": ["ATA", "ATC", "ATT"],
    "H": ["CAC", "CAT"],
    "N": ["AAC", "AAT"],
    "M": ["ATG"],
    "L": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
    "S": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
    "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
    "Q": ["CAA", "CAG"],
    "P": ["CCA", "CCC", "CCG", "CCT"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "T": ["ACA", "ACC", "ACG", "ACT"],
    "Y": ["TAC", "TAT"],
}


def combine_genetic_codons():
    res = {}
    for key in genetic_codons:
        if key == "":
            key = " "
        res[genetic_codons[key]] = res.get(genetic_codons[key], []) + [key]
    print(res)

def compliment_dna(pattern):
    """returns reverse complement for the input pattern"""
    compliment = {"C": "G", "G": "C", "A": "T", "T": "A"}
    output = []
    for key in pattern:
        output.append(compliment[key])
    return output

def peptide_encoding_seq(dna_seq, peptide_seq):
    genetic_codons = {
        "AAA": "K",
        "AAC": "N",
        "AAG": "K",
        "AAT": "N",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AGA": "R",
        "AGC": "S",
        "AGG": "R",
        "AGT": "S",
        "ATA": "I",
        "ATC": "I",
        "ATG": "M",
        "ATT": "I",
        "CAA": "Q",
        "CAC": "H",
        "CAG": "Q",
        "CAT": "H",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "GAA": "E",
        "GAC": "D",
        "GAG": "E",
        "GAT": "D",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "TAA": "*",
        "TAC": "Y",
        "TAG": "*",
        "TAT": "Y",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TGA": "*",
        "TGC": "C",
        "TGG": "W",
        "TGT": "C",
        "TTA": "L",
        "TTC": "F",
        "TTG": "L",
        "TTT": "F",
    }

    temp = []
    i = 0
    while i <= len(dna_seq) - 3:
        key = dna_seq[i : i + 3]
        val = genetic_codons.get(key)
        temp.append(val)
        i = i + 3
    protein_seq = "".join(temp)
    print(protein_seq)
    dna_encoding_seq = []
    i = 0
    while i <= len(protein_seq):
        if peptide_seq == protein_seq[i : i + len(peptide_seq)]:
            dna_encoding_seq.append(dna_seq[i * 3: (i * 3) + len(peptide_seq) * 3])
        i = i + 1
    print('initial dna_encoding_seq', dna_encoding_seq)
    reverse_comp_seq = set()
    for item in dna_encoding_seq:
        reverse_comp_seq.add(reverse_compliment(item))
    s = list(reverse_comp_seq)
    print('reverse_comp_seq', s)
    i = 0
    while i < len(s):
        count = dna_seq.count(s[i])
        if count > 0:
            for _ in range(count):
                dna_encoding_seq.append(s[i])
        i = i+1

    print('final dna_encoding_seq', dna_encoding_seq)

def get_genome_encoding_substrings(dna_seq, peptide_seq):
    genetic_codons = {
        "AAA": "K",
        "AAC": "N",
        "AAG": "K",
        "AAT": "N",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AGA": "R",
        "AGC": "S",
        "AGG": "R",
        "AGT": "S",
        "ATA": "I",
        "ATC": "I",
        "ATG": "M",
        "ATT": "I",
        "CAA": "Q",
        "CAC": "H",
        "CAG": "Q",
        "CAT": "H",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "GAA": "E",
        "GAC": "D",
        "GAG": "E",
        "GAT": "D",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "TAA": "*",
        "TAC": "Y",
        "TAG": "*",
        "TAT": "Y",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TGA": "*",
        "TGC": "C",
        "TGG": "W",
        "TGT": "C",
        "TTA": "L",
        "TTC": "F",
        "TTG": "L",
        "TTT": "F",
    }

    temp = []
    i = 0
    while i <= len(dna_seq) - 3:
        key = dna_seq[i : i + 3]
        val = genetic_codons.get(key)
        temp.append(val)
        i = i + 3
    protein_seq = "".join(temp)
    print(protein_seq)
    dna_encoding_seq = []
    i = 0
    while i <= len(protein_seq):
        if peptide_seq == protein_seq[i : i + len(peptide_seq)]:
            dna_encoding_seq.append(dna_seq[i * 3: (i * 3) + len(peptide_seq) * 3])
        i = i + 1
    return dna_encoding_seq



    # compliment_strand = compliment_dna(dna_seq)
    # print(''.join(compliment_strand))
    # i = 0
    # temp_1 = []
    # while i <= len(reverse_seq) - 3:
    #     key = reverse_seq[i: i + 3]
    #     val = genetic_codons.get(key)
    #     temp_1.append(val)
    #     i = i + 3
    # protein_seq_1 = "".join(temp)
    # print(protein_seq_1)
    # dna_encoding_seq_1 = []
    # i = 0
    # while i <= len(protein_seq_1):
    #     if peptide_seq == protein_seq_1[i: i + len(peptide_seq)]:
    #         dna_encoding_seq_1.append(reverse_seq[i * 3: (i * 3) + len(peptide_seq) * 3])
    #     i = i + 1
    # print(dna_encoding_seq_1)
    # dna_encoding_seq_set = set(dna_encoding_seq)
    # for item in dna_encoding_seq_set:
    #     dna_encoding_seq.append(reverse_compliment(item))
    # print(dna_encoding_seq)
    # reverse_seq = reverse_compliment(dna_seq)
    # print(reverse_seq)
    # total_seq_len = 1
    # for ch in peptide_seq:
    #     seq_len = len(genetic_codons_combined[ch])
    #     total_seq_len *= seq_len
    # print(total_seq_len)
    #
    # all_seq = [''] * total_seq_len
    # for (ch, i) in peptide_seq:
    #     count = len(genetic_codons_combined[ch])
    #     reps = total_seq_len/count
    #     # Repeat count times
    #     current_seq_list = genetic_codons_combined[ch]
    #     for i in range(0, count*reps):
    #         all_seq[i] += current_seq_list[]

    def linear_spectrum(peptide, alphabet, amino_acid_mass):
        prefix_mass = [0]
        i = 0
        while i <= len(peptide):
            for symbol in alphabet:
                if symbol == peptide[i]:
                    prefix_mass[i] = prefix_mass[i-1] + amino_acid_mass[symbol]
            i = i +1
            j = 0
        linear_spectrum = [0]
        while i < len(peptide) - 1:
            while j < len(peptide) + i + 1:
                linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
                j = j+1
            i = i+1
        return linear_spectrum.sort()




# LinearSpectrum(Peptide, Alphabet, AminoAcidMass)
#     PrefixMass(0) ← 0
#     for i ← 1 to |Peptide|
#         for every symbol s in Alphabet
#             if s = i-th amino acid in Peptide
#                 PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass[s]
#     LinearSpectrum ← a list consisting of the single integer 0
#     for i ← 0 to |Peptide| − 1
#         for j ← i + 1 to |Peptide|
#             add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
#     return sorted list LinearSpectrum

