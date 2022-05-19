from ori_finder import reverse_compliment

def protein_translation (rna_sequence):
    genetic_codons = {'AAA': 'K',
                      'AAC': 'N',
                      'AAG': 'K',
                      'AAU': 'N',
                      'ACA': 'T',
                      'ACC': 'T',
                      'ACG': 'T',
                      'ACU': 'T',
                      'AGA': 'R',
                      'AGC': 'S',
                      'AGG': 'R',
                      'AGU': 'S',
                      'AUA': 'I',
                      'AUC': 'I',
                      'AUG': 'M',
                      'AUU': 'I',
                      'CAA': 'Q',
                      'CAC': 'H',
                      'CAG': 'Q',
                      'CAU': 'H',
                      'CCA': 'P',
                      'CCC': 'P',
                      'CCG': 'P',
                      'CCU': 'P',
                      'CGA': 'R',
                      'CGC': 'R',
                      'CGG': 'R',
                      'CGU': 'R',
                      'CUA': 'L',
                      'CUC': 'L',
                      'CUG': 'L',
                      'CUU': 'L',
                      'GAA': 'E',
                      'GAC': 'D',
                      'GAG': 'E',
                      'GAU': 'D',
                      'GCA': 'A',
                      'GCC': 'A',
                      'GCG': 'A',
                      'GCU': 'A',
                      'GGA': 'G',
                      'GGC': 'G',
                      'GGG': 'G',
                      'GGU': 'G',
                      'GUA': 'V',
                      'GUC': 'V',
                      'GUG': 'V',
                      'GUU': 'V',
                      'UAA': '',
                      'UAC': 'Y',
                      'UAG': '',
                      'UAU': 'Y',
                      'UCA': 'S',
                      'UCC': 'S',
                      'UCG': 'S',
                      'UCU': 'S',
                      'UGA': '',
                      'UGC': 'C',
                      'UGG': 'W',
                      'UGU': 'C',
                      'UUA': 'L',
                      'UUC': 'F',
                      'UUG': 'L',
                      'UUU': 'F'}
    i = 0
    protein_seq = []
    while i <= len(rna_sequence)-3:
        key = rna_sequence[i:i+3]
        val = genetic_codons.get(key)
        protein_seq.append(val)
        i = i+3
    return ''.join(protein_seq)


def peptide_encoding_seq(dna_seq, peptide_seq):
    genetic_codons = {'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT' : 'N',
    'ACA' : 'T',
    'ACC' : 'T',
    'ACG' : 'T',
    'ACT' : 'T',
    'AGA' : 'R',
    'AGC' : 'S',
    'AGG' : 'R',
    'AGT' : 'S',
    'ATA' : 'I',
    'ATC' : 'I',
    'ATG' : 'M',
    'ATT' : 'I',
    'CAA' : 'Q',
    'CAC' : 'H',
    'CAG' : 'Q',
    'CAT' : 'H',
    'CCA' : 'P',
    'CCC' : 'P',
    'CCG' : 'P',
    'CCT' : 'P',
    'CGA' : 'R',
    'CGC' : 'R',
    'CGG' : 'R',
    'CGT' : 'R',
    'CTA' : 'L',
    'CTC' : 'L',
    'CTG' : 'L',
    'CTT' : 'L',
    'GAA' : 'E',
    'GAC' : 'D',
    'GAG' : 'E',
    'GAT' : 'D',
    'GCA' : 'A',
    'GCC' : 'A',
    'GCG' : 'A',
    'GCT' : 'A',
    'GGA' : 'G',
    'GGC' : 'G',
    'GGG' : 'G',
    'GGT' : 'G',
    'GTA' : 'V',
    'GTC' : 'V',
    'GTG' : 'V',
    'GTT' : 'V',
    'TAA' : '',
    'TAC' : 'Y',
    'TAG' : '',
    'TAT' : 'Y',
    'TCA' : 'S',
    'TCC' : 'S',
    'TCG' : 'S',
    'TCT' : 'S',
    'TGA' : '',
    'TGC' : 'C',
    'TGG' : 'W',
    'TGT' : 'C',
    'TTA' : 'L',
    'TTC' : 'F',
    'TTG' : 'L',
    'TTT' : 'F'}
    reverse_seq = reverse_compliment(dna_seq)
    print(reverse_seq)
    i = 0
    seq = []
    while i <= len(dna_seq)-3:
        key = dna_seq[i:i+3]
        val = genetic_codons.get(key)
        full_sequence = ''
        for item in peptide_seq:
            if item == val:
                seq.append(key)
        i = i+3
    print(seq)





