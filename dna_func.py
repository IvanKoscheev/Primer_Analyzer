from tables import complement_table, iupac_codes
import itertools
def reverse_complement(dna):
    '''This function reverses the DNA chain and changes the nucleotides to complementary ones.'''
    rev_dna = ''
    for base in dna[::-1]:
        rev_base = complement_table[base]
        rev_dna += rev_base
    return rev_dna
def degenerate_seq_diff(dna):
    '''Make from sequencce which contains degenerate nucleotides 2 sequences: one with highest Tm and second with lowest Tm'''
    low_seq = ''
    high_seq = ''
    for base in dna:
        if base in 'ATGC':
            low_seq += base
            high_seq += base
        elif base in iupac_codes:
            options = sorted(iupac_codes[base],key=lambda x: 1 if x in "GC" else 0)
            low_seq += options[0]
            high_seq += options[-1]
    return low_seq, high_seq
def gc_percent(dna):
    gc_comp = (dna.count('G')+dna.count('C'))/len(dna)
    return gc_comp
def expand_iupac(dna):
    '''Make a list of all possible sequences from one with degenerate nucleotides'''
    options = [iupac_codes.get(base, [base]) for base in dna]
    all_combinations = list(itertools.product(*options))
    return ["".join(variant) for variant in all_combinations]