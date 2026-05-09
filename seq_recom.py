'''This block of function made for searching issues in primer and giving warnings'''
from dna_func import reverse_complement, expand_iupac
from tables import iupac_codes
import re
class Primer_Analyze:
    def hairpins_search(self,seq):
        '''Find all hairpins in sequence'''
        hairpins = []
        loops = []
        rev_hairpins = []
        for dna in expand_iupac(seq):
            i = 0
            while i < len(dna) - 7:
                n = 0
                piece = dna[i:i+4]
                rev_piece = reverse_complement(piece)
                if rev_piece in dna[i+7:]:
                    while (i+7+n+1 <= len(dna)) and (reverse_complement(dna[i:i+4+n+1]) in dna[i+7+n+1:]):
                        n += 1
                    final_piece = dna[i:i+4+n]
                    final_rev_piece = reverse_complement(final_piece)
                    pos_rev = dna.find(final_rev_piece, i+7+n)
                    len_loop = pos_rev - (i + 4 + n)
                    hairpins.append(final_piece)
                    loops.append(len_loop)
                    rev_hairpins.append(final_rev_piece)
                    i += (4 + n)
                else:
                    i += 1
            return zip(hairpins, loops, rev_hairpins)
    def primer_length_recomendation(self,dna):
        if len(dna) < 15:
            return f'Primer too small ({len(dna)}), recomended length: 15-45 nucleotides'
        elif len(dna) > 45:
            return f'Primer too long ({len(dna)}), recomended length: 15-45 nucleotides'
    def poly_sites_search(self,dna):
        '''Find polynucleotide sites'''
        poly_sites = []
        if re.search(r'([ATGC])\1{4,}', dna):
            poly_sites.append("polynucleotide")
        if re.search(r'([ATGC]{2})\1{2,}', dna):
            poly_sites.append("polyDInucleotide")
        if re.search(r'([ATGC]{3})\1{2,}', dna):
            poly_sites.append("polyTRInucleotide")
        if poly_sites:
            return f"found {', '.join(set(poly_sites))} site(s)."
    def tail_gc(self, dna):
        tail = dna[-5:]
        gc_count = sum(tail.count(nucl) for nucl in ['G','C','S','V','B','N'])
        if gc_count > 3:
            return 'Very sticky tail!'
        elif gc_count == 0:
            return 'Very weak tail!'