import pandas as pd
import argparse
import re
from dna_func import gc_percent, degenerate_seq_diff
from thermo import Thermo
from seq_recom import Primer_Analyze
from tables import iupac_codes
class Core_Report:
    def __init__(self, primer_conc, ion_one_plus_conc, mg_conc, dntp_conc,target_temp):
        self.primer_conc = primer_conc * 10**(-9)
        self.ion_one_plus_conc = ion_one_plus_conc* 10**(-3)
        self.mg_conc = mg_conc*10**(-3)
        self.dntp_conc = dntp_conc*10**(-3)
        self.mg_eff = float(self.mg_conc) - float(self.dntp_conc)
        self.target_temp = target_temp
        self.thermo = Thermo(primer_conc, ion_one_plus_conc, mg_conc, dntp_conc,target_temp)
        self.recom = Primer_Analyze()
    def open_file(self, file_txt):
        '''Open file and read sequences line by line and make them clean(without unnecessary simbols)'''
        primers=[]
        full_seq=[]
        with open(file_txt,'r') as file:
            for line in file:
                seq = line
                clean_line_1 = re.sub(r'[\(\[].*?[\)\]]', '', line)
                clean_line_2 = re.sub(r'[^a-zA-Z]', '', clean_line_1)
                if clean_line_2:
                    primers.append(clean_line_2.upper())
                if seq:
                    full_seq.append(seq)
        return primers, full_seq
    def analyze_primer(self,dna,full_seq):
        '''Make a table with all needed values '''
        warnings = []
        low_seq, high_seq = degenerate_seq_diff(dna)
        if low_seq == high_seq:
            tm_min = round(self.thermo.melting_point(dna), 3)
            tm_max = tm_min
            gc = round(gc_percent(dna)*100, 2)
            diff_temp = round(abs(tm_min - self.target_temp),3)
        else:
            tm_min = round(self.thermo.melting_point(low_seq), 3)
            tm_max = round(self.thermo.melting_point(high_seq), 3)
            gc = f'{round(gc_percent(low_seq)*100, 2)}-{round(gc_percent(high_seq)*100, 2)}'
            diff_temp = round(abs(tm_min - self.target_temp),3)
        len_rec = self.recom.primer_length_recomendation(dna)
        seq_rec = self.recom.poly_sites_search(dna)
        tail_gc = self.recom.tail_gc(dna)
        warnings = [w for w in (len_rec,seq_rec,tail_gc) if w]
        ghair = []
        hairpins = []
        if sum(dna.count(deg_base) for deg_base in iupac_codes) <= 3:
            if self.recom.hairpins_search(dna):
                for hairpin, loop, rev_hairpin in self.recom.hairpins_search(dna):
                    dg_hairpin = self.thermo.delta_g_hairpin(hairpin,loop,tm_min)
                    ghair.append(round(dg_hairpin,2))
                    hairpins.append(f'{hairpin}...{rev_hairpin}')
        else:
            warnings.append('a lot of degenerative nucleotides, hairpins is not considered')
        if '[' and ']' in full_seq and dna[0] == 'G':
            warnings.append('G on 5\' end of probe')
        report =  {
            "Sequence": full_seq.strip(),
            "Tm min": tm_min,
            "Tm max": tm_max,
            "GC comp": gc,
            "Status": "OK" if not warnings else " | ".join(warnings),
            "Hairpins": "None" if not hairpins else " | ".join(hairpins),
            "Hairpins dG": "-" if not ghair else ghair,
            "Diff": diff_temp
        }
        return report
    def get_report(self, file_txt):
        print("Reading file")
        primers, names = self.open_file(file_txt)
        data = [self.analyze_primer(p,n) for p,n in zip(primers, names)]
        df = pd.DataFrame(data)
        df_sorted = df.sort_values(by="Diff")
        print("The program is complete")
        return df_sorted
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Primer processing tool')
    parser.add_argument('-p',"--primer_conc",type=float, help="Primers concentration in nM")
    parser.add_argument('-mi',"--monovalent",type=float, help="Monovalent ions concentration in mM")
    parser.add_argument('-m',"--magnesium",type=float, help="Magnesium ions concentration in mM")
    parser.add_argument('-d',"--dntp",type=float, help="dNTP concentration in mM")
    parser.add_argument('-t',"--temp",type=float, help="Polymerase work temperature, melting points will sorted by deviation from this value")
    parser.add_argument('-f',"--file_path",help="Path to the file with primers")
    parser.add_argument('-s',"--spreadsheet",help="Path to the spreadsheet")
    args = parser.parse_args()
    my_lab = Core_Report(primer_conc = args.primer_conc,ion_one_plus_conc = args.monovalent,mg_conc=args.magnesium,dntp_conc=args.dntp,target_temp = args.temp)
    total = my_lab.get_report(args.file_path)
    total.to_csv(args.spreadsheet, index=False, sep=';', encoding='utf-8-sig')
        