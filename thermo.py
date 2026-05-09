from dna_func import reverse_complement, gc_percent
import math
from tables import thermodinam_table, termo_loop_table
class Thermo:
    def __init__(self, primer_conc, ion_one_plus_conc, mg_conc, dntp_conc,target_temp):
        self.primer_conc = primer_conc * 10**(-9)
        self.ion_one_plus_conc = ion_one_plus_conc* 10**(-3)
        self.mg_conc = mg_conc*10**(-3)
        self.dntp_conc = dntp_conc*10**(-3)
        self.mg_eff = float(self.mg_conc) - float(self.dntp_conc)
        self.target_temp = target_temp
    def magnesium_delta(self, dna):
        '''Magnesium correction'''
        ln_mag = math.log(self.mg_eff)
        a, b, c, d = 3.92*10**(-5), -9.11*10**(-6), 6.26*10**(-5), 1.42*10**(-5)
        delta_mag_temp = a + b * ln_mag + gc_percent(dna) * (c + d * ln_mag)
        return delta_mag_temp
    def mono_ion_delta(self, dna):
        '''Monovalent ion corrrection'''
        a, b, c = 4.29*10**(-5), 3.95*10**(-5), 9.40*10**(-6)
        ln_salt = math.log(self.ion_one_plus_conc)
        delta_mono_ion_temp = (a * gc_percent(dna) - b) * ln_salt + c * (ln_salt**2)
        return delta_mono_ion_temp
    def term_correction(self,dna,delta_h_one,delta_s_one):
        '''Terminal nucleotides correction'''
        for a in [0,-1]:
            end_nucleotide = dna[a]
            if end_nucleotide in ['A','T']:
                delta_h_one += 2.3
                delta_s_one += 4.1
            elif end_nucleotide in ['G','C']:
                delta_h_one += 0.1
                delta_s_one -= 2.8
        return delta_h_one, delta_s_one
    def delta_h_s(self,dna):
        '''Nearest neighbours method'''
        delta_s = -5.7
        delta_h = 0.2
        for i in range(len(dna) - 1):
            pair = dna[i:i+2]
            h_pair, s_pair = thermodinam_table[pair]
            delta_s += s_pair
            delta_h += h_pair
        delta_h, delta_s = self.term_correction(dna,delta_h,delta_s)
        return delta_h, delta_s
    def melting_point(self, dna):
        '''Main formula (Owczarzy), which consider Mg, K, Na ions and concentrations of primer and dNTP'''
        delta_h, delta_s = self.delta_h_s(dna)
        if dna == reverse_complement(dna):
            delta_s += (-1.4)
        else:
            self.primer_conc = self.primer_conc/4
        stand_inv_melt_temp = (delta_s+1.987*math.log(self.primer_conc)) / (delta_h*1000)
        if self.mg_eff**(1/2) / self.ion_one_plus_conc < 0.22:
            delta_inv_melt_point = self.mono_ion_delta(dna)
        else:
            delta_inv_melt_point = self.magnesium_delta(dna)
        melt_temp = 1/(stand_inv_melt_temp + delta_inv_melt_point) - 273.15
        return melt_temp
    def delta_g_hairpin(self,hairpin,len_loop,temp):
        '''Calculate delta G for hairpins in primer'''
        delta_h, delta_s = self.delta_h_s(hairpin)
        delta_s += 5.7
        delta_h -= 0.2
        if int(len_loop) <= 6:
            add_dh, add_ds = termo_loop_table[len_loop]
            delta_h += add_dh
            delta_s += add_ds
            delta_g = delta_h - (temp + 273.15)*delta_s/1000
        else:
            dh_loop = 5.5 
            ds_loop = -12.7 - 2.44 * math.log(len_loop / 10)
            delta_h += dh_loop
            delta_s += ds_loop
            delta_g = delta_h - (temp + 273.15) * delta_s / 1000
        return delta_g
        