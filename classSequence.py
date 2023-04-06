'''
The Python code in this class Sequence was recycled from the group work previously done by Cristiana Martins and Maria Couto in previous lective years.
'''

class Sequence:
    
    def __init__(self, seq: str):
        
        '''
         This function verifies and capitalizes the sequence inputted by the user. 
         
         Parameters
         ----------
         seq : str
             Sequence inputted by the user.
        
         Raises
         ------
         TypeError
             Tells if the sequence is invalid.
        '''
        
        self.seq = seq.upper()
        dna = "ACGT-"
        rna = "ACGU-"
        amino = "FLIMVSPTAY_HQNKDECWRG-"
        
        for i in self.seq:
            if i in dna:
                self.check = "DNA"
            elif i in rna:
                self.check = "RNA"
            elif i in amino:
                self.check = "Amino Acid"
            else:
                raise TypeError("Invalid sequence.")
            
    
    def __str__(self):
    
        '''
        This function gives a string with the sequence and its type.
         
        Returns
        ------------
        str
            String with the information of the sequence.
        '''
    
        return f"This is a {self.check} sequence: {self.seq}"
    
    
    def vog_cons(self) -> str:
        
        '''
        This function counts the sequence vowels and consonants.
        
        Returns
        ------------
        str
            String with the number of vowels and consonantes of the sequence
        '''
        
        v = "AEIOU"
        c = "BCDFGHJKLMNPQRSTVWYZ"
        vog = 0
        con = 0

        for i in self.seq:
            if i in v: 
                vog += 1
            elif i in c: 
                con += 1

        return f'Vogals: {vog}, Consonants: {con}'


    def number(self) -> str:
        
        '''
        This function gives the number of the sequence bases (in case it's a DNA or RNA sequence) or amino acids (in case it's an Amino Acid sequence) and errors (if existed).
        
        Returns
        ------------
        str
            String with the characters and respective number of occurences.
        '''
        
        err = 0
        if self.check == 'DNA':
            dna = ['A', 'C', 'T', 'G']
            a = self.seq.count('A')
            c = self.seq.count('C')
            t = self.seq.count('T')
            g = self.seq.count('G')
            for i in self.seq:
                if i not in dna:
                    err += 1
            return f'A: {a}, C: {c}, T: {t}, G: {g}, Errors: {err}'
        
        elif self.check == 'RNA':
            rna = ['A', 'C', 'U', 'G']
            a = self.seq.count('A')
            c = self.seq.count('C')
            t = self.seq.count('T')
            u = self.seq.count('U')
            for i in self.seq:
                if i not in rna:
                    err += 1
            return f'A: {a}, C: {c}, T: {t}, U: {u}, Errors: {err}'
        
        elif self.check == 'Amino Acid':
            amino_acid = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G']
            f = self.seq.count('F')
            l = self.seq.count('L')
            i = self.seq.count('I')
            m = self.seq.count('M')
            v = self.seq.count('V')
            s = self.seq.count('S')
            p = self.seq.count('P')
            t = self.seq.count('T')
            a = self.seq.count('A')
            y = self.seq.count('Y')
            h = self.seq.count('H')
            q = self.seq.count('Q')
            n = self.seq.count('N')
            k = self.seq.count('K')
            d = self.seq.count('D')
            e = self.seq.count('E')
            c = self.seq.count('C')
            w = self.seq.count('W')
            r = self.seq.count('R')
            g = self.seq.count('G')
            for i in self.seq:
                if i not in amino_acid:
                    err += 1
            return f'F: {f}, L: {l}, I: {i}, M: {m}, V: {v}, S: {s}, P: {p}, T: {t}, A: {a}, Y: {y}, H: {h}, Q: {q}, N: {n}, K: {k}, D: {d}, E: {e}, C: {c}, W: {w}, R: {r}, G: {g}, Errors: {err}'
  
      
    def percentage(self) -> str:
        
        '''
        This function gives the elements of the sequence and corresponding percentage.
        
        Returns
        ------------
        str
            String with the characters and respective percentage.
        '''
        
        assert self.check == ('DNA' or 'RNA'), 'The sequence must be DNA or RNA.'
        leng = len(self.seq)
        
        if self.check == 'DNA':
            a = self.seq.count('A')*100 / leng
            c = self.seq.count('C')*100 / leng
            t = self.seq.count('T')*100 / leng
            g = self.seq.count('G')*100 / leng
            
            return f'A: {a}%, C: {c}%, T: {t}%, G: {g}%'
        
        if self.check == 'RNA':
            a = self.seq.count('A')*100 / leng
            c = self.seq.count('C')*100 / leng
            u = self.seq.count('U')*100 / leng
            g = self.seq.count('G')*100 / leng
            
            return f'A: {a}%, C: {c}%, U: {u}%, G: {g}%'
    
    
    def transcription(self) -> str:
        
        '''
        This function does the following replacements: 
            thymine base 'A' ---> uracil base 'U'
            thymine base 'T' ---> uracil base 'A'
            thymine base 'C' ---> uracil base 'G'
            thymine base 'G' ---> uracil base 'C'.
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        str
            DNA sequence.
        
        Returns
        ------------
        str
            RNA sequence.
        '''
        
        assert self.check == "DNA", "The sequence must be DNA."
        if self.seq: return self.seq.replace("A", "U").replace('C', 'G').replace('G','C').replace('T', 'A')
        else:
             raise Exception('Please, confirm the sequence inputted.')
    
    
    def comp_inverse(self) -> str:
        
        '''
        This function inverses and complements the sequence.
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        seq: str
            DNA sequence.
        
        Returns
        ------------
        str
            String corresponding to the inverse and complement transcription of the sequence inputted.
        '''
    
        comp_str = ''
        for i in self.seq:
            if i == 'A':
                comp_str += 'T'
            elif i == 'T':
                comp_str += 'A'
            elif i == 'C':
                comp_str += 'G'
            elif i == 'G':
                comp_str += 'C'
            else:
                raise Exception('Please, confirm the sequence inputted.') 
        return comp_str[::-1]
    
    
    def get_orfs(self) -> list:
        
        '''
        This function constructs a list with all the 6 orfs.
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        seq: str
            DNA sequence.
               
        Returns
        ------------
        list
            List with each orf.
        '''
        
        assert self.check == "DNA", "The sequence must be DNA."
        seq_rev = Sequence.comp_inverse(self)
    
        orf = []      
        
        orf.append([self.seq[i:i+3] for i in range(0, len(self.seq), 3) if i+3<=len(self.seq)])
        orf.append([self.seq[i+1:i+4] for i in range(0, len(self.seq), 3) if i+4<=len(self.seq)])
        orf.append([self.seq[i+2:i+5] for i in range(0, len(self.seq), 3) if i+5<=len(self.seq)])
        orf.append([seq_rev[i:i+3] for i in range(0, len(seq_rev), 3) if i+3<=len(seq_rev)])
        orf.append([seq_rev[i+1:i+4] for i in range(0, len(seq_rev), 3) if i+4<=len(seq_rev)])
        orf.append([seq_rev[i+2:i+5] for i in range(0, len(seq_rev), 3) if i+5<=len(seq_rev)])
        
        return orf
    
    
    def get_codons(self) -> list:
        
        '''
        This function gives all the codons for the first ORF.
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        seq: str
            DNA sequence.
        
        Returns
        ------------
        list
            List with all the codons.
        '''
        
        assert self.check == "DNA", "The sequence must be DNA."
        cod=[self.seq[i:i+3] for i in range(0, len(self.seq), 3) if i+3<=len(self.seq)]
        return cod 
    
    
    def translation(self) -> str:
        
        '''
        This function gives the sequence codons to build an amino acid chain.
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        seq: str
            DNA sequence.
        
        Returns
        ------------
        str
            Amino acid chain.
        '''
        
        assert self.check == "DNA", "The sequence must be DNA."
        codons = Sequence.get_codons(self) 
             
        amino=""
        if self.seq:
            for i in codons:
                if i == "ATG":
                    amino += "M"
                elif i == "TTT" or i == "TTC":
                    amino += "F"
                elif i[:2] == "TT" or i[:2] == "CT":
                    amino += "L"
                elif i[:2] == "AT":
                    amino += "I"
                elif i[:2] == "GT":
                    amino += "V"
                elif i[:2] == "TC" or i == "AGT" or i == "AGC":
                    amino += "S"
                elif i[:2] == "CC":
                    amino += "P"
                elif i[:2] == "AC":
                    amino += "T"
                elif i[:2] == "GC":
                    amino += "A"
                elif i == "TAA" or i == "TAG" or i == "TGA":
                    amino += "_"
                elif i[:2] == "TA":
                    amino += "Y"
                elif i == "CAT" or i == "CAC":
                    amino += "H"
                elif i[:2] == "CA":
                    amino += "Q"
                elif i == "AAT" or i == "AAC":
                    amino += "N"
                elif i[:2] == "AA":
                    amino += "K"
                elif i == "GAT" or i == "GAC":
                    amino += "D"
                elif i[:2] == "GA":
                    amino += "E"
                elif i == "TGG":
                    amino += "W"
                elif i[:2] == "TG":
                    amino += "C"
                elif i[:2] == "CG" or i[:2] == "AG":
                    amino += "R"
                elif i[:2] == "GG":
                    amino += "G"
            return amino
        else:
             raise Exception('Please, confirm the sequence inputted.')
    
    
    def get_all_prots(self) -> list:
        
        '''
        This function constructs a list with all the possible proteins. 
        Only executes if it is a DNA sequence.
        
        Parameters
        ----------
        seq : str
            DNA or amino acid sequence.

        Returns
        ------------
        list_prots: list
            List with all the proteins.
        '''
        
        assert self.check == ("DNA" or "Amino Acid"), "The sequence must be DNA or Amino Acid."
        
        list_prots = []
        
        if self.check == "DNA":
            amino = Sequence.translation(self)
            proteina = False
            prot = ''
            for i in amino:
                if proteina == False:
                    if i == "M":
                        prot += i
                        proteina = True
                elif proteina == True:
                    if i == "_":
                        prot += i
                        proteina = False
                        list_prots.append(prot)
                        prot = []
                    else:
                        prot += i
        
        elif self.check == "Amino Acid":
            proteina = False
            prot = []
            for i in self.seq:
                if proteina == False:
                    if i == "M":
                        prot += i
                        proteina = True
                elif proteina == True:
                    if i == "_":
                        prot += i
                        proteina = False
                        list_prots.append(prot)
                        prot = []
                    else:
                        prot += i
                
        return list_prots 

    
    def get_bigger_prot(self) -> str:
        
        '''
        This functions gives the bigger protein in the sequence.
        
        Returns
        ------------
        str
            String with the biggest protein.
        '''
        
        all_prots = Sequence.get_all_prots(self)
        bigger_prot = ""
        for i in all_prots:
            if len(i) > len(bigger_prot):
                bigger_prot = i
        return f'Biggest protein: {bigger_prot}, Length: {len(bigger_prot)}'
    
    
    
    