'''
The Python code in this class Motifs was recycled from the group work previously done by Cristiana Martins and Maria Couto in previous lective years.
'''

import math
from classSequence import Sequence

class Motifs:
    def __init__(self, seqs: list) -> None:
        
        '''
        This class allows the creation of probabilistic profiles and other functions.

        Parameters
        ----------
        seqs : list
            Sequences of DNA, RNA or Amino Acids.
        '''
        
        self.seqs = seqs
        self.size = len(seqs[0])
        self.type = self.validate_seq()
        self.pseudo = 0
        self.n = 4
        if self.type == 'DNA':
            self.alphabet = 'ACTG'
        elif self.type == 'RNA':
            self.alphabet = 'ACUG'
        elif self.type == 'Amino Acid':
            self.alphabet = 'FLIMVSPTAY_HQNKDECWRG-'
            self.n = 20

    
    def validate_seq(self, seq: str) -> str:
        
        '''
        This function determines the type of the sequence: DNA, RNA or Amino Acid.
        It gives an error if the sequence is not recognized.

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or Amino Acid.

        Returns
        -------
        str: str
            Type of the sequence ('DNA', 'RNA' or 'Amino Acid')
        '''
        
        test = Sequence(seq)
        return test.check
        
    
    def validate_all_seqs(self) -> str:
        
        '''
        This function validates if all the sequences of a given list are of the same type.
        If any of the sequences is not recognized or if they are the same type, an error is displayed.

        Returns
        -------
        str
            Sequence of DNA, RNA or Amino Acids.

        Raises
        ------
        TypeError
            The sequences are not of the same type.
        '''
        
        validate = Sequence(self.seqs[0])
        val = validate.check
        for i in self.seqs[1:]:           
            valid = self.validate(i)
            if valid != val:
                raise TypeError("The sequences are not of the same type.") 
        return val
            
    
    def create_profile_pwm(self) -> list:
        
        '''
        This funcation calculates the probabilistic PWM (Position Weighted Matrix) of a given list of sequences.
        The default pseudocount is 0.

        Returns
        -------
        matrix : list
            It returns a list with a dictionary of the probabilistic profile.
        '''
        
        total = self.size + self.n*self.pseudo
        matrix = []
        for line in zip(*self.list_seq): matrix.append({k: (line.count(k)+self.pseudo)/total for k in self.alphabet})
        
        
    def create_profile_pssm(self) -> list:
        '''
        This funcation calculates the probabilistic PSSM (Position Specific Scoring Matrix) of a given list of sequences.
        The default pseudocount is 0.
        
        Returns
        -------
        matrix : list
            It returns a list with a dictionary of the probabilistic profile.
        '''
        
        total = self.size + self.n*self.pseudo
        matrix = []
        for line in list(zip(*self.seqs)): matrix.append({k: math.log2(((line.count(k)+self.pseudo)/total)/(1/self.n)) for k in self.alphabet})
        return matrix


    def prob_seq_pwm(self, seq: str) -> float:
        
        '''
        This function calculates and returns the probability of a given sequence by the PWM profile.
        
        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or Amino Acids.

        Returns
        -------
        p : float
            A float number representative of the probability of the given sequence.
        '''
        
        assert len(seq) == len(Motifs.create_profile_pwm(self)) or len(Motifs.create_profile_pssm(self)), 'Sequence size does not match associated profile size'
        if Motifs.create_profile_pwm(self):
            p = 1
            for b, c in zip(seq, Motifs.create_profile_pwm(self)):
                p *= c[b]
    
    
    def prob_seq_pssm(self, seq: str, profile: list) -> float:
        
        '''
        This function calculates and returns the probability of a given sequence by the PSSM profile.
        
        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or Amino Acids.

        Returns
        -------
        p : float
            A float number representative of the probability of the given sequence.
        '''
        
        assert len(seq) == len(Motifs.create_profile_pssm(self)), 'Sequence size does not match associated profile size'
        if Motifs.create_profile_pssm(self):
            p = 0
            for b, c in zip(seq, Motifs.create_profile_pssm(self)):
                p += c[b]
        else:
            raise TypeError("Type of profile invalid.")
        return round(p, 5)
    
    
    def seq_most_probable_pwm(self, seq: str) -> str:
        
        '''
        This function calculates the probability of each subsequence of a given sequence, and returns the subsequence with the highest probability.

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or Amino Acids.

        Returns
        -------
        str
            Subsequence with the highest probability.
        '''
        
        seqs = [seq[I:I+len(Motifs.create_profile_pwm(self))] for I in range(len(seq) - len(Motifs.create_profile_pwm(self)) + 1)]
        probs = [self.prob_seq_pwm(s, Motifs.create_profile_pwm(self)) for s in seqs]
        score_max = max(probs)
        return [seqs[I] for I,p in enumerate(probs) if p == score_max][0]
    
    
    def seq_most_probable_pssm(self, seq: str) -> str:
        
        '''
        This function calculates the probability of each subsequence of a given sequence, and returns the subsequence with the highest probability.

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or Amino Acids.

        Returns
        -------
        str
            Subsequence with the highest probability.
        '''
        
        seqs = [seq[I:I+len(Motifs.create_profile_pssm(self))] for I in range(len(seq) - len(Motifs.create_profile_pssm(self)) + 1)]
        probs = [self.prob_seq_pssm(s, Motifs.create_profile_pssm(self)) for s in seqs]
        score_max = max(probs)
        return [seqs[I] for I,p in enumerate(probs) if p == score_max][0]
    
    
    def consensus(self) -> str:
        
        '''
        This function creates a sequence consensus between two different sequences using the profile created.
        
        Parameters
        ----------
        seq1 : str
            String sequence of DNA or Amino Acid.
        seq2 : str
            String sequence of DNA or Amino Acid.
        
        Returns
        -------
        str
            Sequence consensus.
        '''
        
        profile = Motifs.create_profile_pwm(self)
        cons = ''
        for i in range(self.size):
            d = {n: profile[n][i] for n in self.alphabet}
            m = max(d.values())
            l = [n for n in self.alphabet if d[n] == m]
            cons = [s+[n] for n in l for s in cons]
        
        for s in cons:
            return ''.join(s)