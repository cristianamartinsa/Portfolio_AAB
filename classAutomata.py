'''
The Python code in this class Motifs was recycled from the group work previously done by Cristiana Martins and Maria Couto and classes assisted in previous lective years.
'''

class Automata:
    def __init__(self, alphabet: str, pattern: str):
        
        self.numstates = len(pattern) + 1
        self.alphabet = alphabet
        self.transitionTable = {}
        self.pattern = pattern
        self.buildTransitionTable()   
    
    
    def buildTransitionTable(self):
        
        '''
        This function constructs the transition table for the given pattern.
        '''
        
        self.transitionTable = {(q, a): overlap(self.pattern[0:q] + a, self.pattern) for q in range(self.numstates) for a in self.alphabet}


    def printAutomata(self):
        
        '''
        This function pretty prints the automata data.
        '''
        
        print ("States: " , self.numstates)
        print ("Alphabet: " , self.alphabet)
        print ("Transition table:")
        for k in self.transitionTable.keys():
            print (k[0], ",", k[1], " -> ", self.transitionTable[k])


    def nextState(self, current: int, symbol: str) -> int:
        
        '''
        This function obtains the next state for a give current state "current" and character "symbol"
        
        Parameters
        ----------
        current : int
            Current state.
        symbol : str
            Character of the alphabet.
        
        Returns
        -------
        int
            Next state.
        '''
        
        assert symbol in self.alphabet, "Symbol not found"
        return self.transitionTable[(current, symbol)]
        
    
    def applySeq(self, seq: str) -> list:
        
        '''
        This function returns the states along the given sequence
        
        Parameters
        ----------
        seq : str
            Sequence to apply and verify the state alteration
        Returns
        -------
        list
            List of alteration of states
        '''
        q = 0
        res = [q]
        for i in seq:
            res.append(self.transitionTable[(q, i)])
            q = self.transitionTable[(q, i)]
        return res
       
    
    def occurencesPattern(self, text: str) -> tuple:
        
        '''
        This function gives the occurences of the pattern.
        
        Parameters
        ----------
        text : str
            Model to cross the pattern and identify the occurences.
       
        Returns
        -------
        tuple
            Tuple with the occurences of the pattern in the "text".
        '''
        
        q = 0
        res = []
        for i in range(len(text)):
            q = self.nextState(q, text[i])
            if q == self.numstates - 1:
                res.append(i - self.numstates + 2)
        return res


    def overlap(s1: str, s2: str) -> int:
        
        '''
        This function gives the size of the overlap between two sequences. 
        
        Parameters
        ----------
        s1 : str
            Sequence 1.
        s2 : str
            Sequencxe 2.
        
        Returns
        -------
        int
            Overlap size between the two sequences.
        '''
        
        maxov = min(len(s1), len(s2))
        for i in range(maxov,0,-1):
            if s1[-i:] == s2[:i]: return i
        return 0