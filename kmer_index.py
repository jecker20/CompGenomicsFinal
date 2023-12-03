class KMerIndex:
    def __init__(self, k, sequences):
        self.k = k
        self.data = sequences
        self.index = {}
        self.counts = {}

        for name, seq in sequences.items():
            self.add_sequence(name, seq)

        self.most_shared = list(self.index.keys())
        self.most_shared.sort(key = lambda s: self.counts[s], reverse=True)

    def add_sequence(self, name, seq):
        if self.k == 0 or len(seq) < self.k:
            return
        
        for i in range(len(seq) - self.k + 1):
            s = seq[i:i+self.k]
            if s not in self.index:
                self.index[s] = []
                self.counts[s] = 0
            self.index[s].append((name, i))
            self.counts[s] += 1

    def size(self):
        return len(self.index)

    def get(self, key):
        try:
            return self.index[key]
        except:
            return []

    def contains(self, key):
        return key in self.index
