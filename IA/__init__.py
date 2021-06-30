import numpy as np

class Item:
    __slots__ = ['id', 'fitness', 'values', 'phenotype']
    def __init__(self, fitness=None, values=None, phenotype=None, id=0):
        self.fitness = fitness
        self.values = values
        self.phenotype = phenotype
        self.id = id

    def initialize(self, size: int):
        self.values = np.random.randint(2, size=size)
        self.fitness = float("inf")
        self.phenotype = None

    def evaluate(self, funcion):
        self.fitness = funcion.eval(self)

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __str__(self):
        return f'{self.id}\n{self.fitness}\n{self.phenotype.replace(" ","")}\n{binary_to_decimal(self.values)}'
        #return f'{self.fitness} <- {self.values}\n{self.phenotype}'

    def __repr__(self):
        return f'{self.id}\n{self.fitness}\n{self.phenotype.replace("0x","")}\n{binary_to_decimal(self.values)}'

def binary_to_decimal(values, size=3):
    return np.array([int("".join(str(value) for value in array), 2) for array in np.array_split(values, int(len(values)/size))])
