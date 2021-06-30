import numpy as np
class OneMax:
    def eval(self, value):
        return len(value)-np.sum(value)