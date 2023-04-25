from mip import Model, xsum, minimize, BINARY

class Compressor:
    def __init__(self, name, qmin=None, qmax=None, wmin=None, wmax=None): #TODO: Add checks on min/max values
        self.name = name
        self.qmin = qmin
        self.qmax = qmax
        self.wmin = wmin
        self.wmax = wmax

        self.qopt = None
        self.is_on = None
        self.wopt = None


class Optimizer:
    def __init__(self, compressors = []):
        self.compressors = compressors
        self.ncomps = len(self.compressors)
        self.alpha = []
        self.beta = []
        self.theta = []

        self.q = []

        self.m = Model()
        self.I = range(self.ncomps)


    def setup_problem(self):
        for i, c in enumerate(self.compressors):
            self.alpha.append((c.wmax - c.wmin)/(c.qmax - c.qmin))
            self.beta.append(c.wmin - c.qmin * self.alpha[i])
            self.theta.append(self.m.add_var(var_type=BINARY))
            self.q.append(self.m.add_var())
            #Constraints:
            self.m += self.q[i] <= c.qmax * self.theta[i]
            self.m += self.q[i] >= c.qmin * self.theta[i]
        self.m += xsum(self.q[i] for i in self.I) == 0.75

        self.m.objective = minimize(xsum(self.alpha[i] * self.q[i] + self.beta[i] - self.beta[i] * (1 - self.theta[i]) for i in self.I))

    def find_opt(self, Q):
        #Q = 1.75
        self.m.optimize()
        return [bool(self.theta[i].x) for i in self.I], [self.q[i].x for i in self.I]