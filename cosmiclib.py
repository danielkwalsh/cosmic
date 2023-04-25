from mip import Model, xsum, minimize, BINARY, OptimizationStatus

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

        self.qtot_constr = None


class Optimizer:
    def __init__(self, compressors = []):
        self.compressors = compressors
        self.ncomps = len(self.compressors)
        self.alpha = []
        self.beta = []
        self.theta = []

        self.q = []
        self.w = []

        self.m = Model()
        self.I = range(self.ncomps)

    def replace_Qtot_constr(self, Q):
        constr = self.m.constr_by_name('Qtot')
        if constr is not None:
            self.m.remove(constr)
        self.m += xsum(self.q[i] for i in self.I) == Q, 'Qtot'

    # def compute_W(self):

    def setup_problem(self):
        for i, c in enumerate(self.compressors):
            self.alpha.append((c.wmax - c.wmin)/(c.qmax - c.qmin))
            self.beta.append(c.wmin - c.qmin * self.alpha[i])
            self.theta.append(self.m.add_var(var_type=BINARY))
            self.q.append(self.m.add_var())
            self.w.append(self.m.add_var())
            #Constraints:
            self.m += self.q[i] <= c.qmax * self.theta[i]
            self.m += self.q[i] >= c.qmin * self.theta[i]
            self.m += self.w[i] == self.alpha[i] * self.q[i] + self.beta[i] - self.beta[i] * (1 - self.theta[i])
        self.replace_Qtot_constr(0)

        self.m.objective = minimize(xsum(self.w[i] for i in self.I))

    def find_opt(self, Q):
        self.replace_Qtot_constr(Q)
        status = self.m.optimize()
        q_opt = []
        w_opt = []
        outdict = {'compressors': dict()}
        for i, c in enumerate(self.compressors):
            theta_opt_raw = self.theta[i].x
            theta_opt = None if theta_opt_raw is None else (theta_opt_raw == 1)
            q_opt.append(self.q[i].x)
            w_opt.append(self.w[i].x)
            outdict['compressors'].update({c.name: {'theta': theta_opt, 'q': q_opt[i], 'w': w_opt[i]}})
        outdict['Qtot'] = sum(q_opt)
        outdict['Wtot'] = sum(w_opt)
        outdict['Success'] = (status == OptimizationStatus.OPTIMAL)
        outdict['OptimizationStatus'] = status
        return outdict