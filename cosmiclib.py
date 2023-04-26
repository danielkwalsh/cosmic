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
 
 
# Need to move alpha, beta to the base compressor class   
# Counterfactual stuff goes here
class Counterfactualgenerator:
    def __init__(self, compressors = []):
        self.compressors = compressors
        self.ncomps = len(self.compressors)
        self.alpha = [] #  can be moved to the base class 
        self.beta = [] # can be moved to the base class 
        self.theta = []
        self.q_cf = [None] * self.ncomps
        self.w_cf = []
        self.I = range(self.ncomps)
        self.index_map = {c.name: i for i, c in enumerate(self.compressors)} # Creates an index map (maps compressor name to index)
        self.flag = 1 # indicator of infeasibility (set to 0 when infeasible)

    def find_counterfactual(self, Q, cf_option, sequence_cf = []): # sequence_cf is neeeded when cf_option = 1
        Q_rem = Q
        self.cf_sequence = sequence_cf # Preferred sequence set by the user
        if Q > sum(self.compressors[i].qmax for i, c in enumerate(self.compressors)):
            self.q_cf = [None] * self.ncomps
            self.flag = 0
        else:
            if cf_option == 1: # smart CF (incorporates a preferred sequence)
                for i, c in enumerate(self.cf_sequence):
                    compressor_index = self.index_map[c]
                    qmin_i = self.compressors[compressor_index].qmin
                    qmax_i = self.compressors[compressor_index].qmax
                    if 0 < Q_rem < qmin_i:
                        self.q_cf = [None] * self.ncomps
                        self.flag = 0
                    elif Q_rem > qmax_i:
                        self.q_cf[compressor_index] = qmax_i
                        Q_rem = Q_rem- float(self.q_cf[compressor_index])                    
                    else:
                        self.q_cf[compressor_index] = Q_rem
                        Q_rem = Q_rem- float(self.q_cf[compressor_index])
            if cf_option == 2: # Less smart CF (distributes load evenly among compressors)
                den = sum(self.compressors[i].qmax for i, c in enumerate(self.compressors)) # Sum of max capacitiies of all compressors 
                infeasible_compressors = []
                for i, c in enumerate(self.compressors): # This loop sets initial assignment 
                    self.q_cf[i] = Q*c.qmax/den # Initial assignment of compressor load proportional to its max capacity
                    if self.q_cf[i] < c.qmin: # checks feasibility of assignment 
                        self.q_cf[i] = 0 # sets assigned load to zero for infeasible compressors
                        infeasible_compressors.append(i) # Keeps track of infeasible compressors 
                        den = den - c.qmax # adjusts denominator by removing infeasible compressors from the equation 
                for i, c in enumerate(self.compressors): # This loop does reassignment accounting for infeasibility 
                    if i in infeasible_compressors:
                        self.q_cf[i] = Q*c.qmax/den
            else:
                print('wrong counterfactual option selected')

        #return self.q_cf
    
    def post_process_counterfactual(self):
        self. w_cf =  [None] * self.ncomps
        self.theta_cf = [None] * self.ncomps
        outdict_cf = {'compressors': dict()}
        for i, c in enumerate(self.compressors):
            self.alpha.append((c.wmax - c.wmin)/(c.qmax - c.qmin))
            self.beta.append(c.wmin - c.qmin * self.alpha[i])
            self.theta_cf[i] = None if self.q_cf[i] is None else (self.q_cf[i]>0)
            self.w_cf[i] = None if self.q_cf[i] is None else self.alpha[i] * self.q_cf[i] + self.beta[i] - self.beta[i] * (1 - self.theta_cf[i])
            outdict_cf['compressors'].update({c.name: {'theta': self.theta_cf[i], 'q': self.q_cf[i], 'w': self.w_cf[i]}})
        outdict_cf['Qtot'] = None if self.flag == 0 else sum(self.q_cf)
        outdict_cf['Wtot'] = None if self.flag == 0 else sum(self.w_cf)
        outdict_cf['Success'] = (self.flag ==1)
        outdict_cf['OptimizationStatus'] = None
        return outdict_cf


            
    