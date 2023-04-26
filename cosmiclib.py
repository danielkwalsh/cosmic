from mip import Model, xsum, minimize, BINARY, OptimizationStatus
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

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
        self.m.verbose = 0

        self.I = range(self.ncomps)

    def replace_Qtot_constr(self, Q):
        constr = self.m.constr_by_name('Qtot')
        if constr is not None:
            self.m.remove(constr)
        self.m += xsum(self.q[i] for i in self.I) == Q, 'Qtot'


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
        if None in q_opt:
            outdict['Qtot'] = None
            outdict['Wtot'] = None
        else:
            outdict['Qtot'] = sum(q_opt)
            outdict['Wtot'] = sum(w_opt)
        outdict['Success'] = (status == OptimizationStatus.OPTIMAL)
        outdict['OptimizationStatus'] = status
        return outdict

def generate_sequence_plot(opt_obj):
    def make_row(Qtot, sol_dict):
        dict_data = {'Qtot': [Qtot]}
        dict_data.update({f'{c}: q': [sol_dict['compressors'][c]['q']] for c in comp_names})
        dict_data.update({f'{c}: w': [sol_dict['compressors'][c]['w']] for c in comp_names})
        return pd.DataFrame(dict_data)
    compressors = opt_obj.compressors
    Qmax = sum([c.qmax for c in compressors])
    comp_names = [c.name for c in compressors]
    my_optimizer = Optimizer(compressors=compressors)
    my_optimizer.setup_problem()

    df = None
    Qmax_plot = Qmax * 1.1
    for Qtot in np.arange(0,Qmax_plot,0.01):
        solution = my_optimizer.find_opt(Qtot)
        if df is None:
            df = make_row(Qtot, solution)
        else:
            df = pd.concat([df, make_row(Qtot, solution)])

    def get_color(i):
        colors = px.colors.qualitative.Plotly
        size = len(colors)
        meta = i // size
        sub = i % size
        col = px.colors.hex_to_rgb(colors[sub])
        fac = 0.9**meta
        return f'rgb({int(fac*col[0])}, {int(fac*col[1])}, {int(fac*col[2])})'
    # Create a stacked area chart
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    for i, c in enumerate(comp_names):
        fig.add_trace(go.Scatter(x=df['Qtot'], y=df[f'{c}: q'], stackgroup='1', name=c, fillcolor=get_color(i), mode= 'none'), row=1, col=1)
        fig.add_trace(go.Scatter(x=df['Qtot'], y=df[f'{c}: w'], stackgroup='1', name=c, fillcolor=get_color(i), showlegend=False, mode= 'none'), row=2, col=1)
    # Customize the chart layout
    fig.update_layout(title='Optimal Compressor Sequencing', xaxis2_title='Qtot', yaxis_title='q', yaxis2_title='w', height=600)
    # Display the chart
    #fig.show()
    return fig

