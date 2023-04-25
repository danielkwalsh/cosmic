class Compressor:
    def __init__(self, name, qmin=None, qmax=None, wmin=None, wmax=None):
        self.name = name
        self.qmin = qmin
        self.qmax = qmax
        self.wmin = wmin
        self.wmax = wmax

        self.qopt = None
        self.is_on = None
        self.wopt = None
