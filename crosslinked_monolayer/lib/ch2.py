import mbuild as mb


class CH2(mb.Compound):
    """A methylene bridge with two open, opposite ports.

    Ports
    -----
    'down' : orientation=(0,0,-1), separation=0.07, anchor=carbon
    'up' : orientation=(0,0,1), separation=0.07, anchor=carbon

    """
    def __init__(self):
        super(CH2, self).__init__()

        mb.load('ch2.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, 0, -1], separation=0.07),
                 label='down')
        self.add(mb.Port(anchor=self[0], orientation=[0, 0, 1], separation=0.07),
                 label='up')
