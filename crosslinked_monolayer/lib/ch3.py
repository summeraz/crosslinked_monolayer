import mbuild as mb


class CH3(mb.Compound):
    """A methyl group with an open port.

    Ports
    -----
    'down' : orientation=(0,0,-1), separation=0.07, anchor=carbon

    """
    def __init__(self):
        super(CH3, self).__init__()

        mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
        self.translate(-self[0].pos)

        self.add(mb.Port(anchor=self[0], orientation=[0, 0, -1], separation=0.07),
                 label='down')
