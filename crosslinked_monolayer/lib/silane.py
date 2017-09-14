import mbuild as mb


class Silane(mb.Compound):
    """An Si(OH)2 group with two open, opposite ports.

    Ports
    -----
    'down' : orientation=(0,0,-1), shift=0.07, anchor=silicon
    'up' : orientation=(0,0,1), shift=0.07, anchor=silicon

    """
    def __init__(self):
        super(Silane, self).__init__()

        mb.load('silane.pdb', compound=self, relative_to_module=self.__module__)
        mb.x_axis_transform(self, new_origin=self[0], point_on_x_axis=self[1])

        self.add(mb.Port(anchor=self[0], orientation=[0,0,-1], shift=0.07),
                 name='down')
