import mbuild as mb
import numpy as np


class Hydrogen(mb.Compound):
    """A hydrogen atom with a single port.

    Ports
    -----
    'up' : Oriented along (0, 0, 1), shifted 0.07nm from the atomic center

    """
    def __init__(self):
        super(Hydrogen, self).__init__()
        self.add(mb.Particle(name='H'))

        self.add(mb.Port(anchor=self[0], orientation=[0, 0, 1], shift=0.07),
                 name='up')
