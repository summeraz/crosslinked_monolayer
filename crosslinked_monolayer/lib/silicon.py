import mbuild as mb
import numpy as np


class Silicon(mb.Compound):
    """A silicon atom with a four ports.

    Ports
    -----
    'up' : Oriented along (0, 0, 1), shifted 0.07nm from the atomic center
    'down' : Oriented along (0, 0, -1), shifted 0.07nm from the atomic center
    'right' : Oriented along (1, 0, 0), shifted 0.07nm from the atomic center
    'left' : Oriented along (-1, 0, 0), shifted 0.07nm from the atomic center

    """
    def __init__(self):
        super(Silicon, self).__init__()
        self.add(mb.Particle(name='Si'))

        self.add(mb.Port(anchor=self[0], orientation=[0, 0, 1], separation=0.07),
                 label='up')
        self.add(mb.Port(anchor=self[0], orientation=[0, 0, -1], separation=0.07),
                 label='down')
        self.add(mb.Port(anchor=self[0], orientation=[1, 0, 0], separation=0.07),
                 label='right')
        self.add(mb.Port(anchor=self[0], orientation=[-1, 0, 0], separation=0.07),
                 label='left')
