from warnings import warn

import mbuild as mb

from crosslinked_monolayer.lib import CH2, CH3


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port.

    """
    def __init__(self, n, cap_front=True, cap_end=True):
        if n < 2:
            raise ValueError('n must be 2 or more')
        super(Alkane, self).__init__()

        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.Polymer(CH2(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(CH3(), 'methyl_front')
            mb.force_overlap(move_this=self['chain'],
                             from_positions=self['chain']['up'],
                             to_positions=self['methyl_front']['up'])
        else:
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'], self['chain']['down'])
        else:
            self.add(chain['down'], 'down', containment=False)
