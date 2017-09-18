import os
from pkg_resources import resource_filename

def get_forcefield(name):
    forcefield = resource_filename('crosslinked_monolayer',
        os.path.join('forcefield', '{}.xml'.format(name)))
    if not os.path.exists(forcefield):
        raise IOError("Forcefield '{}' does not exist!".format(name))
    return forcefield
