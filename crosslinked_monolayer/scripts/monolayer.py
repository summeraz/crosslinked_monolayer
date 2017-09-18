from __future__ import division

from copy import deepcopy
import random
from warnings import warn

import mbuild as mb
import networkx as nx
import numpy as np

from crosslinked_monolayer.lib import Hydroxyl, Silicon


class CrosslinkedMonolayer(mb.Compound):
    """A recipe for constructing an alkylsilane monolayer with crosslinks.

    Crosslinked monolayers created using this recipe will feature two types
    of chains: surface-bound chains and chains attached only via crosslinks.
    A crosslinked monolayer is constructed through the following three stages:
        1. Insert surface-bound chains. This procedure uses a distance-criterion
           to prevent steric overlaps. Each chain is considered to be a cylinder
           with a user specified diameter (`spacing`). Attachment site locations
           are randomly chosen and chains are attached to the surface until no
           available attachment sites remain (i.e. placing chains at any of the
           remaining sites would lead to steric overlap).
        2. Insert crosslinked chains. Crosslinked chains are added to the monolayer
           by choosing random locations in the surface plane for insertion attempts.
           Steric overlap (again, dictating by `spacing`) is checked for each
           insertion attempt, and if no overlap will occur, then a chain is placed
           at this location. This process continues until a user-specified number
           of consecutive failed insertion attempts has occurred
           (`max_failed_attempts`).
        3. This stage takes the chains added to the system in Stage 2 and connects
           them to neighboring chains via crosslinks.

    By default this class will construct the monolayer using a purely steric-based
    approach, dictated by the provided `spacing`. However, if `n_chemisorbed` is
    specified, then (provided this does not lead to steric overlaps) then the
    number of chains attached directly to the surface can be controlled. Note that
    `n_chemisorbed` must be smaller than the maximum number of chemisorbed chains
    that can be placed without overlaps. If a larger value of `n_chemisorbed` is
    provided, then the maximum number of chemisorbed chains without overlaps will
    be included and the user will be provided with a warning message.

    Parameters
    ----------
    chain : mb.Compound
        Prototype of chain to attach to the surface. Should feature an available
        `Port`.
        Note: Silane groups will be added to chain prototypes during monolayer
              creation, so chains should not already feature a silane group.
    surface : mb.Compound
        Surface on which the monolayer will be built. Available binding sites
        should feature `Ports`.
    spacing : float
        Minimum spacing (in nm) allowed between chains. This should usually be the
        Van der Waals diameter of the chain (e.g. 0.42nm for alkylsilanes)
    backfill : mb.Compound, optional, default=None
        Compound to place at vacant surface sites.
    n_chemisorbed : int, optional, default=None
        The exact number of chains that should be chemisorbed to the surface.
        If `None`, which is the default, then the number of chemisorbed chains is
        determined from `spacing`.
    seed : int, optional, default=12345
        Seed for the random number generator 
    chain_port_name : string, optional, default='down'
        Name of the port on `chain` to be used to create bonds to the surface
    backfill_port_name : string, optional, default='down'
        Name of the port on `backfill` to be used to create bonds to the surface
    max_failed_attempts : int, optional, default=2.5e3
        The number of consecutive failed insertion attempts for crosslinked chains
        before the monolayer is considered "complete" (i.e. that no additional
        locations exist where chains could be placed).
    verbose : bool, optional, default=False
        Create the monolayer in verbose mode, where additional information is
        outputted to the screen to inform the user of the status of monolayer
        construction.

    Attributes
    ----------
    crosslink_graph : nx.Graph
        Graph-like object that stores information on the crosslinking network

    To-do
    -----
    Features
    --------
    - More robust routine to ensure chains are oriented in the desired direction
    - Check to make sure chain prototype passed as the `chain` argument does not
      already contain a silane (could be as easy as seeing if the chain contains
      any Si atoms, and just providing a warning).
    - Provide updated chain density in verbose mode

    """
    def __init__(self, chain, surface, spacing, backfill=None, n_chemisorbed=None,
                 seed=12345, chain_port_name='down', backfill_port_name='down',
                 max_failed_attempts=2.5e3, verbose=False):
        super(CrosslinkedMonolayer, self).__init__()

        random.seed(seed)
        np.random.seed(seed)

        self.crosslink_graph = nx.Graph()

        self.add(surface, 'surface')

        # Add a silicon atom to the bottom of the chain prototype
        silicon = Silicon()
        mb.force_overlap(silicon, silicon['up'], chain[chain_port_name])
        chain.add(silicon, 'silicon')
        #mb.force_overlap(chain, chain[chain_port_name], silicon['up'])
        chain.add(chain['silicon']['down'], chain_port_name, containment=False,
            replace=True)

        self._add_chemisorbed_chains(chain, spacing, n_chemisorbed, chain_port_name,
            verbose)
        self._add_crosslinked_chains(chain, spacing, chain_port_name,
            max_failed_attempts, verbose)
        self._determine_crosslink_network(verbose)
        self._create_crosslinks()
        self._add_backfill(backfill, backfill_port_name)

    def _add_chemisorbed_chains(self, chain, spacing, n_chemisorbed,
                                chain_port_name, verbose):
        """Attach chains to the surface...
        """
        available_sites = np.array(self['surface'].available_ports())
        if not n_chemisorbed:
            max_chains = len(available_sites)
        else:
            max_chains = n_chemisorbed
        added_chains = 0
        while len(available_sites) > 0 and added_chains < max_chains:
            binding_site = random.choice(available_sites)
            new_chain = mb.clone(chain)
            self.crosslink_graph.add_node(new_chain.available_ports()[0].anchor,
                pos=(binding_site.pos[0], binding_site.pos[1]),
                surface_bound=True)
            si = list(new_chain.particles_by_name('Si'))[0]
            self.add(new_chain, 'chemisorbed_chain[$]')
            mb.force_overlap(new_chain, new_chain[chain_port_name], binding_site)
            available_sites = [site for site in available_sites
                               if (site != binding_site and
                               self.min_periodic_distance(np.array([site.pos[0],
                               site.pos[1],0]), np.array([binding_site.pos[0],
                               binding_site.pos[1],0])) > spacing)]
            added_chains += 1
            if verbose:
                print('Added chemisorbed chain {}'
                      ''.format(len(self['chemisorbed_chain'])))
                print('{} available sites remaining'.format(len(available_sites)))

        if n_chemisorbed and added_chains < n_chemisorbed:
            warn('Specified number of chemisorbed chains would cause steric '
                 'overlap! Only adding {} chemisorbed chains.'.format(added_chains))
        if len(available_sites) > 0:
            warn('Only adding {} chemisorbed chains; however, additional sites are '
                 'available!'.format(n_chemisorbed))

    def _add_crosslinked_chains(self, chain, spacing, chain_port_name,
                                max_failed_attempts, verbose):
        """Add crosslinked chains... """
        # Add hydroxyl to prototype for crosslinked chain
        chain.spin(np.pi/2, [1,0,0])
        hydroxyl = Hydroxyl()
        mb.force_overlap(hydroxyl, hydroxyl['down'], chain[chain_port_name])
        chain.add(hydroxyl)

        surface_level = max(self['surface'].xyz[:,2])

        failed_attempts = 0
        # Keep adding chains until too many consecutive failed attempts
        while failed_attempts < max_failed_attempts:
            # Choose random site in the xy plane
            site = np.append(np.random.rand(1,2)[0], 0.0) * self.periodicity

            # Determine locations of chains already in the monolayer
            chemisorbed_chain_locations = np.array([chain.pos for chain
                                                    in self['chemisorbed_chain']])
            if 'crosslinked_chain' in self.labels:
                crosslinked_chain_locations = np.array([chain.pos for chain
                                                        in self['crosslinked_chain']])
            try:
                chain_locations = np.vstack((chemisorbed_chain_locations,
                                             crosslinked_chain_locations))
            except UnboundLocalError:
                chain_locations = chemisorbed_chain_locations
            chain_locations[:,2] = 0
            dists = [self.min_periodic_distance(site, loc) for loc in chain_locations]
            if min(dists) > spacing:
                failed_attempts = 0
                new_chain = mb.clone(chain)
                new_chain.translate_to(site)
                new_chain.translate([0.0, 0.0,
                    surface_level + (new_chain.pos[2]-min(new_chain.xyz[:,2]))+0.15])
                self.add(new_chain, 'crosslinked_chain[$]')
                new_chain_si = list(new_chain.particles_by_name('Si'))[0]
                self.crosslink_graph.add_node(new_chain_si, pos=(site[0], site[1]),
                    surface_bound=False)
                if verbose:
                    print('Added crosslinked chain {} ({} total chains)'
                          ''.format(len(self['crosslinked_chain']),
                          (len(self['crosslinked_chain']) +
                          len(self['chemisorbed_chain']))))
            else:
                failed_attempts += 1
                if verbose and failed_attempts % 100 == 0:
                    print('{} consecutive failed insertions (max {})'
                          ''.format(failed_attempts, int(max_failed_attempts)))

    def _determine_crosslink_network(self, verbose):
        """Determine crosslinks between monolayer chains.
        """
        pos = nx.get_node_attributes(self.crosslink_graph, 'pos')
        attachment = nx.get_node_attributes(self.crosslink_graph, 'surface_bound')
        nodes = self.crosslink_graph.nodes()
        for node in nodes:
            nodes_clone = nodes[:]
            if (not attachment[node] and
                    len(self.crosslink_graph.neighbors(node)) < 2):
                # Find the closest node that does not already have two edges
                attached = False
                while not attached:
                    attached_to = self._find_closest_node(node, nodes_clone, pos)
                    if len(self.crosslink_graph.neighbors(attached_to)) < 2:
                        self.crosslink_graph.add_edge(node, attached_to)
                        attached = True
                    else:
                        nodes_clone.remove(attached_to)

        '''
        Make sure that through the crosslink network, all chains can be traced
        to an attachment site.
          1. Loop over all subgraphs
          2. If subgraph cannot be traced to an attachment site, choose an end of
             the subgraph and attach a crosslink to the nearest node.
          3. Re-determine all subgraphs and repeat
          4. Note: Likely want some sort of distance criterion to make sure the
                   bonds being drawn aren't super crazy.
        '''
        all_connected = False
        while not all_connected:
            networks = list(nx.connected_components(self.crosslink_graph))
            if verbose:
                print('Found {} total crosslinking networks'.format(len(networks)))
            if all(any(attachment[site] for site in list(network))
                    for network in networks):
                all_connected = True
            else:
                if verbose:
                    print('Not all chains can be traced to a surface site, adding '
                          'additional crosslinks.')
                for network in networks:
                    network = list(network)
                    if not any([attachment[site] for site in network]):
                        nodes_clone = nodes[:]
                        n_neighbors = [len(self.crosslink_graph.neighbors(site))
                                       for site in network]
                        site = network[n_neighbors.index(1)]
                        attached = False
                        while not attached:
                            attached_to = self._find_closest_node(site, nodes_clone,
                                                                  pos)
                            if (len(self.crosslink_graph.neighbors(attached_to)) < 2
                                    and attached_to not in network):
                                self.crosslink_graph.add_edge(site, attached_to)
                                attached = True
                            else:
                                nodes_clone.remove(attached_to)
                        break

    def _create_crosslinks(self):
        """Creates the crosslinks between monolayer chains. """
        crosslinks = self.crosslink_graph.edges()
        pos = nx.get_node_attributes(self.crosslink_graph, 'pos')
        attachment = nx.get_node_attributes(self.crosslink_graph, 'surface_bound')
        for crosslink in crosslinks:
            atom1 = crosslink[0]
            atom2 = crosslink[1]
            loc1 = deepcopy(atom1.pos)
            loc2 = deepcopy(atom2.pos)
            for i in range(2):
                if loc1[i] - loc2[i] > self.periodicity[i]/2:
                    loc1[i] -= self.periodicity[i]
                elif loc1[i] - loc2[i] < -self.periodicity[i]/2:
                    loc1[i] += self.periodicity[i]
            # Add an oxygen atom midway between the two silicon atoms
            midpoint = (loc1 + loc2) / 2
            oxygen = mb.Particle(name='O', pos=midpoint)
            self.add(oxygen, 'crosslink_oxygen[$]')
            self.add_bond((atom1, oxygen))
            self.add_bond((atom2, oxygen))

        # Add hydroxyls to dangling bonds of silane silicons
        for node in self.crosslink_graph.nodes():
            neighbors = len(self.bond_graph.neighbors(node))
            if neighbors == 2:
                self._add_hydroxyls(node, 2)
            elif neighbors == 3:
                self._add_hydroxyls(node, 1)

    def _add_hydroxyls(self, node, n):
        if n == 2:
            hydroxyl = Hydroxyl()
            oxygen = list(hydroxyl.particles_by_name('O'))[0]
            hydroxyl.translate((node.pos - oxygen.pos) + np.array([0.15, 0, 0]))
            self.add(hydroxyl, 'hydroxyl[$]')
            self.add_bond((node, oxygen))
        hydroxyl = Hydroxyl()
        bond_vectors = [neighbor.pos - node.pos
                        for neighbor in self.bond_graph.neighbors(node)
                        if neighbor.name != 'C']
        crosslink_bond_vector = bond_vectors[np.argmin([np.dot(vector, [0,0,-1])
                                             for vector in bond_vectors])]
        for i in range(2):
            if abs(crosslink_bond_vector[i]) > self.periodicity[i]/2:
                crosslink_bond_vector[i] *= -1
        hydroxyl_bond_vector = ((crosslink_bond_vector * -1) /
                                np.linalg.norm(crosslink_bond_vector))
        oxygen = list(hydroxyl.particles_by_name('O'))[0]
        hydroxyl.translate((node.pos - oxygen.pos) + (hydroxyl_bond_vector)*0.15)
        hydrogen = list(hydroxyl.particles_by_name('H'))[0]
        hydrogen.translate_to(oxygen.pos + (hydroxyl_bond_vector) * 0.1)
        self.add(hydroxyl, 'hydroxyl[$]')
        self.add_bond((node, oxygen))

    def _find_closest_node(self, node, nodes, pos):
        dists = [self.min_periodic_distance(np.array([pos[node][0],pos[node][1],0.0]),
                 np.array([pos[a_node][0], pos[a_node][1], 0.0]))
                 if a_node != node else 1000.0 for a_node in nodes]
        return nodes[np.argmin(dists)]

    def _add_backfill(self, backfill, backfill_port_name):
        for port in self['surface'].available_ports():
            backfill_clone = mb.clone(backfill)
            mb.force_overlap(backfill_clone, backfill_clone[backfill_port_name],
                             port)
            self.add(backfill_clone)

    def draw_crosslink_network(self, filename):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt 
        pos = nx.get_node_attributes(self.crosslink_graph, 'pos')
        attachment = nx.get_node_attributes(self.crosslink_graph, 'surface_bound')
        surface_bound = [node for node, att in attachment.items()
                         if att]
        crosslinked = [node for node, att in attachment.items()
                       if not att]
        node_collection = nx.draw_networkx_nodes(self.crosslink_graph, pos,
            nodelist=surface_bound, node_color='black', linewidths=5.0,
            node_size=100)
        node_collection = nx.draw_networkx_nodes(self.crosslink_graph, pos,
            nodelist=surface_bound, node_color='red', node_size=100)
        node_collection = nx.draw_networkx_nodes(self.crosslink_graph, pos,
            nodelist=crosslinked, node_color='black', linewidths=5.0,
            node_size=100)
        node_collection = nx.draw_networkx_nodes(self.crosslink_graph, pos,
            nodelist=crosslinked, node_color='cyan', node_size=100)
        for edge in self.crosslink_graph.edges():
            point1 = list(pos[edge[0]])
            point2 = list(pos[edge[1]])
            fake_point1 = deepcopy(point1)
            fake_point2 = deepcopy(point2)
            for i in range(2):
                if point1[i] - point2[i] > self.periodicity[i]/2:
                    fake_point1[i] -= self.periodicity[i]
                    fake_point2[i] += self.periodicity[i]
                elif point1[i] - point2[i] < -self.periodicity[i]/2:
                    fake_point1[i] += self.periodicity[i]
                    fake_point2[i] -= self.periodicity[i]
            if fake_point1 == point1:
                plt.plot([point1[0], point2[0]], [point1[1], point2[1]],
                    marker='None', linestyle='-', color='black',
                    linewidth=1.0)
            else:
                plt.plot([point1[0], fake_point2[0]], [point1[1], fake_point2[1]],
                    marker='None', linestyle='-', color='black',
                    linewidth=1.0)
                plt.plot([point2[0], fake_point1[0]], [point2[1], fake_point1[1]],
                    marker='None', linestyle='-', color='black',
                    linewidth=1.0)
        plt.xlim(0.0, self.periodicity[0])
        plt.ylim(0.0, self.periodicity[1])
        plt.savefig(filename)

if __name__ == "__main__":
    from mbuild.examples import Alkane
    from mbuild.lib.atoms import H
    from mbuild.lib.bulk_materials import AmorphousSilica
    from mbuild.recipes import SilicaInterface

    seed = 12345
    surface = mb.SilicaInterface(bulk_silica=AmorphousSilica(),
        thickness=1.2, seed=seed)
    xlinked_monolayer = CrosslinkedMonolayer(Alkane(10, cap_end=False), surface,
        0.42, backfill=H(), max_failed_attempts=1e2, verbose=True,
        backfill_port_name='up')
    xlinked_monolayer.save('test.mol2')
    xlinked_monolayer.draw_crosslink_network('test.pdf')
